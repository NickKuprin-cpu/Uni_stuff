"""Построение осциллограмм и расчёт величин для отчёта по ВОЛС."""

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "raw data"
FIGURES = ROOT / "figures"
C_M_S = 299_792_458
GROUP_INDEX = 1.50
REFERENCE_GROUP_INDEX = 1.4682
PEAK_SPLIT_NS = 1_000.0
BASELINE_END_NS = -0.60
WINDOW_NS = (-1.5, 8.0)


@dataclass(frozen=True)
class Trace:
    time_ns: np.ndarray
    signal: np.ndarray

    @property
    def step_ns(self) -> float:
        return float(np.median(np.diff(self.time_ns)))


@dataclass(frozen=True)
class Pulse:
    trace: Trace
    background: float
    background_std: float
    background_count: int
    front_ns: float


def read_trace(filename: str) -> Trace:
    with (DATA / filename).open(encoding="ascii") as source:
        for line in source:
            if line.strip() == "Time,Ampl":
                break
        else:
            raise ValueError(f"В {filename} нет заголовка Time,Ampl")
        samples = np.loadtxt(source, delimiter=",", ndmin=2)

    if samples.shape[0] < 3 or samples.shape[1] != 2 or not np.isfinite(samples).all():
        raise ValueError(f"Некорректные отсчёты в {filename}")
    time_s, signal = samples.T
    if np.any(np.diff(time_s) <= 0):
        raise ValueError(f"Время в {filename} должно возрастать")
    return Trace(time_s * 1e9, signal)


def prepare_pulse(
    trace: Trace, cutoff_ns: float = BASELINE_END_NS, *,
    median: bool = False, background_shift: float = 0.0,
) -> Pulse:
    before_front = trace.signal[trace.time_ns <= cutoff_ns]
    if before_front.size < 2:
        raise ValueError("Недостаточно отсчётов до фронта")
    background = float(np.median(before_front) if median else np.mean(before_front))
    background += background_shift
    background_std = float(np.std(before_front, ddof=1))
    signal = trace.signal - background
    peak_index = int(np.argmax(signal))
    half = signal[peak_index] / 2
    crossings = np.flatnonzero(
        (signal[:peak_index] <= half) & (signal[1 : peak_index + 1] > half)
    )
    if half <= 0 or crossings.size == 0:
        raise ValueError("В записи отсутствует передний фронт импульса")
    index = crossings[-1]
    front = float(np.interp(half, signal[index : index + 2], trace.time_ns[index : index + 2]))
    return Pulse(
        Trace(trace.time_ns - front, signal), background,
        background_std, before_front.size, front,
    )


def pulse_area(
    pulse: Pulse, window_ns: tuple[float, float], *, coarse_phase: int | None = None
) -> float:
    start, stop = window_ns
    t, signal = pulse.trace.time_ns, pulse.trace.signal
    if not t[0] <= start < stop <= t[-1]:
        raise ValueError(f"Окно {window_ns} нс выходит за пределы записи")
    if coarse_phase is not None:
        # Крайние точки сохраняют общий диапазон интерполяции.
        retained = np.unique(np.r_[0, np.arange(coarse_phase, len(t), 2), len(t) - 1])
        t, signal = t[retained], signal[retained]
    window_times = np.r_[start, t[(t > start) & (t < stop)], stop]
    values = np.interp(window_times, t, signal)
    return float(np.sum(np.diff(window_times) * (values[:-1] + values[1:]) / 2))


def ru(value: float, digits: int) -> str:
    return f"{value:.{digits}f}".replace(".", ",")


def processing_loss_range(
    raw_pulses: tuple[Trace, Trace], nominal_pulses: tuple[Pulse, Pulse]
) -> tuple[float, float]:
    variants = [
        (WINDOW_NS[1], BASELINE_END_NS, False),
        (5.0, BASELINE_END_NS, False),
        (6.0, BASELINE_END_NS, False),
        (WINDOW_NS[1], BASELINE_END_NS, True),
        (WINDOW_NS[1], -1.0, False),
    ]
    losses = []
    for stop, cutoff, median in variants:
        area_ranges = []
        for trace, nominal in zip(raw_pulses, nominal_pulses):
            areas = []
            noise = nominal.background_std
            for shift in (-noise, 0.0, noise):
                pulse = prepare_pulse(trace, cutoff, median=median, background_shift=shift)
                areas.extend(
                    pulse_area(pulse, (WINDOW_NS[0], stop), coarse_phase=phase)
                    for phase in (None, 0, 1)
                )
            if min(areas) <= 0:
                raise ValueError("Площади импульсов должны быть положительными")
            area_ranges.append((min(areas), max(areas)))
        (short_min, short_max), (long_min, long_max) = area_ranges
        losses.extend((10 * np.log10(short_min / long_max), 10 * np.log10(short_max / long_min)))
    return min(losses), max(losses)


def write_results(
    overview: Trace, peaks_ns: tuple[float, float], raw_pulses: tuple[Trace, Trace],
    pulses: tuple[Pulse, Pulse],
) -> None:
    delay_ns = peaks_ns[1] - peaks_ns[0]
    length_m = C_M_S * delay_ns * 1e-9 / GROUP_INDEX
    delay_grid_ns = 2 * overview.step_ns
    length_grid_m = C_M_S * delay_grid_ns * 1e-9 / GROUP_INDEX
    reference_length_m = C_M_S * delay_ns * 1e-9 / REFERENCE_GROUP_INDEX
    areas = tuple(pulse_area(pulse, WINDOW_NS) for pulse in pulses)
    if min(areas) <= 0:
        raise ValueError("Площади импульсов должны быть положительными")
    ratio = areas[1] / areas[0]
    loss = -10 * np.log10(ratio)
    relative_loss = loss / (length_m / 1_000)
    loss_min, loss_max = processing_loss_range(raw_pulses, pulses)
    relative_min = loss_min / ((length_m + length_grid_m) / 1_000)
    relative_max = loss_max / ((length_m - length_grid_m) / 1_000)

    values = {
        "OverviewStep": ru(overview.step_ns, 0),
        "FineStep": ru(raw_pulses[0].step_ns, 1),
        "PeakDelay": ru(delay_ns / 1_000, 2),
        "DelayGrid": ru(delay_grid_ns, 0),
        "GroupIndex": ru(GROUP_INDEX, 2),
        "LengthMeters": ru(length_m, 0),
        "LengthGrid": ru(length_grid_m, 0),
        "ReferenceIndex": ru(REFERENCE_GROUP_INDEX, 4),
        "ReferenceLength": ru(reference_length_m, 0),
        "WindowStart": ru(WINDOW_NS[0], 1),
        "WindowEnd": ru(WINDOW_NS[1], 1),
        "ShortArea": ru(areas[0], 3),
        "LongArea": ru(areas[1], 3),
        "ShortNoise": ru(pulses[0].background_std, 4),
        "LongNoise": ru(pulses[1].background_std, 4),
        "PathLoss": ru(loss, 2),
        "LossLower": ru(loss_min, 2),
        "LossUpper": ru(loss_max, 2),
        "EffectiveLoss": ru(relative_loss, 2),
        "EffectiveLower": ru(relative_min, 2),
        "EffectiveUpper": ru(relative_max, 2),
    }
    lines = ["% Численные результаты из analysis.py"]
    for name, value in values.items():
        number = value.replace(",", "{,}")
        lines.append(rf"\newcommand{{\{name}}}{{{number}}}")
    (ROOT / "results.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Задержка {delay_ns / 1_000:.2f} мкс, длина {length_m:.2f} м, вклад сетки ±{length_grid_m:.2f} м")
    for label, pulse in zip(("Короткое плечо", "Плечо с катушкой"), pulses):
        print(
            f"{label}: фон {pulse.background:.6f} по {pulse.background_count} отсчётам, "
            f"СКО {pulse.background_std:.6f}, полувысота переднего фронта {pulse.front_ns:.3f} нс"
        )
    print(f"S0={areas[0]:.6f}, SL={areas[1]:.6f}, SL/S0={ratio:.6f}")
    print(f"Ослабление {loss:.4f} дБ, эффективное удельное значение {relative_loss:.4f} дБ/км")
    print(f"Диапазон ослабления по обработке: {loss_min:.4f}–{loss_max:.4f} дБ")
    print(f"Диапазон удельного ослабления: {relative_min:.4f}–{relative_max:.4f} дБ/км")


def finish_figure(fig: plt.Figure, ax: plt.Axes, filename: str) -> None:
    decimal_comma = FuncFormatter(lambda value, _: f"{value:g}".replace(".", ","))
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(MaxNLocator(nbins=7))
        axis.set_major_formatter(decimal_comma)
    ax.grid(alpha=0.23, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("Сигнал над фоном, ед. Ampl")
    fig.savefig(FIGURES / filename, metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


def plot_traces(
    overview: Trace, peak_indices: tuple[int, int], pulses: tuple[Pulse, Pulse]
) -> None:
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 10, "pdf.fonttype": 42})

    time_us = overview.time_ns / 1_000
    overview_background = np.mean(overview.signal[overview.time_ns < -100])
    signal = overview.signal - overview_background
    fig, ax = plt.subplots(figsize=(6.2, 2.3), layout="constrained")
    ax.plot(time_us, signal, color="#2b5c84", linewidth=0.9)
    for index, offset in zip(peak_indices, ((9, 3), (-9, 8))):
        t_peak = time_us[index]
        ax.scatter(t_peak, signal[index], color="#173c5c", s=17, zorder=3)
        ax.annotate(
            f"{ru(t_peak, 2)} мкс", (t_peak, signal[index]), xytext=offset,
            textcoords="offset points", va="bottom",
            ha="left" if offset[0] > 0 else "right",
        )
    ax.set_ylim(-0.025, 0.51)
    ax.set_xlabel("Время общей записи, мкс")
    finish_figure(fig, ax, "delay.pdf")

    fig, ax = plt.subplots(figsize=(6.2, 2.9), layout="constrained")
    ax.axvspan(*WINDOW_NS, color="#d8e3ed", alpha=0.55, linewidth=0)
    ax.axvline(0, color="#687582", linewidth=0.8, linestyle="--")
    for pulse, color, label in zip(
        pulses, ("#225e87", "#b24827"), ("Короткое плечо", "Плечо с катушкой")
    ):
        ax.plot(
            pulse.trace.time_ns, pulse.trace.signal, color=color, label=label,
            linewidth=1.1, marker=".", markersize=2.4,
        )
    ax.set_xlim(-2.2, 9.0)
    ax.set_ylim(-0.06, 1.24)
    ax.set_xlabel("Время от полувысоты переднего фронта, нс")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    finish_figure(fig, ax, "pulses.pdf")


def main() -> None:
    overview = read_trace("C2300000.txt")
    split = int(np.searchsorted(overview.time_ns, PEAK_SPLIT_NS))
    peak_indices = (
        int(np.argmax(overview.signal[:split])),
        split + int(np.argmax(overview.signal[split:])),
    )
    peaks_ns = tuple(float(overview.time_ns[index]) for index in peak_indices)
    raw_pulses = (read_trace("C2300001.txt"), read_trace("C2300002.txt"))
    pulses = tuple(prepare_pulse(trace) for trace in raw_pulses)
    write_results(overview, peaks_ns, raw_pulses, pulses)
    plot_traces(overview, peak_indices, pulses)


if __name__ == "__main__":
    main()
