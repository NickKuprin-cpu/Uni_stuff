# Physics study and research

This workspace is primarily for learning physics through problem solving, derivations, retrieval practice, and computational experiments. It also holds side projects and research that can be picked up and reproduced later.

## Start here

[Quantum Mechanics 1](Quantum_Mechanics_1/) contains the [homework sheets](Quantum_Mechanics_1/homework/) and the [original course description](Quantum_Mechanics_1/Description). The sheets retain their original filenames. Quote their paths in shell commands to handle spaces and parentheses.

## Where work goes

Give each subject or course a stable top-level folder, using `Quantum_Mechanics_1` as the naming pattern. Add folders as work arrives. Inside a subject, use these names consistently.

| Folder | Contents |
| --- | --- |
| `notes/` | Lecture notes, topic summaries, and derivations. |
| `homework/` | Assigned problem sheets and corresponding work. |
| `controls/` | Quizzes, in-term tests, and corrections. |
| `exams/` | Past papers, mock exams, and exam reviews. |
| `practice/` | Self-selected problems, retrieval prompts, and worked examples. |
| `labs/` | Lab instructions, measurements, and reports. |
| `projects/` | Side projects tied to the subject, one folder per project. |
| `references/` | Subject-wide references and reading notes. |

Use top-level `research/<project>/` for research that spans subjects or has its own question and results. Link back to the relevant course notes or projects from the research README. Keep sources used only by one project with that project.

Keep each item in one place and link to it from other contexts. Preserve original handouts and raw data. Put your own attempts, solutions, and processed results in separate files. Name new files for what they contain, such as `hw-03-attempt.md` or `harmonic-oscillator.md`. When repeated terms make names ambiguous, use year or term subfolders inside the relevant category, such as `exams/2026/`.

## Pick up where you left off

- For homework and practice, link to the original problem and record your attempt and where you got stuck. After resolving it, write a short correction and a prompt you can answer from memory later.
- For controls and exams, record the missed concept, why the first approach failed, and a follow-up problem that tests the same idea.
- When a subject has several active threads, add a short subject `README.md` with links to active work, open questions, and the next concrete step. Update it when the next step changes.
- For a project, keep a short `README.md` with the question, input and source locations, dependencies, the command or notebook steps to reproduce the result, and any checks or limitations. Use relative paths in commands and links. Keep large reproducible outputs out of Git and record how to regenerate them.

Keep LaTeX source and its figures with the relevant coursework or project. Prefer pdfLaTeX and use the [LaTeX writing and build guide](AGENTS.md#write-latex-documents) for short, readable documents.

The [agent instructions](AGENTS.md) also cover tutoring, checking physics, and keeping the workspace easy to resume.
