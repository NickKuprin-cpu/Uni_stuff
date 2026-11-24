# Working in this repository

Use this workspace primarily to help the user study physics and to support reproducible research. Work in the existing checkout. Follow the [repository layout](README.md#where-work-goes) and create category folders as needed.

## Teach and verify physics

- Match the help to the request. For guided study, build on the user's attempt and give a focused hint or question. For a requested solution, show the reasoning and enough intermediate steps to reproduce it. For exam prep, offer retrieval prompts or timed problems when relevant.
- Keep the problem statement and the user's attempt distinct from any correction or solution you add. Use the course's notation where practical and identify the source handout and problem number when working from one.
- State assumptions and define symbols when needed. Check units, signs, and relevant limits or symmetries. For numerical work, compare with a known case or an independent calculation when feasible. Mark unverified results as such.
- For PDFs, use `pdftotext -layout "<file>" -` for a first pass. Inspect the rendered page when equations, subscripts, or diagrams matter. Cite the file and page or problem number so the work is easy to find again.

## Write LaTeX documents

- Keep documents short and readable. Use the required template or a standard document class, a small preamble, and only the headings, equations, and figures that help the reader follow the physics. Keep enough derivation to justify the result.
- For edits, preserve the document's structure and style. Keep changes focused on the requested content and the formatting needed for it. Add packages or macros only when they solve a concrete problem.
- Prefer pdfLaTeX when it supports the document's template, fonts, and language. If another engine is required, record the reason and build command alongside the source. For a standard pdfLaTeX build, run `latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex` from the document folder, using the actual entry filename. Inspect the PDF for layout issues and check the build for unresolved references and citations.

## Make work easy to resume

- Save persistent work with its subject and purpose, or in its research project folder. Keep original files intact and link to them with relative paths. Use Markdown with LaTeX for written derivations where practical.
- Record the next concrete step in the relevant note or project README when a session leaves work open. For research and simulations, record sources, assumptions, parameters, and how to reproduce the result. Distinguish hypotheses from findings.
- Document dependencies and run commands in the project that needs them. Use its existing tools where practical and verify the workflow with a representative calculation or run. Report which checks actually ran.
