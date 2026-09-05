# AGENTS.md — Advanced scRNA-seq analysis course

## Scope

Maintain a graduate-level, semester-capable scRNA-seq analysis course. This repository is the **public student edition**. It must never contain instructor notes, solutions, answer keys, expected answers, private rubrics, grading guidance, or hidden annotations.

Only Modules 1–7 are implemented in the current release:

1. count matrices and measurement;
2. quality control;
3. normalization;
4. highly variable genes;
5. scaling and regression;
6. PCA;
7. choosing PCA dimensions.

Do not create placeholder modules or empty navigation entries for later topics.

## Pedagogical sequence

Each module should follow:

> **Concept → computation → visualization → perturbation → interpretation → biological consequence**

Every major concept should include a definition before jargon is used, a controlled example, runnable R code, a side-by-side comparison that changes one decision, interpretation questions, and limitations. Students should explain results in words rather than merely execute a pipeline.

## Scientific requirements

- Distinguish reads, UMIs, raw counts, normalized expression, and scaled expression.
- Explain biological zeros, sampling zeros, and structural/technical absence without claiming they can be perfectly identified from one observed zero.
- Treat QC thresholds as sample-, tissue-, protocol-, and question-dependent decisions.
- Distinguish stressed but viable cells from dying/low-quality cells.
- Explain that QC can preferentially remove biological populations and alter downstream conclusions.
- Distinguish normalization from scaling, batch correction, regression, imputation, and differential expression.
- Do not present SCTransform as automatically superior to log normalization.
- Treat donor/sample as the biological replicate for condition-level inference.
- Never imply that UMAP defines clusters or that a computational cluster is automatically a cell type.

## Code and data standards

- Use clear, current R and Seurat code; avoid unnecessary abstraction.
- Keep paths relative to the repository root.
- Keep controlled data small enough for an ordinary laptop.
- Use deterministic random seeds in simulations and stochastic analyses.
- Preserve raw integer UMI counts; never overwrite counts with normalized values.
- Keep required dependencies in `scripts/00_setup.R` and `renv.lock` synchronized.
- Optional dependencies must be labeled and tested when locally available.
- Do not commit large raw sequencing files, credentials, or tokens.

## Quarto standards

- Render without executing R by default (`execute: enabled: false`).
- Code shown on the site must be copy/paste runnable from the project root.
- Use `$...$` for inline math and `$$...$$` for display math.
- Preserve heading hierarchy, callouts, code-copy controls, and descriptive alt text.
- Do not add links to pages or files that do not exist.

## Student/instructor synchronization

Student-facing files are copied unchanged into the private instructor repository. Instructor-only material belongs only in the private repository's `instructor/`, `solutions/`, `rubrics/`, `teaching-notes/`, and `expected-results/` directories. Before publishing, compare shared files byte-for-byte and scan the public tree for instructor-only language and paths.

## Before finishing

1. Run `Rscript scripts/check_environment.R`.
2. Run `Rscript scripts/generate_synthetic_data.R --check`.
3. Run `Rscript scripts/smoke_test_modules_01_03.R`.
4. Run `python3 scripts/validate_repository.py`.
5. Render the complete site with `quarto render`.
6. Check internal links, images, code fences, and MathJax output.
7. Confirm public/instructor separation and shared-file synchronization.
8. Summarize modified files, tests, and untested assumptions.
