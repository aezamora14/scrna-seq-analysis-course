# Advanced scRNA-seq Analysis Course

This is the public student edition of a graduate-level computational course in single-cell RNA sequencing analysis. The course emphasizes how analytical decisions change measurements, visualizations, retained biology, and scientific conclusions.

The current release fully implements three foundational modules:

1. **What scRNA-seq count matrices actually measure**
2. **Quality control as a biological and statistical decision**
3. **Normalization and interpretation of expression values**

Later modules have intentionally not been scaffolded yet. Depth, executable exercises, and validation take priority over placeholder content.

## Learning approach

Each module follows:

> **Concept → computation → visualization → perturbation → interpretation → biological consequence**

Students repeatedly compare analyses that differ by one choice, then answer:

- What changed?
- What did not change?
- Why?
- Does the biological conclusion change?

## Repository structure

- `modules/` — fully developed student lessons for Modules 1–3
- `labs/` — runnable student lab scripts
- `data/` — deterministic, controlled teaching datasets
- `scripts/` — setup, data generation, figure generation, validation, and smoke tests
- `figures/` — reproducibly generated teaching figures
- `setup/` — local and Codespaces setup instructions
- `resources/` — glossary and course-wide reference material
- `references/` — citations and software documentation links
- `.devcontainer/` — Codespaces/devcontainer configuration
- `.github/workflows/` — automated validation and GitHub Pages publishing

This public repository contains no instructor notes or solutions. Those are maintained in a separate private instructor repository.

## Local setup

1. Install R 4.5 or later, RStudio, and Quarto.
2. Clone the repository:

   ```bash
   git clone https://github.com/aezamora14/scrna-seq-analysis-course.git
   cd scrna-seq-analysis-course
   ```

3. Restore the reproducible R environment:

   ```bash
   Rscript scripts/00_setup.R
   ```

4. Verify the environment and teaching data:

   ```bash
   Rscript scripts/check_environment.R
   Rscript scripts/generate_synthetic_data.R --check
   Rscript scripts/smoke_test_modules_01_03.R
   python3 scripts/validate_repository.py
   ```

5. Preview the course website:

   ```bash
   quarto preview
   ```

6. Open the corresponding script in `labs/` while reading each module.

## Data statement

The bundled dataset is synthetic and deliberately contains known cell populations, empty droplets, low-quality profiles, stressed cells, ambient RNA contamination, and simulated doublets. It is designed for controlled teaching—not benchmarking, clinical inference, or biological discovery.

## Reproducibility

- R dependencies are pinned in `renv.lock`.
- Synthetic data and figures are generated from version-controlled scripts with fixed seeds.
- Quarto rendering does not execute analysis code.
- Continuous integration checks the environment, data integrity, R workflows, internal links, and rendered site.

## License and disclaimer

Course code and synthetic data are for education. The synthetic data do not represent patients or clinical samples and must not be interpreted as biological evidence.
