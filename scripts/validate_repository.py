from __future__ import annotations

import csv
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

REQUIRED_FILES = [
    "AGENTS.md",
    "README.md",
    "COURSE_EDITION",
    "_quarto.yml",
    "index.qmd",
    "schedule-quarter.qmd",
    "schedule-semester.qmd",
    "renv.lock",
    ".devcontainer/devcontainer.json",
    ".github/workflows/validate.yml",
    ".github/workflows/publish.yml",
    "setup/index.qmd",
    "modules/01-count-matrix/index.qmd",
    "modules/02-quality-control/index.qmd",
    "modules/03-normalization/index.qmd",
    "labs/01-count-matrix-lab.R",
    "labs/02-quality-control-lab.R",
    "labs/03-normalization-lab.R",
    "scripts/generate_synthetic_data.R",
    "scripts/generate_figures.R",
    "scripts/smoke_test_modules_01_03.R",
    "data/synthetic_qc_counts.rds",
    "data/synthetic_qc_metadata.csv",
    "data/synthetic_data_manifest.csv",
]

REQUIRED_FILES += ['modules/04-variable-features/index.qmd', 'modules/05-scaling-regression/index.qmd', 'modules/06-pca/index.qmd', 'modules/07-choosing-dimensions/index.qmd', 'labs/04-variable-features-lab.R', 'labs/05-scaling-regression-lab.R', 'labs/06-pca-lab.R', 'labs/07-choosing-dimensions-lab.R', 'scripts/smoke_test_modules_04_07.R', 'scripts/generate_representation_figures.R']

PUBLIC_FORBIDDEN_PATHS = [
    "instructor",
    "solutions",
    "rubrics",
    "teaching-notes",
    "expected-results",
]

missing = [path for path in REQUIRED_FILES if not (ROOT / path).exists()]
if missing:
    raise SystemExit(f"Missing required files: {missing}")

module_dirs = sorted(
    path.name for path in (ROOT / "modules").iterdir() if path.is_dir()
)
expected_modules = [
    "01-count-matrix",
    "02-quality-control",
    "03-normalization",
    "04-variable-features",
    "05-scaling-regression",
    "06-pca",
    "07-choosing-dimensions",

]
if module_dirs != expected_modules:
    raise SystemExit(
        f"Expected only Modules 1–7; found module directories: {module_dirs}"
    )

edition = (ROOT / "COURSE_EDITION").read_text(encoding="utf-8").strip()
if edition not in {"student", "instructor"}:
    raise SystemExit(f"Unknown COURSE_EDITION value: {edition!r}")

if edition == "student":
    for path in PUBLIC_FORBIDDEN_PATHS:
        if (ROOT / path).exists():
            raise SystemExit(f"Instructor-only path exists in public repository: {path}")

student_sources = list((ROOT / "modules").rglob("*.qmd")) + list(
    (ROOT / "labs").rglob("*.R")
)
forbidden_phrases = re.compile(
    r"instructor note|answer key|expected answer|grading key",
    flags=re.IGNORECASE,
)
for path in student_sources:
    text = path.read_text(encoding="utf-8")
    match = forbidden_phrases.search(text)
    if match:
        raise SystemExit(
            f"Possible instructor-only phrase in {path.relative_to(ROOT)}: "
            f"{match.group(0)!r}"
        )

qmd_files = list(ROOT.rglob("*.qmd"))
for path in qmd_files:
    text = path.read_text(encoding="utf-8")
    relative = path.relative_to(ROOT)
    if not text.startswith("---\n"):
        raise SystemExit(f"Missing YAML front matter in {relative}")
    if text.count("```") % 2:
        raise SystemExit(f"Unbalanced fenced code blocks in {relative}")
    if text.count("$$") % 2:
        raise SystemExit(f"Unbalanced display-math delimiters in {relative}")
    if "\\[" in text or "\\]" in text:
        raise SystemExit(f"Unsupported display-math delimiter in {relative}")
    if text.count("\n:::") % 2:
        raise SystemExit(f"Unbalanced Quarto div/callout fences in {relative}")

link_pattern = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
for path in qmd_files:
    text = path.read_text(encoding="utf-8")
    for raw_target in link_pattern.findall(text):
        target = raw_target.split()[0].strip("<>")
        if target.startswith(("http://", "https://", "mailto:", "#")):
            continue
        target_path = target.split("#", 1)[0]
        if not target_path:
            continue
        resolved = (path.parent / target_path).resolve()
        if not resolved.exists():
            raise SystemExit(
                f"Broken local link in {path.relative_to(ROOT)}: {target}"
            )

manifest_path = ROOT / "data/synthetic_data_manifest.csv"
with manifest_path.open(newline="", encoding="utf-8") as handle:
    manifest = list(csv.DictReader(handle))

for row in manifest:
    path = ROOT / row["file"]
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != row["sha256"]:
        raise SystemExit(f"Checksum mismatch for {path.relative_to(ROOT)}")
    if path.stat().st_size != int(float(row["bytes"])):
        raise SystemExit(f"File-size mismatch for {path.relative_to(ROOT)}")

ignored_generated_dirs = {".git", ".quarto", "_site", "outputs", "library"}
large_files = [
    (path.relative_to(ROOT), path.stat().st_size)
    for path in ROOT.rglob("*")
    if path.is_file()
    and not ignored_generated_dirs.intersection(path.relative_to(ROOT).parts)
    and path.stat().st_size > 20 * 1024 * 1024
]
if large_files:
    raise SystemExit(f"Files larger than 20 MiB found: {large_files}")

print(
    "Repository validation passed: only Modules 1–7 are present; "
    f"{len(qmd_files)} Quarto sources, data checksums, links, math, and "
    f"{edition}-edition content boundaries are valid."
)
