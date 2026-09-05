required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "ggplot2",
  "patchwork",
  "dplyr",
  "tidyr",
  "tibble",
  "scales",
  "digest",
  "sctransform"
)

optional_packages <- c(
  "SingleCellExperiment",
  "scDblFinder"
)

optional_accelerators <- c(
  "glmGamPoi"
)

missing_required <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_required)) {
  stop(
    "Missing required R packages: ",
    paste(missing_required, collapse = ", "),
    ". Run Rscript scripts/00_setup.R."
  )
}

cat("R:", R.version.string, "\n")
for (package in required_packages) {
  cat(package, as.character(packageVersion(package)), "\n")
}

for (package in optional_packages) {
  status <- if (requireNamespace(package, quietly = TRUE)) {
    as.character(packageVersion(package))
  } else {
    "not installed (optional exercise will be skipped)"
  }
  cat(package, status, "\n")
}

for (package in optional_accelerators) {
  status <- if (requireNamespace(package, quietly = TRUE)) {
    paste(as.character(packageVersion(package)), "(optional accelerator)")
  } else {
    "not installed (SCTransform uses its slower built-in fitting path)"
  }
  cat(package, status, "\n")
}

quarto_candidates <- c(
  unname(Sys.which("quarto")),
  "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto",
  "C:/Program Files/Quarto/bin/quarto.exe"
)
quarto_candidates <- quarto_candidates[nzchar(quarto_candidates)]
quarto <- quarto_candidates[file.exists(quarto_candidates)][1]

if (is.na(quarto)) {
  warning("Quarto is not on PATH. R checks can run, but use an installed Quarto binary to render.")
} else {
  cat("Quarto:", system2(quarto, "--version", stdout = TRUE)[1], "\n")
}

message("Required environment check passed.")
