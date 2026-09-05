suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

source("scripts/course_helpers.R")
set.seed(20260904)

# Module 1 ------------------------------------------------------------------

obj <- add_qc_metrics(load_teaching_object())
counts <- LayerData(obj, assay = "RNA", layer = "counts")

stopifnot(
  inherits(counts, "sparseMatrix"),
  identical(dim(counts), c(260L, 700L)),
  all(Matrix::colSums(counts) == obj$nCount_RNA),
  all(Matrix::colSums(counts > 0) == obj$nFeature_RNA)
)

sparsity <- 1 - length(counts@x) / (nrow(counts) * ncol(counts))
stopifnot(sparsity > 0.4, sparsity < 1)

toy <- read.csv(
  "data/toy_counts.csv",
  row.names = 1,
  check.names = FALSE
)
toy <- as.matrix(toy)
stopifnot(
  colSums(toy)[["Cell_A"]] == 100,
  colSums(toy)[["Cell_B"]] == 200,
  all(toy[, "Cell_B"] == 2 * toy[, "Cell_A"])
)

set.seed(20260904)
candidate <- which(
  obj$truth_qc_class == "high_quality" &
    obj$nCount_RNA > median(obj$nCount_RNA)
)[1]
original <- as.integer(counts[, candidate])
undersampled <- rbinom(length(original), original, 0.35)
stopifnot(
  sum(undersampled) < sum(original),
  sum(undersampled > 0) <= sum(original > 0),
  any(original > 0 & undersampled == 0)
)

message("Module 1 smoke checks passed.")

# Module 2 ------------------------------------------------------------------

flags <- qc_strategy_flags(obj)
stopifnot(
  identical(names(flags), c("permissive", "moderate", "aggressive")),
  all(vapply(flags, length, integer(1)) == ncol(obj)),
  sum(flags$permissive) > sum(flags$moderate),
  sum(flags$moderate) > sum(flags$aggressive)
)

truth <- obj$truth_cell_type
aggressive_retention <- tapply(flags$aggressive, truth, mean)
stopifnot(
  aggressive_retention[["NK"]] < aggressive_retention[["Monocyte"]],
  aggressive_retention[["B"]] < aggressive_retention[["Monocyte"]],
  aggressive_retention[["Doublet"]] == 1
)

filtered_objects <- lapply(flags, function(keep) run_qc_pipeline(obj, keep))
stopifnot(
  all(vapply(filtered_objects, function(x) "pca" %in% Reductions(x), logical(1))),
  all(vapply(filtered_objects, function(x) "umap" %in% Reductions(x), logical(1))),
  all(vapply(filtered_objects, function(x) length(levels(Idents(x))) > 1, logical(1)))
)

message("Module 2 required smoke checks passed.")

if (
  requireNamespace("SingleCellExperiment", quietly = TRUE) &&
  requireNamespace("scDblFinder", quietly = TRUE)
) {
  suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(scDblFinder)
  })
  set.seed(20260904)
  sce <- as.SingleCellExperiment(obj)
  sce <- scDblFinder(
    sce,
    samples = "sample_id",
    verbose = FALSE,
    BPPARAM = BiocParallel::SerialParam(RNGseed = 20260904)
  )
  stopifnot(
    "scDblFinder.class" %in% colnames(colData(sce)),
    all(colData(sce)$scDblFinder.class %in% c("singlet", "doublet"))
  )
  message("Module 2 optional scDblFinder smoke check passed.")
} else {
  message("Module 2 optional scDblFinder check skipped: package not installed.")
}

# Module 3 ------------------------------------------------------------------

counts_per_100 <- sweep(toy, 2, colSums(toy), "/") * 100
manual_log <- log1p(counts_per_100)

stopifnot(
  counts_per_100["CD3D", "Cell_A"] == 20,
  counts_per_100["CD3D", "Cell_B"] == 20,
  isTRUE(all.equal(
    manual_log["CD3D", "Cell_A"],
    log(21),
    tolerance = 1e-12
  ))
)

toy_obj <- CreateSeuratObject(
  counts = toy,
  min.cells = 0,
  min.features = 0
)
toy_obj <- NormalizeData(
  toy_obj,
  normalization.method = "LogNormalize",
  scale.factor = 100,
  verbose = FALSE
)

stopifnot(
  all(as.matrix(LayerData(toy_obj, layer = "counts")) == toy),
  isTRUE(all.equal(
    as.matrix(LayerData(toy_obj, layer = "data")),
    manual_log,
    tolerance = 1e-10
  ))
)

obj_log <- NormalizeData(obj, scale.factor = 10000, verbose = FALSE)
stopifnot(
  "data" %in% Layers(obj_log[["RNA"]]),
  all(LayerData(obj_log, layer = "data") >= 0)
)

obj_sct <- SCTransform(
  obj,
  assay = "RNA",
  new.assay.name = "SCT",
  vst.flavor = "v2",
  return.only.var.genes = FALSE,
  verbose = FALSE
)
stopifnot(
  "SCT" %in% names(obj_sct@assays),
  "scale.data" %in% as.character(Layers(obj_sct[["SCT"]])),
  "NKG7" %in% rownames(LayerData(obj_sct, assay = "SCT", layer = "scale.data"))
)

message("Module 3 smoke checks passed.")

# Execute the complete student lab scripts in isolated environments --------

plot_device <- tempfile(fileext = ".pdf")
pdf(plot_device)
on.exit({
  invisible(dev.off())
  unlink(plot_device)
}, add = TRUE)

for (lab in c(
  "labs/01-count-matrix-lab.R",
  "labs/02-quality-control-lab.R",
  "labs/03-normalization-lab.R"
)) {
  sys.source(lab, envir = new.env(parent = globalenv()))
}

message("All student lab scripts executed successfully.")
message("Modules 1–3 smoke test passed.")
