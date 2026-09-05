# Module 3 lab: normalization
# Run from the repository root.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

source("scripts/course_helpers.R")
set.seed(20260904)
dir.create("outputs", showWarnings = FALSE)

# Part 1: hand-scale matrix -------------------------------------------------

toy <- read.csv(
  "data/toy_counts.csv",
  row.names = 1,
  check.names = FALSE
)
toy <- as.matrix(toy)

library_sizes <- colSums(toy)
counts_per_100 <- sweep(toy, 2, library_sizes, "/") * 100
log_normalized <- log1p(counts_per_100)

print(toy)
print(library_sizes)
print(round(counts_per_100, 3))
print(round(log_normalized, 3))

# Interpretation prompts:
# - Why do Cell_A and Cell_B become equal after library-size correction?
# - Which biological differences remain?
# - Why are the log-normalized values not UMI counts?

# Part 2: verify Seurat's implementation -----------------------------------

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

seurat_counts <- as.matrix(LayerData(toy_obj, layer = "counts"))
seurat_normalized <- as.matrix(LayerData(toy_obj, layer = "data"))

stopifnot(all(seurat_counts == toy))
stopifnot(isTRUE(all.equal(
  seurat_normalized,
  log_normalized,
  tolerance = 1e-10
)))

# Part 3: visualize raw, relative, and log-normalized ----------------------

matrix_to_long <- function(matrix, stage) {
  as.data.frame(matrix) |>
    tibble::rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "cell", values_to = "value") |>
    mutate(stage = stage)
}

comparison_long <- bind_rows(
  matrix_to_long(toy, "Raw UMI counts"),
  matrix_to_long(counts_per_100, "Counts per 100"),
  matrix_to_long(log_normalized, "log(1 + counts per 100)")
)

normalization_plot <- ggplot(
  comparison_long,
  aes(cell, gene, fill = value)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 1))) +
  facet_wrap(~stage, nrow = 1, scales = "free") +
  scale_fill_viridis_c() +
  labs(x = "Cell", y = "Gene")

print(normalization_plot)
ggsave(
  "outputs/module-03-hand-normalization.png",
  normalization_plot,
  width = 12,
  height = 5,
  dpi = 150
)

# Part 4: change only the scale factor -------------------------------------

scale_factors <- c(100, 1000, 10000)
normalized_by_scale <- lapply(
  scale_factors,
  function(scale_factor) {
    log1p(sweep(toy, 2, colSums(toy), "/") * scale_factor)
  }
)
names(normalized_by_scale) <- paste0("scale_", scale_factors)

print(lapply(
  normalized_by_scale,
  function(matrix) round(matrix["CD3D", ], 3)
))

# Interpretation prompts:
# - Which numerical values changed?
# - Which zeros or rankings changed?
# - Does the scale factor have a biological unit?

# Part 5: larger controlled dataset ----------------------------------------

obj <- add_qc_metrics(load_teaching_object())
obj_log <- NormalizeData(obj, scale.factor = 10000, verbose = FALSE)

raw_nkg7 <- as.numeric(LayerData(obj_log, layer = "counts")["NKG7", ])
log_nkg7 <- as.numeric(LayerData(obj_log, layer = "data")["NKG7", ])

depth_long <- data.frame(
  nCount_RNA = obj$nCount_RNA,
  truth_cell_type = obj$truth_cell_type,
  raw_NKG7 = raw_nkg7,
  normalized_NKG7 = log_nkg7
) |>
  filter(truth_cell_type %in% c("CD8 T", "NK")) |>
  pivot_longer(
    cols = c(raw_NKG7, normalized_NKG7),
    names_to = "expression_scale",
    values_to = "NKG7"
  )

depth_plot <- ggplot(
  depth_long,
  aes(nCount_RNA, NKG7, color = truth_cell_type)
) +
  geom_point(alpha = 0.55) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_log10() +
  facet_wrap(~expression_scale, scales = "free_y") +
  labs(
    x = "Observed library size",
    y = "NKG7 value",
    color = "Controlled cell type"
  )

print(depth_plot)
ggsave(
  "outputs/module-03-depth-comparison.png",
  depth_plot,
  width = 10,
  height = 5,
  dpi = 150
)

# Part 6: optional SCTransform comparison ----------------------------------

obj_sct <- SCTransform(
  obj,
  assay = "RNA",
  new.assay.name = "SCT",
  vst.flavor = "v2",
  return.only.var.genes = FALSE,
  verbose = FALSE
)

log_gzmk <- as.numeric(
  LayerData(obj_log, assay = "RNA", layer = "data")["GZMK", ]
)

sct_residual <- as.numeric(
  LayerData(
    obj_sct,
    assay = "SCT",
    layer = "scale.data"
  )["GZMK", ]
)

method_comparison <- data.frame(
  nCount_RNA = obj$nCount_RNA,
  truth_cell_type = obj$truth_cell_type,
  LogNormalize = log_gzmk,
  SCTransform_residual = sct_residual
) |>
  filter(truth_cell_type %in% c("CD8 T", "NK")) |>
  pivot_longer(
    cols = c(LogNormalize, SCTransform_residual),
    names_to = "method",
    values_to = "GZMK_value"
  )

sct_plot <- ggplot(
  method_comparison,
  aes(nCount_RNA, GZMK_value, color = truth_cell_type)
) +
  geom_point(alpha = 0.55) +
  scale_x_log10() +
  facet_wrap(~method, scales = "free_y") +
  labs(
    x = "Observed library size",
    y = "Method-specific GZMK value"
  )

print(sct_plot)
ggsave(
  "outputs/module-03-sctransform-comparison.png",
  sct_plot,
  width = 10,
  height = 5,
  dpi = 150
)

# Final interpretation prompts:
# - What changed between raw counts and LogNormalize?
# - What did not change?
# - Which problems remain after normalization?
# - Why are log-normalized values and Pearson residuals not directly comparable?
# - Does either method automatically produce a more valid biological result?

message("Module 3 student lab executed successfully.")
