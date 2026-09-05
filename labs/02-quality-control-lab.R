# Module 2 lab: quality-control sensitivity analysis
# Run from the repository root.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
})

source("scripts/course_helpers.R")
set.seed(20260904)
dir.create("outputs", showWarnings = FALSE)

obj <- add_qc_metrics(load_teaching_object())
qc_metadata <- obj[[]]

# Part 1: visualize QC distributions ---------------------------------------

qc_long <- qc_metadata |>
  select(
    sample_id,
    nCount_RNA,
    nFeature_RNA,
    percent.mt,
    percent.ribo,
    log10_genes_per_umi
  ) |>
  pivot_longer(
    -sample_id,
    names_to = "metric",
    values_to = "value"
  )

violin_plot <- ggplot(
  qc_long,
  aes(sample_id, value, fill = sample_id)
) +
  geom_violin(scale = "width", show.legend = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  labs(x = "Sample", y = "Observed value")

scatter_plot <- ggplot(
  qc_metadata,
  aes(nCount_RNA, nFeature_RNA, color = percent.mt)
) +
  geom_point(alpha = 0.65) +
  scale_x_log10() +
  scale_color_viridis_c() +
  labs(
    x = "Observed UMIs",
    y = "Detected genes",
    color = "% mitochondrial"
  )

print(violin_plot)
print(scatter_plot)

ggsave(
  "outputs/module-02-qc-violins.png",
  violin_plot,
  width = 10,
  height = 7,
  dpi = 150
)
ggsave(
  "outputs/module-02-qc-scatter.png",
  scatter_plot,
  width = 6,
  height = 5,
  dpi = 150
)

# Interpretation prompts:
# - Which sample has the largest typical library size?
# - Which sample would lose the most cells under a universal count minimum?
# - Which joint metric patterns are consistent with empty droplets or damage?

# Part 2: compare three threshold strategies -------------------------------

flags <- qc_strategy_flags(obj)

retention <- data.frame(
  cell_id = Cells(obj),
  sample_id = obj$sample_id,
  truth_cell_type = obj$truth_cell_type,
  truth_qc_class = obj$truth_qc_class,
  permissive = flags$permissive,
  moderate = flags$moderate,
  aggressive = flags$aggressive
) |>
  pivot_longer(
    cols = c(permissive, moderate, aggressive),
    names_to = "strategy",
    values_to = "retained"
  )

overall_retention <- retention |>
  group_by(strategy) |>
  summarise(
    removed = sum(!retained),
    retained = sum(retained),
    .groups = "drop"
  )

type_retention <- retention |>
  group_by(strategy, truth_cell_type) |>
  summarise(
    total = n(),
    retained_fraction = mean(retained),
    retained = sum(retained),
    .groups = "drop"
  )

qc_class_retention <- retention |>
  group_by(strategy, truth_qc_class) |>
  summarise(
    total = n(),
    retained_fraction = mean(retained),
    retained = sum(retained),
    .groups = "drop"
  )

print(overall_retention)
print(type_retention)
print(qc_class_retention)

retention_plot <- ggplot(
  type_retention,
  aes(truth_cell_type, retained_fraction, fill = strategy)
) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(x = "Controlled truth label", y = "Fraction retained") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(retention_plot)
ggsave(
  "outputs/module-02-retention.png",
  retention_plot,
  width = 9,
  height = 5,
  dpi = 150
)

# Interpretation prompts:
# - Which controlled biological populations are preferentially lost?
# - Which technical class survives all simple threshold strategies?
# - Why is "aggressive" not synonymous with "accurate"?

# Part 3: hold the pipeline fixed while changing retained cells -------------

filtered_objects <- lapply(
  flags,
  function(keep) run_qc_pipeline(obj, keep)
)

downstream_summary <- data.frame(
  strategy = names(filtered_objects),
  retained_cells = vapply(filtered_objects, ncol, numeric(1)),
  clusters = vapply(
    filtered_objects,
    function(filtered) length(levels(Idents(filtered))),
    integer(1)
  )
)
print(downstream_summary)

pca_plots <- Map(
  function(filtered, strategy) {
    DimPlot(
      filtered,
      reduction = "pca",
      group.by = "truth_cell_type",
      raster = FALSE
    ) +
      ggtitle(paste("PCA:", strategy))
  },
  filtered_objects,
  names(filtered_objects)
)

cluster_plots <- Map(
  function(filtered, strategy) {
    DimPlot(
      filtered,
      reduction = "umap",
      group.by = "seurat_clusters",
      label = TRUE,
      repel = TRUE,
      raster = FALSE
    ) +
      ggtitle(paste("Clusters:", strategy)) +
      NoLegend()
  },
  filtered_objects,
  names(filtered_objects)
)

pca_comparison <- wrap_plots(pca_plots, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")
cluster_comparison <- wrap_plots(cluster_plots, ncol = 3)

print(pca_comparison)
print(cluster_comparison)

ggsave(
  "outputs/module-02-pca-comparison.png",
  pca_comparison,
  width = 11,
  height = 5.5,
  dpi = 150
)
ggsave(
  "outputs/module-02-cluster-comparison.png",
  cluster_comparison,
  width = 11,
  height = 4.5,
  dpi = 150
)

# Part 4: optional doublet detection ---------------------------------------

if (
  requireNamespace("SingleCellExperiment", quietly = TRUE) &&
  requireNamespace("scDblFinder", quietly = TRUE)
) {
  suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(scDblFinder)
  })

  # Make the optional classifier as reproducible as practical. Its fitted
  # model can still change across package/platform versions, so interpret the
  # confusion table rather than expecting one exact set of barcode calls.
  set.seed(20260904)
  sce <- as.SingleCellExperiment(obj)
  sce <- scDblFinder(
    sce,
    samples = "sample_id",
    BPPARAM = BiocParallel::SerialParam(RNGseed = 20260904)
  )

  print(table(
    predicted = colData(sce)$scDblFinder.class,
    controlled_truth = obj$truth_qc_class == "doublet"
  ))
} else {
  message("Optional scDblFinder exercise skipped: package is not installed.")
}

# Final interpretation prompts:
# - What changed across QC strategies?
# - What did not change in the function parameters?
# - Why did the PCA and graph communities nevertheless change?
# - Which strategy would you defend for a stated biological question?
# - Which sensitivity results would you include in a report?

message("Module 2 student lab executed successfully.")
