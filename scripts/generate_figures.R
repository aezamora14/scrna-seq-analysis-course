suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
})

source("scripts/course_helpers.R")
set.seed(20260904)

dir.create("figures/module-01", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/module-02", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/module-03", recursive = TRUE, showWarnings = FALSE)

save_plot <- function(filename, plot, width, height) {
  ggsave(
    filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 160,
    bg = "white"
  )
}

# Module 1: sparse matrix and undersampling ---------------------------------

obj <- add_qc_metrics(load_teaching_object())
counts <- LayerData(obj, assay = "RNA", layer = "counts")

sparse_summary <- tibble(
  representation = c("Sparse dgCMatrix", "Dense matrix"),
  memory_bytes = c(
    as.numeric(object.size(counts)),
    as.numeric(object.size(as.matrix(counts)))
  )
)

memory_plot <- ggplot(
  sparse_summary,
  aes(representation, memory_bytes, fill = representation)
) +
  geom_col(width = 0.65, show.legend = FALSE) +
  scale_y_continuous(labels = label_bytes()) +
  labs(
    title = "Same counts, different storage",
    x = NULL,
    y = "Object size"
  ) +
  theme_minimal(base_size = 12)

candidate <- which(
  obj$truth_qc_class == "high_quality" &
    obj$nCount_RNA > median(obj$nCount_RNA)
)[1]
original <- as.integer(counts[, candidate])
undersampled <- rbinom(length(original), size = original, prob = 0.35)

undersampling_summary <- tibble(
  representation = rep(c("Original", "35% sampled"), each = 2),
  metric = rep(c("Library size", "Detected genes"), 2),
  value = c(
    sum(original),
    sum(original > 0),
    sum(undersampled),
    sum(undersampled > 0)
  )
)

undersampling_plot <- ggplot(
  undersampling_summary,
  aes(representation, value, fill = representation)
) +
  geom_col(width = 0.65, show.legend = FALSE) +
  facet_wrap(~metric, scales = "free_y") +
  labs(
    title = "Undersampling creates more observed zeros",
    x = NULL,
    y = "Observed value"
  ) +
  theme_minimal(base_size = 12)

save_plot(
  "figures/module-01/sparse-and-undersampling.png",
  memory_plot + undersampling_plot,
  width = 10,
  height = 4.2
)

selected_genes <- c(
  "CD3D", "TRAC", "NKG7", "GNLY", "MS4A1", "CD79A",
  "LST1", "LYZ", "S100A8", "MALAT1", "MT-CO1", "RPLP0"
)
selected_cells <- unlist(lapply(
  c("CD4 T", "CD8 T", "NK", "B", "Monocyte"),
  function(cell_type) head(Cells(obj)[obj$truth_cell_type == cell_type], 8)
))

matrix_long <- as.data.frame(as.matrix(counts[selected_genes, selected_cells])) |>
  tibble::rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "cell_id", values_to = "count") |>
  left_join(obj[[]] |> as_tibble(), by = "cell_id")

matrix_plot <- ggplot(
  matrix_long,
  aes(cell_id, gene, fill = log1p(count))
) +
  geom_tile() +
  facet_grid(~truth_cell_type, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(name = "log(1 + UMI)") +
  labs(
    title = "A sparse gene-by-cell matrix still contains biological structure",
    x = "Individual barcodes",
    y = "Gene"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

save_plot(
  "figures/module-01/sparse-matrix-structure.png",
  matrix_plot,
  width = 11,
  height = 5.5
)

# Module 2: QC distributions and threshold consequences --------------------

metadata <- obj[[]] |> as_tibble()

qc_long <- metadata |>
  select(cell_id, sample_id, nCount_RNA, nFeature_RNA, percent.mt) |>
  pivot_longer(
    cols = c(nCount_RNA, nFeature_RNA, percent.mt),
    names_to = "metric",
    values_to = "value"
  )

qc_violin <- ggplot(qc_long, aes(sample_id, value, fill = sample_id)) +
  geom_violin(scale = "width", trim = TRUE, show.legend = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, show.legend = FALSE) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  labs(
    title = "QC distributions differ across samples",
    x = "Sample",
    y = "Observed value"
  ) +
  theme_minimal(base_size = 11)

save_plot(
  "figures/module-02/qc-distributions-by-sample.png",
  qc_violin,
  width = 11,
  height = 4.2
)

flags <- qc_strategy_flags(obj)
retention <- tibble(
  cell_id = Cells(obj),
  truth_cell_type = obj$truth_cell_type,
  permissive = flags$permissive,
  moderate = flags$moderate,
  aggressive = flags$aggressive
) |>
  pivot_longer(
    cols = c(permissive, moderate, aggressive),
    names_to = "strategy",
    values_to = "retained"
  ) |>
  group_by(strategy, truth_cell_type) |>
  summarise(retained_fraction = mean(retained), .groups = "drop") |>
  mutate(strategy = factor(strategy, levels = c("permissive", "moderate", "aggressive")))

retention_plot <- ggplot(
  retention,
  aes(truth_cell_type, retained_fraction, fill = strategy)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(
    title = "Threshold choices remove populations unequally",
    x = "Controlled truth label",
    y = "Barcodes retained",
    fill = "QC strategy"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_plot(
  "figures/module-02/qc-retention-by-cell-type.png",
  retention_plot,
  width = 9,
  height = 5.2
)

pipeline_objects <- lapply(flags, function(keep) run_qc_pipeline(obj, keep))

truth_colours <- c(
  "B" = "#F8766D",
  "CD4 T" = "#C49A00",
  "CD8 T" = "#53B400",
  "Doublet" = "#00C094",
  "Empty droplet" = "#00B6EB",
  "Monocyte" = "#619CFF",
  "NK" = "#F564E3"
)

truth_plots <- Map(
  function(filtered, strategy) {
    DimPlot(
      filtered,
      reduction = "pca",
      group.by = "truth_cell_type",
      raster = FALSE
    ) +
      scale_colour_manual(
        values = truth_colours,
        limits = names(truth_colours),
        drop = FALSE
      ) +
      ggtitle(paste("PCA:", strategy))
  },
  pipeline_objects,
  names(pipeline_objects)
)

truth_plots[2:3] <- lapply(
  truth_plots[2:3],
  function(plot) plot + guides(colour = "none")
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
      guides(colour = "none", fill = "none")
  },
  pipeline_objects,
  names(pipeline_objects)
)

save_plot(
  "figures/module-02/qc-downstream-comparison.png",
  wrap_plots(c(truth_plots, cluster_plots), ncol = 3, guides = "collect") &
    theme(legend.position = "bottom"),
  width = 12,
  height = 8
)

# Module 3: hand normalization and depth dependence -------------------------

toy <- read.csv("data/toy_counts.csv", row.names = 1, check.names = FALSE)
toy <- as.matrix(toy)
counts_per_100 <- sweep(toy, 2, colSums(toy), "/") * 100
log_normalized <- log1p(counts_per_100)

matrix_to_long <- function(x, stage) {
  as.data.frame(x) |>
    tibble::rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "cell", values_to = "value") |>
    mutate(stage = stage)
}

normalization_long <- bind_rows(
  matrix_to_long(toy, "Raw UMI counts"),
  matrix_to_long(counts_per_100, "Counts per 100"),
  matrix_to_long(log_normalized, "log(1 + counts per 100)")
) |>
  mutate(stage = factor(
    stage,
    levels = c("Raw UMI counts", "Counts per 100", "log(1 + counts per 100)")
  ))

normalization_plot <- ggplot(normalization_long, aes(cell, gene, fill = value)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = round(value, 1)), size = 3) +
  facet_wrap(~stage, nrow = 1, scales = "free") +
  scale_fill_viridis_c() +
  labs(
    title = "Raw counts, library-size correction, and log transformation",
    x = "Cell",
    y = "Gene",
    fill = "Value"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(
  "figures/module-03/normalization-hand-example.png",
  normalization_plot,
  width = 12,
  height = 5
)

obj_log <- NormalizeData(obj, scale.factor = 10000, verbose = FALSE)
nkg7_raw <- as.numeric(counts["NKG7", ])
nkg7_log <- as.numeric(LayerData(obj_log, assay = "RNA", layer = "data")["NKG7", ])

depth_data <- tibble(
  nCount_RNA = obj$nCount_RNA,
  truth_cell_type = obj$truth_cell_type,
  raw = nkg7_raw,
  log_normalized = nkg7_log
) |>
  filter(truth_cell_type %in% c("CD8 T", "NK")) |>
  pivot_longer(
    cols = c(raw, log_normalized),
    names_to = "stage",
    values_to = "NKG7"
  ) |>
  mutate(stage = factor(stage, levels = c("raw", "log_normalized")))

depth_plot <- ggplot(
  depth_data,
  aes(nCount_RNA, NKG7, color = truth_cell_type)
) +
  geom_point(alpha = 0.55, size = 1.3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~stage, scales = "free_y") +
  scale_x_log10() +
  labs(
    title = "Normalization reduces, but does not erase, depth relationships",
    x = "Observed library size (log scale)",
    y = "NKG7 value",
    color = "Controlled cell type"
  ) +
  theme_minimal(base_size = 11)

save_plot(
  "figures/module-03/depth-before-after-normalization.png",
  depth_plot,
  width = 10,
  height = 4.8
)

message("Teaching figures generated successfully.")
