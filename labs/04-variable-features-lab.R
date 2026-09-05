# Module 4 lab: highly variable genes
# Run from the repository root.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

source("scripts/course_helpers.R")
set.seed(20260907)
dir.create("outputs", showWarnings = FALSE)

# Part 1: controlled mean and variance -------------------------------------

toy_hvg <- rbind(
  stable_high_housekeeping = c(8.0, 8.1, 7.9, 8.0, 8.1, 7.9, 8.0, 8.0),
  lineage_marker = c(0.0, 0.1, 0.0, 0.2, 5.0, 5.2, 4.9, 5.1),
  low_noisy = c(0, 0, 0, 1, 0, 0, 0, 0),
  genuinely_variable = c(0.5, 1.0, 1.5, 2.0, 3.5, 4.0, 4.5, 5.0),
  important_but_stable = c(1.0, 1.0, 1.1, 0.9, 1.0, 1.0, 0.9, 1.1)
)
colnames(toy_hvg) <- paste0("Cell_", 1:8)

toy_summary <- data.frame(
  gene_pattern = rownames(toy_hvg),
  mean = rowMeans(toy_hvg),
  variance = apply(toy_hvg, 1, var),
  row.names = NULL
)
print(toy_summary)

x <- toy_hvg["lineage_marker", ]
manual_mean <- sum(x) / length(x)
manual_variance <- sum((x - manual_mean)^2) / (length(x) - 1)
stopifnot(
  isTRUE(all.equal(manual_mean, mean(x))),
  isTRUE(all.equal(manual_variance, var(x)))
)

toy_plot <- ggplot(toy_summary, aes(mean, variance, label = gene_pattern)) +
  geom_point(size = 3) +
  geom_text(nudge_y = 0.25, check_overlap = TRUE) +
  labs(x = "Mean normalized expression", y = "Variance across cells")
print(toy_plot)

# Interpretation prompts:
# - Which pattern separates the first four from the last four cells?
# - Why is high expression not equivalent to high variability?
# - Why might one isolated low observation be an unstable signal?
# - Can a stable gene still be biologically essential?

# Part 2: empirical mean--variance relationship ---------------------------

obj <- load_representation_object()
obj <- NormalizeData(obj, normalization.method = "LogNormalize", verbose = FALSE)
normalized <- as.matrix(LayerData(obj, assay = "RNA", layer = "data"))

gene_summary <- data.frame(
  gene = rownames(normalized),
  mean = rowMeans(normalized),
  variance = apply(normalized, 1, var),
  row.names = NULL
)

mean_variance_plot <- ggplot(gene_summary, aes(mean, variance)) +
  geom_point(alpha = 0.25, size = 0.7) +
  scale_x_continuous(trans = "log1p") +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = "Mean log-normalized expression (log1p axis)",
    y = "Variance across cells (log1p axis)"
  )

print(mean_variance_plot)
ggsave(
  "outputs/module-04-mean-variance.png",
  mean_variance_plot,
  width = 7,
  height = 5,
  dpi = 150
)

# Part 3: Seurat variable features ----------------------------------------

obj_2000 <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)
hvg_2000 <- VariableFeatures(obj_2000)

stopifnot(length(hvg_2000) == 2000L)
print(head(hvg_2000, 20))

markers <- c("CD3D", "MS4A1", "NKG7", "LST1")
marker_status <- data.frame(
  gene = markers,
  selected_at_2000 = markers %in% hvg_2000,
  rank = match(markers, hvg_2000)
)
print(marker_status)

hvg_plot <- LabelPoints(
  plot = VariableFeaturePlot(obj_2000),
  points = markers,
  repel = TRUE
)
print(hvg_plot)

# Interpretation prompts:
# - Which familiar markers entered the set?
# - For any omitted marker, what can and cannot be concluded?
# - How does standardized variability differ from raw variance?

# Part 4: change only nfeatures --------------------------------------------

hvg_counts <- c(500, 1000, 2000, 3000, 5000)
hvg_objects <- lapply(hvg_counts, function(n) {
  FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = n,
    verbose = FALSE
  )
})
names(hvg_objects) <- paste0("hvg_", hvg_counts)
hvg_sets <- lapply(hvg_objects, VariableFeatures)

stopifnot(all(vapply(hvg_sets, length, integer(1)) == hvg_counts))

pair_grid <- expand.grid(
  set_1 = names(hvg_sets),
  set_2 = names(hvg_sets),
  stringsAsFactors = FALSE
)
pair_grid$jaccard <- mapply(
  function(a, b) jaccard_similarity(hvg_sets[[a]], hvg_sets[[b]]),
  pair_grid$set_1,
  pair_grid$set_2
)

jaccard_plot <- ggplot(pair_grid, aes(set_1, set_2, fill = jaccard)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", jaccard))) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(x = NULL, y = NULL, fill = "Jaccard") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

print(jaccard_plot)
ggsave(
  "outputs/module-04-hvg-jaccard.png",
  jaccard_plot,
  width = 7,
  height = 6,
  dpi = 150
)

# Part 5: carry 500 versus 5,000 features into PCA -------------------------

run_hvg_pca <- function(feature_object) {
  feature_object <- ScaleData(
    feature_object,
    features = VariableFeatures(feature_object),
    verbose = FALSE
  )
  RunPCA(
    feature_object,
    features = VariableFeatures(feature_object),
    npcs = 20,
    seed.use = 20260907,
    verbose = FALSE
  )
}

pca_500 <- run_hvg_pca(hvg_objects$hvg_500)
pca_5000 <- run_hvg_pca(hvg_objects$hvg_5000)

pca_plot <-
  DimPlot(pca_500, reduction = "pca", group.by = "truth_cell_type", raster = FALSE, combine = FALSE)[[1]] +
    ggtitle("500 HVGs") +
  DimPlot(pca_5000, reduction = "pca", group.by = "truth_cell_type", raster = FALSE, combine = FALSE)[[1]] +
    ggtitle("5,000 HVGs") +
  plot_layout(guides = "collect")

print(pca_plot)
ggsave(
  "outputs/module-04-hvg-pca.png",
  pca_plot,
  width = 11,
  height = 5,
  dpi = 150
)

score_correlations <- c(
  PC1 = cor(Embeddings(pca_500, "pca")[, 1], Embeddings(pca_5000, "pca")[, 1]),
  PC2 = cor(Embeddings(pca_500, "pca")[, 2], Embeddings(pca_5000, "pca")[, 2])
)
print(score_correlations)
print(rank_pca_loadings(pca_500, pcs = 1:3, n = 5))
print(rank_pca_loadings(pca_5000, pcs = 1:3, n = 5))

# Final interpretation prompts:
# - What changed when nfeatures changed, and what was held fixed?
# - Which broad relationships persist across the two PCA representations?
# - Which loading programs or lower-level relationships change?
# - Does an opposite PC sign indicate contradictory biology?
# - Recommend one feature count and state both its benefit and its risk.

message("Module 4 student lab executed successfully.")


# Explicit overlap and newly admitted genes, relative to the preceding cutoff.
overlap_counts <- data.frame(nfeatures = hvg_counts,
  overlap_with_500 = vapply(hvg_sets, function(x) length(intersect(x, hvg_sets[[1]])), integer(1)),
  added_since_previous = c(length(hvg_sets[[1]]), vapply(2:5, function(i)
    length(setdiff(hvg_sets[[i]], hvg_sets[[i - 1]])), integer(1))))
print(overlap_counts)
# Label top-ranked features as well as canonical markers.
top_hvg_plot <- LabelPoints(VariableFeaturePlot(obj_2000), points = head(hvg_2000, 5),
  repel = TRUE, xnudge = 0, ynudge = 0)
# Schematic, not Seurat's fitted curve: compare genes at the same expected mean.
schematic_data <- data.frame(mean = c(1, 1, 5, 5), variance = c(1, 4, 5, 8))
schematic_plot <- ggplot(schematic_data, aes(mean, variance)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + geom_point(size = 3) +
  labs(title = "Schematic: variation relative to an expected trend",
    subtitle = "Dashed line is illustrative, not a fitted Seurat model",
    x = "Mean", y = "Variance")
print(schematic_plot)
