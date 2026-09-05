# Module 5 lab: scaling and regression
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

# Part 1: gene-wise z scores -----------------------------------------------

toy_scaled <- rbind(
  NKG7 = c(0.5, 0.8, 2.1, 4.6),
  CD3D = c(2.7, 2.9, 3.0, 3.2),
  MALAT1 = c(6.8, 6.9, 6.8, 6.9)
)
colnames(toy_scaled) <- paste0("Cell_", LETTERS[1:4])

x <- toy_scaled["NKG7", ]
manual_mean <- sum(x) / length(x)
manual_sd <- sqrt(sum((x - manual_mean)^2) / (length(x) - 1))
manual_z <- (x["Cell_D"] - manual_mean) / manual_sd

manual_scaled <- t(apply(toy_scaled, 1, function(values) {
  (values - mean(values)) / sd(values)
}))
base_scaled <- t(scale(t(toy_scaled)))

stopifnot(
  isTRUE(all.equal(unname(manual_z), unname(base_scaled["NKG7", "Cell_D"]))),
  identical(dimnames(manual_scaled), dimnames(base_scaled)),
  isTRUE(all.equal(as.vector(manual_scaled), as.vector(base_scaled), tolerance = 1e-12))
)
print(round(base_scaled, 3))

matrix_long <- function(matrix, stage) {
  as.data.frame(matrix) |>
    tibble::rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "cell", values_to = "value") |>
    mutate(stage = stage)
}

scale_long <- bind_rows(
  matrix_long(toy_scaled, "Normalized expression"),
  matrix_long(base_scaled, "Gene-wise scaled expression")
)

normalized_heatmap <- ggplot(subset(scale_long, stage == "Normalized expression"), aes(cell, gene, fill = value)) +
  geom_tile(color = "white") + geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient(low = "white", high = "skyblue3", name = "Normalized") +
  labs(title = "Normalized expression", x = "Cell", y = "Gene")
scaled_heatmap <- ggplot(subset(scale_long, stage == "Gene-wise scaled expression"), aes(cell, gene, fill = value)) +
  geom_tile(color = "white") + geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient2(low = "skyblue3", mid = "white", high = "salmon", midpoint = 0, name = "Scaled") +
  labs(title = "Gene-wise scaled expression", x = "Cell", y = "Gene")
scale_plot <- normalized_heatmap + scaled_heatmap

print(scale_plot)
ggsave(
  "outputs/module-05-normalized-scaled.png",
  scale_plot,
  width = 10,
  height = 4.5,
  dpi = 150
)

# Interpretation prompts:
# - What do positive and negative scaled values mean?
# - Why is a scaled NKG7 value of +2 not four times a CD3D value of +0.5?
# - Which original data layers would you preserve?

# Part 2: PCA with and without variance scaling ----------------------------

unscaled_example <- data.frame(
  cell = paste0("Cell_", 1:10),
  group = rep(c("naive-like", "cytotoxic"), each = 5),
  large_range_gene = c(20, 80, 35, 95, 45, 25, 85, 40, 100, 50),
  coherent_state_gene = c(1.0, 1.2, 0.8, 1.1, 0.9, 4.8, 5.2, 5.0, 4.9, 5.1)
)

unscaled_pca <- prcomp(
  unscaled_example[, c("large_range_gene", "coherent_state_gene")],
  center = TRUE,
  scale. = FALSE
)
scaled_pca <- prcomp(
  unscaled_example[, c("large_range_gene", "coherent_state_gene")],
  center = TRUE,
  scale. = TRUE
)

unscaled_example$PC1_unscaled <- unscaled_pca$x[, 1]
unscaled_example$PC1_scaled <- scaled_pca$x[, 1]

print(unscaled_pca$rotation)
print(scaled_pca$rotation)

p_unscaled <- ggplot(unscaled_example, aes(group, PC1_unscaled, color = group)) +
  geom_point(size = 3) +
  ggtitle("Without variance scaling")
p_scaled <- ggplot(unscaled_example, aes(group, PC1_scaled, color = group)) +
  geom_point(size = 3) +
  ggtitle("With variance scaling")
pca_scale_plot <- p_unscaled + p_scaled + plot_layout(guides = "collect")

print(pca_scale_plot)
ggsave(
  "outputs/module-05-pca-scaling.png",
  pca_scale_plot,
  width = 10,
  height = 4.5,
  dpi = 150
)

# Interpretation prompts:
# - Which gene dominates PC1 without variance scaling, and why?
# - Which coherent group signal becomes more visible after scaling?
# - Does scaling prove that either variable is biologically desirable?

# Part 3: reference-population dependence ---------------------------------

reference_obj <- NormalizeData(load_representation_object(), verbose = FALSE)
reference_features <- c("NKG7", "CD3D", "IL7R")
whole_scaled <- ScaleData(
  reference_obj,
  features = reference_features,
  verbose = FALSE
)

t_cells <- Cells(reference_obj)[reference_obj$truth_cell_type %in% c("CD4 T", "CD8 T")]
t_reference <- subset(reference_obj, cells = t_cells)
t_scaled <- ScaleData(
  t_reference,
  features = reference_features,
  verbose = FALSE
)
same_cell <- t_cells[1]

reference_comparison <- data.frame(
  gene = reference_features,
  normalized_value = as.numeric(
    LayerData(reference_obj, layer = "data")[reference_features, same_cell]
  ),
  scaled_whole_immune = as.numeric(
    LayerData(whole_scaled, layer = "scale.data")[reference_features, same_cell]
  ),
  scaled_T_cells_only = as.numeric(
    LayerData(t_scaled, layer = "scale.data")[reference_features, same_cell]
  )
)
print(reference_comparison)

stopifnot(any(abs(
  reference_comparison$scaled_whole_immune -
    reference_comparison$scaled_T_cells_only
) > 1e-6))

reference_long <- reference_comparison |>
  select(-normalized_value) |>
  pivot_longer(-gene, names_to = "reference", values_to = "scaled_value")

reference_plot <- ggplot(
  reference_long,
  aes(gene, scaled_value, fill = reference)
) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(
    title = paste("Same cell:", same_cell),
    x = "Gene",
    y = "Scaled value"
  )
print(reference_plot)

# Interpretation prompts:
# - What changed when the reference population changed?
# - What remained identical about the cell and its normalized expression?
# - Why is a z score contextual rather than intrinsic?

# Part 4: no regression versus two regression choices ---------------------

base_obj <- load_representation_object()
base_obj <- NormalizeData(base_obj, verbose = FALSE)
base_obj <- FindVariableFeatures(
  base_obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)
features <- VariableFeatures(base_obj)

run_scaled_pca <- function(obj, regress = NULL) {
  obj <- ScaleData(
    obj,
    features = features,
    vars.to.regress = regress,
    verbose = FALSE
  )
  RunPCA(
    obj,
    features = features,
    npcs = 30,
    seed.use = 20260907,
    verbose = FALSE
  )
}

no_regression <- run_scaled_pca(base_obj)
regress_mito <- run_scaled_pca(base_obj, "percent.mt")
regress_cycle <- run_scaled_pca(base_obj, "cell_cycle_score")

regression_plot <- wrap_plots(
  list(
    DimPlot(no_regression, reduction = "pca", group.by = "cell_cycle_phase") +
      ggtitle("No regression"),
    DimPlot(regress_mito, reduction = "pca", group.by = "cell_cycle_phase") +
      ggtitle("Regress percent.mt"),
    DimPlot(regress_cycle, reduction = "pca", group.by = "cell_cycle_phase") +
      ggtitle("Regress cell-cycle score")
  ),
  ncol = 3,
  guides = "collect"
)
print(regression_plot)
ggsave(
  "outputs/module-05-regression-pca.png",
  regression_plot,
  width = 13,
  height = 4.6,
  dpi = 150
)

pc_associations <- data.frame(
  analysis = c("none", "percent.mt", "cell_cycle_score"),
  PC1_cor_percent_mt = c(
    cor(Embeddings(no_regression, "pca")[, 1], no_regression$percent.mt),
    cor(Embeddings(regress_mito, "pca")[, 1], regress_mito$percent.mt),
    cor(Embeddings(regress_cycle, "pca")[, 1], regress_cycle$percent.mt)
  ),
  PC1_cor_cell_cycle = c(
    cor(Embeddings(no_regression, "pca")[, 1], no_regression$cell_cycle_score),
    cor(Embeddings(regress_mito, "pca")[, 1], regress_mito$cell_cycle_score),
    cor(Embeddings(regress_cycle, "pca")[, 1], regress_cycle$cell_cycle_score)
  )
)
print(pc_associations)

loading_comparison <- bind_rows(
  rank_pca_loadings(no_regression, pcs = 1:5, n = 5) |>
    mutate(analysis = "none"),
  rank_pca_loadings(regress_mito, pcs = 1:5, n = 5) |>
    mutate(analysis = "percent.mt"),
  rank_pca_loadings(regress_cycle, pcs = 1:5, n = 5) |>
    mutate(analysis = "cell_cycle_score")
)
print(loading_comparison)

stopifnot(
  all(c("counts", "data", "scale.data") %in% Layers(no_regression[["RNA"]])),
  all(LayerData(no_regression, layer = "counts") == LayerData(base_obj, layer = "counts")),
  all(LayerData(no_regression, layer = "data") == LayerData(base_obj, layer = "data"))
)

# Final interpretation prompts:
# - Which associations and loading programs changed under each regression?
# - Which broad lineage relationships remained?
# - Under what analysis goal would cell-cycle regression be defensible?
# - Under what analysis goal would it remove the response of interest?
# - Is percent.mt necessarily technical in real tissue data?
# - Recommend no regression or one regression and defend the remaining estimand.

message("Module 5 student lab executed successfully.")


# Calculate three values manually, preserving sample-SD arithmetic.
manual_three <- (toy_scaled["NKG7", 1:3] - manual_mean) / manual_sd
print(manual_three)
stopifnot(max(abs(manual_three - base_scaled["NKG7", 1:3])) < 1e-12)
# A transparent regression example: fit, predict, subtract, and scale residuals.
regression_toy <- data.frame(q = 0:5, expression = c(1, 3, 4, 7, 8, 10))
fit <- lm(expression ~ q, data = regression_toy)
regression_toy$fitted <- predict(fit)
regression_toy$residual <- regression_toy$expression - regression_toy$fitted
print(regression_toy)
print(scale(regression_toy$residual))
# Compare distributions and value mappings using identical cells and genes.
gene_values <- data.frame(
  gene = rep(reference_features, each = ncol(whole_scaled)),
  normalized = as.vector(t(as.matrix(LayerData(whole_scaled, layer = "data")[reference_features, ]))),
  scaled = as.vector(t(LayerData(whole_scaled, layer = "scale.data")[reference_features, ])))
mapping_plot <- ggplot(gene_values, aes(normalized, scaled)) + geom_point(alpha = 0.25) +
  facet_wrap(~gene) + labs(x = "Normalized expression", y = "Gene-wise scaled value")
distribution_data <- pivot_longer(gene_values, c(normalized, scaled), names_to = "stage", values_to = "value")
distribution_plot <- ggplot(distribution_data, aes(value, color = gene)) +
  geom_density() + facet_wrap(~stage, scales = "free") + labs(y = "Density")
print(mapping_plot)
print(distribution_plot)
# Residual gene distributions across all three preprocessing choices.
regression_distributions <- bind_rows(lapply(
  list(none = no_regression, percent_mt = regress_mito, cycle = regress_cycle),
  function(x) matrix_long(LayerData(x, layer = "scale.data")[c("MKI67", "NKG7"), ], "scaled")),
  .id = "analysis")
regression_distribution_plot <- ggplot(regression_distributions, aes(value, color = analysis)) +
  geom_density() + facet_wrap(~gene) + labs(x = "Scaled/residual value", y = "Density")
print(regression_distribution_plot)
