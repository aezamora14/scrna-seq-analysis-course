# Module 6: PCA geometry, loadings, scores, and population dependence.
# Run from the repository root; simulation labels are known only in teaching data.
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(patchwork)})
source("scripts/course_helpers.R")
set.seed(20260907)
dir.create("outputs", showWarnings = FALSE)

# Part 1: six genes measured in eight cells. Values are normalized expression.
toy_pca <- rbind(
  CCR7 = c(5, 4.6, 4.2, 3.8, 2, 1.5, 1, 0.5),
  IL7R = c(4.8, 4.5, 4, 3.9, 2.1, 1.4, 1.2, 0.4),
  TCF7 = c(4.4, 4.8, 3.8, 4, 1.8, 1.6, 0.8, 0.7),
  NKG7 = c(0.5, 0.9, 1.3, 1.8, 3.8, 4.2, 4.8, 5),
  GNLY = c(0.4, 1.1, 1.2, 1.5, 3.5, 4.5, 4.7, 5.2),
  PRF1 = c(0.8, 0.7, 1.6, 1.4, 4, 4.1, 5.1, 4.8)
)
colnames(toy_pca) <- paste0("Cell_", 1:8)
print(toy_pca)
Z <- t(scale(t(toy_pca)))
print(round(Z, 2))
# prcomp expects cells as rows; genes are columns.
toy_fit <- prcomp(t(Z), center = FALSE, scale. = FALSE)
L <- toy_fit$rotation
S <- toy_fit$x
print(round(L, 3))
print(round(S, 3))
# A score is a weighted sum over genes, not the loading of one gene.
contributions <- Z[, "Cell_1"] * L[, "PC1"]
print(contributions)
stopifnot(abs(sum(contributions) - S["Cell_1", "PC1"]) < 1e-10)
stopifnot(max(abs(t(Z) %*% L - S)) < 1e-10)
stopifnot(max(abs(crossprod(L) - diag(ncol(L)))) < 1e-10)
stopifnot(max(abs(cov(S) - diag(diag(cov(S))))) < 1e-10)

# Two-dimensional geometry first: arrows are perpendicular loading directions.
two <- t(Z[c("CCR7", "NKG7"), ])
two_fit <- prcomp(two, center = FALSE, scale. = FALSE)
arrows <- data.frame(x = 0, y = 0,
  xend = two_fit$rotation[1, ] * 2,
  yend = two_fit$rotation[2, ] * 2,
  PC = c("PC1", "PC2"))
geometry_plot <- ggplot(as.data.frame(two), aes(CCR7, NKG7)) +
  geom_point(size = 3) +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend, color = PC),
    inherit.aes = FALSE, arrow = grid::arrow(length = grid::unit(0.15, "inches"))) +
  coord_equal(xlim = c(-2, 2), ylim = c(-2, 2)) +
  labs(title = "Two genes, two orthogonal directions", x = "Scaled CCR7", y = "Scaled NKG7")
print(geometry_plot)
ggsave("outputs/module-06-geometry.png", geometry_plot, width = 7, height = 5)

# Perturbation: reverse one PC's sign. Scores and loadings must reverse together.
flipped_L <- L
flipped_S <- S
flipped_L[, 1] <- -L[, 1]
flipped_S[, 1] <- -S[, 1]
stopifnot(max(abs(S %*% t(L) - flipped_S %*% t(flipped_L))) < 1e-10)
# Compare reconstruction with one PC and all PCs.
reconstruction_1 <- S[, 1, drop = FALSE] %*% t(L[, 1, drop = FALSE])
reconstruction_all <- S %*% t(L)
print(c(one_PC_error = sum((t(Z) - reconstruction_1)^2),
        all_PC_error = sum((t(Z) - reconstruction_all)^2)))

# Part 2: inspect the actual Seurat matrices rather than assigning genes to PCs.
whole <- prepare_representation_object(nfeatures = 2000, npcs = 50)
pc_loadings <- Loadings(whole, reduction = "pca")
pc_scores <- Embeddings(whole, reduction = "pca")
print(dim(pc_loadings))
print(dim(pc_scores))
print(pc_loadings[1:6, 1:5])
print(pc_scores[1:6, 1:5])
input <- LayerData(whole, layer = "scale.data")[rownames(pc_loadings), , drop = FALSE]
projected <- t(input) %*% pc_loadings
stopifnot(max(abs(projected - pc_scores)) < 1e-4)

# Rank both ends independently for PCs 1--5.
loading_tables <- do.call(rbind, lapply(1:5, function(k) {
  weights <- pc_loadings[, k]
  positive <- head(sort(weights[weights > 0], decreasing = TRUE), 6)
  negative <- head(sort(weights[weights < 0]), 6)
  data.frame(PC = paste0("PC", k), gene = c(names(positive), names(negative)),
    loading = c(positive, negative), direction = rep(c("positive", "negative"),
      c(length(positive), length(negative))))
}))
print(loading_tables)
write.csv(loading_tables, "outputs/module-06-loadings.csv", row.names = FALSE)
loading_plot <- ggplot(loading_tables, aes(loading, reorder(gene, loading), fill = direction)) +
  geom_col() + facet_wrap(~PC, scales = "free_y", ncol = 2) +
  labs(x = "Gene loading", y = "Gene") + theme(legend.position = "bottom")
print(loading_plot)
ggsave("outputs/module-06-loadings.png", loading_plot, width = 10, height = 11)

# Total variance uses every PCA-input gene, not just the 50 retained PCs.
variance_whole <- pca_variance_table(whole)
print(head(variance_whole))
stopifnot(sum(variance_whole$percent_variance) <= 100 + 1e-6)
variance_plot <- ggplot(variance_whole, aes(PC, percent_variance)) +
  geom_line() + geom_point() + labs(y = "% of total PCA-input variance")
print(variance_plot)
ggsave("outputs/module-06-variance.png", variance_plot, width = 7, height = 4.5)

# Part 3: population perturbation. Relearn HVGs, scaling, and PCA within T cells.
t_cells <- Cells(whole)[whole$truth_cell_type %in% c("CD4 T", "CD8 T")]
t_only <- prepare_representation_object(subset(load_representation_object(), cells = t_cells),
  nfeatures = 2000, npcs = 50)
print(c(HVG_overlap = length(intersect(VariableFeatures(whole), VariableFeatures(t_only))),
        HVG_jaccard = jaccard_similarity(VariableFeatures(whole), VariableFeatures(t_only))))
print(rank_pca_loadings(t_only, pcs = 1:5, n = 6))
print(head(pca_variance_table(t_only)))
# A cross-PC correlation table allows axes to move or reverse sign.
score_comparison <- cor(Embeddings(whole, "pca")[t_cells, 1:5],
                        Embeddings(t_only, "pca")[t_cells, 1:5])
print(round(score_comparison, 2))
subset_plot <- DimPlot(subset(whole, cells = t_cells), reduction = "pca", group.by = "truth_t_state", combine = FALSE)[[1]] +
  ggtitle("T cells on whole-immune axes") +
  DimPlot(t_only, reduction = "pca", group.by = "truth_t_state", combine = FALSE)[[1]] +
  ggtitle("T cells on relearned T-cell axes") + plot_layout(guides = "collect")
print(subset_plot)
ggsave("outputs/module-06-subset.png", subset_plot, width = 11, height = 5)
stopifnot(identical(LayerData(whole, layer = "counts")[, t_cells],
                    LayerData(t_only, layer = "counts")))
# Interpret: which axes changed, which programs persisted, and why?
# Does orthogonality imply independent biological pathways? Are PC numbers portable?
message("Module 6 student lab executed successfully.")
