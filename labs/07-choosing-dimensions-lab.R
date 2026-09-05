# Module 7: hold the PCA object and resolution fixed; vary only retained dimensions.
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(patchwork)})
source("scripts/course_helpers.R")
set.seed(20260907)
dir.create("outputs", showWarnings = FALSE)
obj <- prepare_representation_object(nfeatures = 2000, npcs = 50)
original_scores <- Embeddings(obj, "pca")
original_counts <- LayerData(obj, layer = "counts")
print(ElbowPlot(obj, ndims = 50))
print(rank_pca_loadings(obj, pcs = c(10, 20, 30, 40), n = 6))

# A tiny distance example: a difference on PC3 is invisible if only PC1:2 are used.
coordinates <- rbind(A = c(0, 0, 0), B = c(0.1, 0, 4), C = c(1, 0, 0))
print(as.matrix(dist(coordinates[, 1:2])))
print(as.matrix(dist(coordinates[, 1:3])))

dims_values <- c(10, 20, 30, 40, 50)
fixed_resolution <- 0.5
# Keep this loop explicit: students can see which graph each partition uses.
for (n_dims in dims_values) {
  nn_name <- paste0("nn_dims_", n_dims)
  snn_name <- paste0("snn_dims_", n_dims)
  cluster_name <- paste0("clusters_dims_", n_dims)
  obj <- FindNeighbors(obj, reduction = "pca", dims = seq_len(n_dims),
    graph.name = c(nn_name, snn_name), nn.method = "rann", verbose = FALSE)
  obj <- FindClusters(obj, graph.name = snn_name, resolution = fixed_resolution,
    cluster.name = cluster_name, random.seed = 20260907, verbose = FALSE)
}
# Fix one display embedding so only cluster colors change across panels.
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, seed.use = 20260907,
  reduction.name = "umap_display", verbose = FALSE)
plots <- lapply(dims_values, function(n) {
  DimPlot(obj, reduction = "umap_display", group.by = paste0("clusters_dims_", n)) +
    ggtitle(paste("Graph PCs 1–", n, " | resolution 0.5", sep = "")) + NoLegend()
})
dims_plot <- wrap_plots(plots, ncol = 2) +
  plot_annotation(caption = "Same cells and display coordinates in every panel; cluster numbers are local labels.")
print(dims_plot)
ggsave("outputs/module-07-dims.png", dims_plot, width = 11, height = 12)
stopifnot(identical(original_scores, Embeddings(obj, "pca")),
          identical(original_counts, LayerData(obj, layer = "counts")))

cluster_summary <- data.frame(dims = dims_values, resolution = fixed_resolution,
  clusters = vapply(dims_values, function(n) length(unique(obj[[paste0("clusters_dims_", n)]][, 1])), integer(1)))
print(cluster_summary)
write.csv(cluster_summary, "outputs/module-07-cluster-counts.csv", row.names = FALSE)
contingency <- table(PC10 = obj$clusters_dims_10, PC40 = obj$clusters_dims_40)
print(contingency)
print(prop.table(contingency, 1))
for (n in dims_values) {
  print(n)
  print(table(cluster = obj[[paste0("clusters_dims_", n)]][, 1], lineage = obj$truth_cell_type))
}
# Optional: ARI compares partitions without assuming numeric labels correspond.
ari <- outer(dims_values, dims_values, Vectorize(function(a, b) {
  adjusted_rand_index(obj[[paste0("clusters_dims_", a)]][, 1],
                      obj[[paste0("clusters_dims_", b)]][, 1])
}))
dimnames(ari) <- list(dims_values, dims_values)
print(round(ari, 3))
stopifnot(max(abs(diag(ari) - 1)) < 1e-10)

# Connect graph changes to gene programs using normalized expression, not cluster names.
markers <- c("CD3D", "MS4A1", "NKG7", "LST1", "ISG15", "IFIT1", "MKI67")
Idents(obj) <- "clusters_dims_40"
marker_plot <- DotPlot(obj, features = markers) + RotatedAxis() +
  labs(title = "Marker evidence for the 40-PC partition")
print(marker_plot)
ggsave("outputs/module-07-markers.png", marker_plot, width = 9, height = 5)
# Identify any split by its contingency membership; then inspect marker coherence.
# Do not assume an interferon population must appear at PC40 in this simulation.
# Submit a justified cutoff, an alternative, and evidence that could change your choice.
# Complete the HVG -> scaling -> PCA -> dims -> graph input/output/decision table.
message("Module 7 student lab executed successfully.")
