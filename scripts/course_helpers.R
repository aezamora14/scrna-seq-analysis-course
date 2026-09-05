suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

load_teaching_object <- function() {
  counts <- readRDS("data/synthetic_qc_counts.rds")
  metadata <- read.csv(
    "data/synthetic_qc_metadata.csv",
    check.names = FALSE,
    na.strings = ""
  )
  rownames(metadata) <- metadata$cell_id

  stopifnot(identical(colnames(counts), rownames(metadata)))

  CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    min.cells = 0,
    min.features = 0,
    project = "controlled_qc_course"
  )
}

load_representation_object <- function() {
  counts <- readRDS("data/synthetic_representation_counts.rds")
  metadata <- read.csv(
    "data/synthetic_representation_metadata.csv",
    check.names = FALSE,
    na.strings = ""
  )
  rownames(metadata) <- metadata$cell_id

  stopifnot(identical(colnames(counts), rownames(metadata)))

  obj <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    min.cells = 0,
    min.features = 0,
    project = "controlled_representation_course"
  )
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
}

prepare_representation_object <- function(
  obj = load_representation_object(),
  nfeatures = 2000,
  npcs = 50,
  vars.to.regress = NULL,
  seed = 20260907
) {
  set.seed(seed)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", verbose = FALSE)
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = nfeatures,
    verbose = FALSE
  )
  obj <- ScaleData(obj, vars.to.regress = vars.to.regress, verbose = FALSE)
  obj <- RunPCA(
    obj,
    features = VariableFeatures(obj),
    npcs = npcs,
    seed.use = seed,
    verbose = FALSE
  )
  obj
}

jaccard_similarity <- function(x, y) {
  union_size <- length(union(x, y))
  if (union_size == 0L) {
    return(NA_real_)
  }
  length(intersect(x, y)) / union_size
}

pca_variance_table <- function(obj, reduction = "pca") {
  standard_deviation <- Stdev(obj, reduction = reduction)
  variance <- standard_deviation^2
  genes <- rownames(Loadings(obj, reduction = reduction))
  scaled <- LayerData(obj, assay = DefaultAssay(obj[[reduction]]), layer = "scale.data")
  total_variance <- sum(apply(scaled[genes, , drop = FALSE], 1, var))
  data.frame(
    PC = seq_along(standard_deviation),
    standard_deviation = standard_deviation,
    variance = variance,
    percent_variance = 100 * variance / total_variance,
    cumulative_percent = 100 * cumsum(variance) / total_variance,
    row.names = NULL
  )
}

rank_pca_loadings <- function(obj, pcs = 1:5, n = 8, reduction = "pca") {
  loadings <- Loadings(obj, reduction = reduction)
  selected <- intersect(paste0("PC_", pcs), colnames(loadings))
  do.call(
    rbind,
    lapply(selected, function(pc) {
      values <- loadings[, pc]
      positive <- head(order(values, decreasing = TRUE), n)
      negative <- head(order(values, decreasing = FALSE), n)
      data.frame(
        PC = pc,
        direction = rep(c("positive", "negative"), each = n),
        rank = rep(seq_len(n), 2),
        gene = c(names(values)[positive], names(values)[negative]),
        loading = c(values[positive], values[negative]),
        row.names = NULL
      )
    })
  )
}

adjusted_rand_index <- function(x, y) {
  contingency <- table(x, y)
  choose_two <- function(z) z * (z - 1) / 2
  cell_pairs <- sum(choose_two(contingency))
  row_pairs <- sum(choose_two(rowSums(contingency)))
  column_pairs <- sum(choose_two(colSums(contingency)))
  total_pairs <- choose_two(sum(contingency))
  if (total_pairs == 0) return(NA_real_)
  expected <- row_pairs * column_pairs / total_pairs
  denominator <- 0.5 * (row_pairs + column_pairs) - expected
  if (denominator == 0) {
    # Identical all-in-one or all-singleton partitions have perfect agreement.
    return(1)
  }
  (cell_pairs - expected) / denominator
}

run_dims_sensitivity <- function(
  obj,
  dims_values = c(10, 20, 30, 40, 50),
  resolution = 0.5,
  seed = 20260907
) {
  for (n_dims in dims_values) {
    nn_name <- paste0("nn_dims_", n_dims)
    snn_name <- paste0("snn_dims_", n_dims)
    cluster_name <- paste0("clusters_dims_", n_dims)
    reduction_name <- paste0("umap_dims_", n_dims)

    obj <- FindNeighbors(
      obj,
      reduction = "pca",
      dims = seq_len(n_dims),
      graph.name = c(nn_name, snn_name),
      verbose = FALSE
    )
    obj <- FindClusters(
      obj,
      graph.name = snn_name,
      resolution = resolution,
      cluster.name = cluster_name,
      random.seed = seed,
      verbose = FALSE
    )
    obj <- RunUMAP(
      obj,
      reduction = "pca",
      dims = seq_len(n_dims),
      reduction.name = reduction_name,
      reduction.key = paste0("UMAP", n_dims, "_"),
      seed.use = seed,
      verbose = FALSE
    )
  }
  obj
}

add_qc_metrics <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

  metadata <- obj[[]]
  obj$log10_genes_per_umi <-
    log10(pmax(metadata$nFeature_RNA, 1)) /
    log10(pmax(metadata$nCount_RNA, 2))

  obj
}

qc_strategy_flags <- function(obj) {
  metadata <- obj[[]]

  list(
    permissive = with(
      metadata,
      nCount_RNA >= 50 &
        nFeature_RNA >= 20 &
        percent.mt < 35
    ),
    moderate = with(
      metadata,
      nCount_RNA >= 300 &
        nFeature_RNA >= 80 &
        percent.mt < 20
    ),
    aggressive = with(
      metadata,
      nCount_RNA >= 900 &
        nFeature_RNA >= 130 &
        percent.mt < 10
    )
  )
}

run_qc_pipeline <- function(obj, keep, seed = 20260904) {
  set.seed(seed)
  filtered <- subset(obj, cells = Cells(obj)[keep])

  filtered <- NormalizeData(filtered, verbose = FALSE)
  filtered <- FindVariableFeatures(
    filtered,
    selection.method = "vst",
    nfeatures = 180,
    verbose = FALSE
  )
  filtered <- ScaleData(filtered, verbose = FALSE)
  filtered <- RunPCA(filtered, npcs = 20, verbose = FALSE)
  filtered <- FindNeighbors(filtered, dims = 1:15, verbose = FALSE)
  filtered <- FindClusters(
    filtered,
    resolution = 0.5,
    random.seed = seed,
    verbose = FALSE
  )
  filtered <- RunUMAP(
    filtered,
    dims = 1:15,
    seed.use = seed,
    verbose = FALSE
  )

  filtered
}
