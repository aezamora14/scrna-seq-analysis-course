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
