suppressPackageStartupMessages({
  library(Matrix)
})

set.seed(20260904)

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required to generate and verify the data manifest.")
}

data_dir <- "data"
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

toy_counts <- matrix(
  c(
    20, 40,  0,  2,
    10, 20, 50,  0,
     0,  0,  0, 30,
     0,  0, 10,  0,
    60,120, 30, 60,
    10, 20, 10,  8
  ),
  nrow = 6,
  byrow = TRUE,
  dimnames = list(
    c("CD3D", "NKG7", "MS4A1", "LST1", "MALAT1", "MT-CO1"),
    c("Cell_A", "Cell_B", "Cell_C", "Cell_D")
  )
)

mitochondrial_genes <- c(
  "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2",
  "MT-ATP6", "MT-CYB", "MT-ND4", "MT-ND5"
)

ribosomal_genes <- c(
  "RPLP0", "RPL3", "RPL7", "RPL13A", "RPL32", "RPL37A",
  "RPS3", "RPS6", "RPS12", "RPS18", "RPS27A", "RPS29"
)

housekeeping_genes <- c(
  "MALAT1", "ACTB", "GAPDH", "B2M", "EEF1A1",
  "TUBA1B", "FTL", "FTH1", "HSP90AB1", "UBC"
)

immune_genes <- c(
  "CD3D", "CD3E", "TRAC", "CD4", "IL7R", "CCR7", "LTB", "TCF7",
  "CD8A", "CD8B", "NKG7", "GZMK", "CCL5", "CTSW", "PRF1", "GZMB",
  "GNLY", "KLRD1", "FCGR3A", "TRDC",
  "MS4A1", "CD79A", "CD79B", "CD37", "CD74", "HLA-DRA", "CD19", "CD22",
  "LST1", "TYROBP", "FCER1G", "CTSS", "LYZ", "S100A8", "S100A9",
  "LILRB1", "CTSD", "LGALS3", "SAT1", "CTSZ"
)

state_genes <- c(
  "FOS", "JUN", "JUNB", "HSPA1A", "HSPA1B", "DNAJB1",
  "IFIT1", "IFIT3", "ISG15", "MX1", "IFI6", "STAT1"
)

ambient_genes <- c("HBB", "HBA1", "HBA2", "ALB", "APOA1")

named_genes <- unique(c(
  mitochondrial_genes,
  ribosomal_genes,
  housekeeping_genes,
  immune_genes,
  state_genes,
  ambient_genes
))

n_genes <- 260L
generic_genes <- sprintf("GENE%03d", seq_len(n_genes - length(named_genes)))
genes <- c(named_genes, generic_genes)

# A highly skewed background profile keeps the small teaching matrix sparse
# without adding thousands of uninformative rows.
base_weight <- rgamma(length(genes), shape = 0.35, rate = 2) + 0.001
names(base_weight) <- genes
base_weight[housekeeping_genes] <- c(8, 4, 4, 2, 3, 2, 2, 2, 2, 2)
base_weight[ribosomal_genes] <- 1.5
base_weight[mitochondrial_genes] <- 0.25
base_weight[immune_genes] <- 0.01
base_weight[state_genes] <- 0.02
base_weight[ambient_genes] <- 0.003

marker_sets <- list(
  `CD4 T` = c("CD3D", "CD3E", "TRAC", "CD4", "IL7R", "CCR7", "LTB", "TCF7"),
  `CD8 T` = c("CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "NKG7", "GZMK", "CCL5", "CTSW"),
  NK = c("NKG7", "GNLY", "KLRD1", "FCGR3A", "PRF1", "GZMB", "CTSW"),
  B = c("MS4A1", "CD79A", "CD79B", "CD37", "CD74", "HLA-DRA", "CD19", "CD22"),
  Monocyte = c("LST1", "TYROBP", "FCER1G", "CTSS", "LYZ", "S100A8", "S100A9", "LILRB1", "CTSD", "LGALS3")
)

cell_types <- names(marker_sets)
cell_type_prob <- c(0.27, 0.24, 0.16, 0.16, 0.17)
sample_ids <- c("S1", "S2", "S3", "S4")
sample_depth <- c(S1 = 0.75, S2 = 1.00, S3 = 1.45, S4 = 0.65)
cell_type_depth <- c(`CD4 T` = 0.95, `CD8 T` = 0.85, NK = 0.62, B = 0.75, Monocyte = 1.30)

n_singlet <- 640L
singlet_id <- sprintf("BC%04d-1", seq_len(n_singlet))
singlet_sample <- sample(sample_ids, n_singlet, replace = TRUE)
singlet_type <- sample(cell_types, n_singlet, replace = TRUE, prob = cell_type_prob)

qc_class <- rep("high_quality", n_singlet)
qc_class[sample(seq_len(n_singlet), 50)] <- "low_quality"
available <- which(qc_class == "high_quality")
qc_class[sample(available, 45)] <- "stressed"
available <- which(qc_class == "high_quality")
qc_class[sample(available, 45)] <- "ambient_high"

ambient_profile <- base_weight
ambient_profile[] <- 0.01
ambient_profile[c("MALAT1", "LYZ", "HBB", "HBA1", "HBA2", "ALB", "APOA1")] <-
  c(8, 7, 10, 8, 6, 4, 3)
ambient_profile[ribosomal_genes] <- 0.5
ambient_profile <- ambient_profile / sum(ambient_profile)

make_profile <- function(cell_type, sample_id, qc, condition) {
  weights <- base_weight
  weights[marker_sets[[cell_type]]] <- weights[marker_sets[[cell_type]]] + 8

  if (cell_type %in% c("CD4 T", "CD8 T")) {
    weights[c("CD3D", "CD3E", "TRAC")] <- weights[c("CD3D", "CD3E", "TRAC")] + 5
  }
  if (cell_type == "CD8 T") {
    weights[c("NKG7", "CCL5", "CTSW")] <- weights[c("NKG7", "CCL5", "CTSW")] + 3
  }
  if (condition == "stimulated") {
    weights[c("IFIT1", "IFIT3", "ISG15", "MX1", "IFI6", "STAT1")] <-
      weights[c("IFIT1", "IFIT3", "ISG15", "MX1", "IFI6", "STAT1")] + 2
  }
  if (qc == "stressed") {
    weights[state_genes[seq_len(6)]] <- weights[state_genes[seq_len(6)]] + 9
    weights[mitochondrial_genes] <- weights[mitochondrial_genes] * 2.5
  }
  if (qc == "low_quality") {
    weights[mitochondrial_genes] <- weights[mitochondrial_genes] * 16
    weights[ribosomal_genes] <- weights[ribosomal_genes] * 0.4
  }
  weights / sum(weights)
}

singlet_counts <- matrix(
  0L,
  nrow = length(genes),
  ncol = n_singlet,
  dimnames = list(genes, singlet_id)
)

for (i in seq_len(n_singlet)) {
  sample_id <- singlet_sample[i]
  condition <- if (sample_id %in% c("S1", "S2")) "control" else "stimulated"
  qc <- qc_class[i]
  median_library <- switch(
    qc,
    low_quality = 260,
    stressed = 1500,
    ambient_high = 900,
    high_quality = 2100
  )
  library_size <- max(
    40L,
    as.integer(round(rlnorm(
      1,
      log(median_library * sample_depth[sample_id] * cell_type_depth[singlet_type[i]]),
      0.35
    )))
  )
  profile <- make_profile(singlet_type[i], sample_id, qc, condition)
  singlet_counts[, i] <- as.integer(rmultinom(1, library_size, profile))

  if (qc == "ambient_high") {
    ambient_total <- max(40L, as.integer(round(library_size * 0.35)))
    singlet_counts[, i] <- singlet_counts[, i] +
      as.integer(rmultinom(1, ambient_total, ambient_profile))
  }
}

n_doublet <- 30L
doublet_id <- sprintf("BC%04d-1", n_singlet + seq_len(n_doublet))
doublet_pairs <- replicate(
  n_doublet,
  sample(which(qc_class == "high_quality"), 2, replace = FALSE)
)
doublet_counts <- vapply(
  seq_len(n_doublet),
  function(i) singlet_counts[, doublet_pairs[1, i]] + singlet_counts[, doublet_pairs[2, i]],
  integer(length(genes))
)
rownames(doublet_counts) <- genes
colnames(doublet_counts) <- doublet_id

n_empty <- 30L
empty_id <- sprintf("BC%04d-1", n_singlet + n_doublet + seq_len(n_empty))
empty_counts <- vapply(
  seq_len(n_empty),
  function(i) {
    total <- max(15L, as.integer(rpois(1, lambda = 70)))
    as.integer(rmultinom(1, total, ambient_profile))
  },
  integer(length(genes))
)
rownames(empty_counts) <- genes
colnames(empty_counts) <- empty_id

all_counts <- cbind(singlet_counts, doublet_counts, empty_counts)
storage.mode(all_counts) <- "integer"
sparse_counts <- as(all_counts, "dgCMatrix")

singlet_condition <- ifelse(singlet_sample %in% c("S1", "S2"), "control", "stimulated")

metadata <- data.frame(
  cell_id = colnames(sparse_counts),
  sample_id = c(
    singlet_sample,
    singlet_sample[doublet_pairs[1, ]],
    sample(sample_ids, n_empty, replace = TRUE)
  ),
  condition = c(
    singlet_condition,
    singlet_condition[doublet_pairs[1, ]],
    rep("empty", n_empty)
  ),
  truth_cell_type = c(
    singlet_type,
    rep("Doublet", n_doublet),
    rep("Empty droplet", n_empty)
  ),
  truth_qc_class = c(
    qc_class,
    rep("doublet", n_doublet),
    rep("empty_droplet", n_empty)
  ),
  doublet_source_1 = c(
    rep(NA_character_, n_singlet),
    singlet_type[doublet_pairs[1, ]],
    rep(NA_character_, n_empty)
  ),
  doublet_source_2 = c(
    rep(NA_character_, n_singlet),
    singlet_type[doublet_pairs[2, ]],
    rep(NA_character_, n_empty)
  ),
  stringsAsFactors = FALSE
)

expected_files <- c(
  toy = file.path(data_dir, "toy_counts.csv"),
  counts = file.path(data_dir, "synthetic_qc_counts.rds"),
  metadata = file.path(data_dir, "synthetic_qc_metadata.csv")
)

write_outputs <- function() {
  toy_table <- data.frame(gene = rownames(toy_counts), toy_counts, check.names = FALSE)
  write.csv(toy_table, expected_files[["toy"]], row.names = FALSE, na = "")
  saveRDS(sparse_counts, expected_files[["counts"]], version = 3)
  write.csv(metadata, expected_files[["metadata"]], row.names = FALSE, na = "")

  manifest <- data.frame(
    file = unname(expected_files),
    bytes = as.numeric(file.info(unname(expected_files))$size),
    sha256 = vapply(
      unname(expected_files),
      digest::digest,
      character(1),
      algo = "sha256",
      file = TRUE
    ),
    rows = c(nrow(toy_counts), nrow(sparse_counts), nrow(metadata)),
    columns = c(ncol(toy_counts), ncol(sparse_counts), ncol(metadata)),
    stringsAsFactors = FALSE
  )
  write.csv(
    manifest,
    file.path(data_dir, "synthetic_data_manifest.csv"),
    row.names = FALSE,
    na = ""
  )
}

check_outputs <- function() {
  missing <- unname(expected_files)[!file.exists(unname(expected_files))]
  if (length(missing)) {
    stop("Missing generated data files: ", paste(missing, collapse = ", "))
  }

  saved_toy <- read.csv(expected_files[["toy"]], row.names = 1, check.names = FALSE)
  saved_counts <- readRDS(expected_files[["counts"]])
  saved_metadata <- read.csv(
    expected_files[["metadata"]],
    check.names = FALSE,
    na.strings = ""
  )
  required_metadata <- c(
    "cell_id",
    "sample_id",
    "condition",
    "truth_cell_type",
    "truth_qc_class",
    "doublet_source_1",
    "doublet_source_2"
  )

  stopifnot(
    isTRUE(all.equal(as.matrix(saved_toy), toy_counts, check.attributes = FALSE)),
    inherits(saved_counts, "dgCMatrix"),
    identical(dim(saved_counts), c(260L, 700L)),
    identical(rownames(saved_counts), genes),
    identical(colnames(saved_counts), saved_metadata$cell_id),
    identical(colnames(saved_metadata), required_metadata),
    nrow(saved_metadata) == ncol(saved_counts),
    all(is.finite(saved_counts@x)),
    all(saved_counts@x >= 0),
    all(saved_counts@x == floor(saved_counts@x)),
    setequal(
      unique(saved_metadata$truth_qc_class),
      c("high_quality", "low_quality", "stressed", "ambient_high", "doublet", "empty_droplet")
    )
  )

  manifest_path <- file.path(data_dir, "synthetic_data_manifest.csv")
  if (!file.exists(manifest_path)) {
    stop("Missing data manifest: ", manifest_path)
  }
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  current_hashes <- vapply(
    manifest$file,
    digest::digest,
    character(1),
    algo = "sha256",
    file = TRUE
  )
  stopifnot(identical(unname(current_hashes), unname(manifest$sha256)))

  message(
    "Synthetic data check passed: ",
    nrow(saved_counts), " genes × ", ncol(saved_counts), " barcodes."
  )
}

if (check_only) {
  check_outputs()
} else {
  write_outputs()
  check_outputs()
  message("Synthetic teaching data generated deterministically.")
}
