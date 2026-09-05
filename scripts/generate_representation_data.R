suppressPackageStartupMessages({
  library(Matrix)
})

set.seed(20260907)

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required to generate and verify the data manifest.")
}

data_dir <- "data"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# This second controlled dataset begins after QC and normalization. Its larger
# gene dimension makes 500--5,000 feature-selection comparisons possible while
# remaining small enough for an ordinary laptop.
mitochondrial_genes <- c(
  "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP6", "MT-CYB",
  "MT-ND4", "MT-ND5", "MT-ND6", "MT-CO3"
)

housekeeping_genes <- c(
  "MALAT1", "ACTB", "GAPDH", "B2M", "EEF1A1", "RPLP0", "RPL13A",
  "RPS18", "TUBA1B", "UBC"
)

lineage_sets <- list(
  `CD4 T` = c("CD3D", "CD3E", "TRAC", "CD4", "LTB", "MAL"),
  `CD8 T` = c("CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "CTSW"),
  NK = c("NKG7", "GNLY", "KLRD1", "PRF1", "GZMB", "FCGR3A"),
  B = c("MS4A1", "CD79A", "CD79B", "CD37", "CD74", "HLA-DRA"),
  Monocyte = c("LST1", "TYROBP", "FCER1G", "CTSS", "LYZ", "S100A8", "S100A9")
)

naive_genes <- c("CCR7", "IL7R", "TCF7", "LEF1", "LTB", "MAL")
cytotoxic_genes <- c("NKG7", "GNLY", "PRF1", "GZMB", "CCL5", "CTSW")
ifn_genes <- c("IFIT1", "IFIT3", "ISG15", "MX1", "IFI6", "STAT1")
s_genes <- c("MCM5", "PCNA", "TYMS", "MCM2", "MCM4", "RRM1")
g2m_genes <- c("MKI67", "TOP2A", "NUSAP1", "UBE2C", "BIRC5", "CDK1")
stress_genes <- c("FOS", "JUN", "JUNB", "HSPA1A", "HSPA1B", "DNAJB1")

named_genes <- unique(c(
  mitochondrial_genes,
  housekeeping_genes,
  unlist(lineage_sets, use.names = FALSE),
  naive_genes,
  cytotoxic_genes,
  ifn_genes,
  s_genes,
  g2m_genes,
  stress_genes
))

n_genes <- 6000L
generic_genes <- sprintf("GENE%04d", seq_len(n_genes - length(named_genes)))
genes <- c(named_genes, generic_genes)

cell_type_counts <- c(`CD4 T` = 140L, `CD8 T` = 140L, NK = 70L, B = 65L, Monocyte = 65L)
truth_cell_type <- rep(names(cell_type_counts), cell_type_counts)
n_cells <- length(truth_cell_type)
cell_id <- sprintf("RC%04d-1", seq_len(n_cells))

# Randomize row order without changing the controlled group sizes.
cell_order <- sample(seq_len(n_cells))
truth_cell_type <- truth_cell_type[cell_order]

sample_levels <- c(
  "D1_control", "D1_stimulated", "D2_control",
  "D2_stimulated", "D3_control", "D3_stimulated"
)
sample_id <- rep(sample_levels, length.out = n_cells)[sample(seq_len(n_cells))]
donor_id <- sub("_(control|stimulated)$", "", sample_id)
condition <- sub("^(D[0-9]+)_", "", sample_id)

truth_t_state <- rep("not_T", n_cells)
t_cells <- which(truth_cell_type %in% c("CD4 T", "CD8 T"))
truth_t_state[t_cells] <- ifelse(
  truth_cell_type[t_cells] == "CD4 T",
  sample(c("naive-like", "memory-like"), length(t_cells), replace = TRUE, prob = c(0.62, 0.38)),
  sample(c("naive-like", "cytotoxic"), length(t_cells), replace = TRUE, prob = c(0.43, 0.57))
)

ifn_score <- rnorm(n_cells, mean = ifelse(condition == "stimulated", 0.85, -0.25), sd = 0.55)
cycling <- rbinom(n_cells, 1, prob = 0.24)
cell_cycle_score <- pmax(0, rnorm(n_cells, mean = 1.5 * cycling, sd = 0.35))
cell_cycle_phase <- ifelse(
  cycling == 0,
  "G1",
  sample(c("S", "G2M"), n_cells, replace = TRUE)
)
stress_score <- pmax(0, rnorm(n_cells, mean = 0.25, sd = 0.45))
mitochondrial_technical_score <- pmax(0, rnorm(n_cells, mean = 0.35, sd = 0.5))

# Weak independent programs create realistic lower-ranked axes. Some later PCs
# have coherent loadings, while sufficiently late PCs are mostly sampling noise.
n_weak_programs <- 12L
weak_program_size <- 18L
weak_program_genes <- split(
  generic_genes[seq_len(n_weak_programs * weak_program_size)],
  rep(seq_len(n_weak_programs), each = weak_program_size)
)
weak_scores <- matrix(rnorm(n_cells * n_weak_programs), nrow = n_cells)

base_weight <- rgamma(n_genes, shape = 0.28, rate = 1.1) + 0.0005
names(base_weight) <- genes
base_weight[housekeeping_genes] <- c(15, 8, 8, 5, 8, 5, 5, 5, 5, 5)
base_weight[mitochondrial_genes] <- 1.1
base_weight[unique(unlist(lineage_sets, use.names = FALSE))] <- 0.03
base_weight[unique(c(naive_genes, cytotoxic_genes, ifn_genes, s_genes, g2m_genes, stress_genes))] <- 0.02

add_program <- function(weights, genes, strength) {
  weights[genes] <- weights[genes] + pmax(strength, 0)
  weights
}

counts <- matrix(
  0L,
  nrow = n_genes,
  ncol = n_cells,
  dimnames = list(genes, cell_id)
)

for (i in seq_len(n_cells)) {
  weights <- base_weight
  cell_type <- truth_cell_type[i]
  weights <- add_program(weights, lineage_sets[[cell_type]], 14)

  if (cell_type %in% c("CD4 T", "CD8 T")) {
    weights <- add_program(weights, c("CD3D", "CD3E", "TRAC"), 10)
  }
  if (truth_t_state[i] == "naive-like") {
    weights <- add_program(weights, naive_genes, 12)
  }
  if (truth_t_state[i] == "memory-like") {
    weights <- add_program(weights, c("IL7R", "LTB", "MAL", "CTSW"), 6)
  }
  if (truth_t_state[i] == "cytotoxic") {
    weights <- add_program(weights, cytotoxic_genes, 13)
  }

  weights <- add_program(weights, ifn_genes, 5.5 * exp(ifn_score[i] / 2))
  if (cell_cycle_phase[i] == "S") {
    weights <- add_program(weights, s_genes, 10 * cell_cycle_score[i])
  }
  if (cell_cycle_phase[i] == "G2M") {
    weights <- add_program(weights, g2m_genes, 10 * cell_cycle_score[i])
  }
  weights <- add_program(weights, stress_genes, 3.5 * stress_score[i])
  weights[mitochondrial_genes] <- weights[mitochondrial_genes] *
    exp(1.05 * mitochondrial_technical_score[i])

  for (program in seq_len(n_weak_programs)) {
    weights[weak_program_genes[[program]]] <-
      weights[weak_program_genes[[program]]] * exp(0.42 * weak_scores[i, program])
  }

  donor_genes <- generic_genes[220L + (match(donor_id[i], c("D1", "D2", "D3")) - 1L) * 15L + seq_len(15L)]
  weights[donor_genes] <- weights[donor_genes] * 2.1

  library_size <- as.integer(round(rlnorm(1, log(5200), 0.24)))
  counts[, i] <- as.integer(rmultinom(1, library_size, weights / sum(weights)))
}

storage.mode(counts) <- "integer"
sparse_counts <- as(counts, "dgCMatrix")

metadata <- data.frame(
  cell_id = cell_id,
  sample_id = sample_id,
  donor_id = donor_id,
  condition = condition,
  truth_cell_type = truth_cell_type,
  truth_t_state = truth_t_state,
  cell_cycle_phase = cell_cycle_phase,
  cell_cycle_score = round(cell_cycle_score, 6),
  ifn_score = round(ifn_score, 6),
  stress_score = round(stress_score, 6),
  mitochondrial_technical_score = round(mitochondrial_technical_score, 6),
  stringsAsFactors = FALSE
)

expected_files <- c(
  counts = file.path(data_dir, "synthetic_representation_counts.rds"),
  metadata = file.path(data_dir, "synthetic_representation_metadata.csv")
)
manifest_path <- file.path(data_dir, "representation_data_manifest.csv")

write_outputs <- function() {
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
    rows = c(nrow(sparse_counts), nrow(metadata)),
    columns = c(ncol(sparse_counts), ncol(metadata)),
    stringsAsFactors = FALSE
  )
  write.csv(manifest, manifest_path, row.names = FALSE, na = "")
}

check_outputs <- function() {
  missing <- unname(expected_files)[!file.exists(unname(expected_files))]
  if (length(missing)) {
    stop("Missing representation data files: ", paste(missing, collapse = ", "))
  }
  saved_counts <- readRDS(expected_files[["counts"]])
  saved_metadata <- read.csv(expected_files[["metadata"]], check.names = FALSE)

  stopifnot(
    inherits(saved_counts, "dgCMatrix"),
    identical(dim(saved_counts), c(6000L, 480L)),
    identical(rownames(saved_counts), genes),
    identical(colnames(saved_counts), saved_metadata$cell_id),
    identical(colnames(saved_metadata), colnames(metadata)),
    all(saved_counts@x >= 0),
    all(saved_counts@x == floor(saved_counts@x)),
    all(c("CD3D", "MS4A1", "NKG7", "LST1", "CCR7", "GNLY") %in% rownames(saved_counts)),
    setequal(unique(saved_metadata$truth_cell_type), names(cell_type_counts)),
    setequal(unique(saved_metadata$cell_cycle_phase), c("G1", "S", "G2M"))
  )

  if (!file.exists(manifest_path)) {
    stop("Missing representation manifest: ", manifest_path)
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
    "Representation data check passed: ",
    nrow(saved_counts), " genes × ", ncol(saved_counts), " cells."
  )
}

if (check_only) {
  check_outputs()
} else {
  write_outputs()
  check_outputs()
  message("Synthetic representation data generated deterministically.")
}
