# Module 1 lab: count matrices and measurement
# Run from the repository root.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260904)
dir.create("outputs", showWarnings = FALSE)

counts <- readRDS("data/synthetic_qc_counts.rds")
metadata <- read.csv("data/synthetic_qc_metadata.csv", na.strings = "")
rownames(metadata) <- metadata$cell_id

obj <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata,
  min.cells = 0,
  min.features = 0
)

# Part 1: matrix anatomy -----------------------------------------------------

print(class(counts))
print(dim(counts))
print(counts[1:8, 1:6])

nonzero_entries <- length(counts@x)
total_entries <- nrow(counts) * ncol(counts)
sparsity <- 1 - nonzero_entries / total_entries

print(c(
  nonzero_entries = nonzero_entries,
  total_entries = total_entries,
  sparsity = sparsity
))

print(c(
  sparse_bytes = as.numeric(object.size(counts)),
  dense_bytes = as.numeric(object.size(as.matrix(counts)))
))

# Interpretation prompt:
# Why do sparse and dense representations contain the same scientific values
# while requiring different amounts of memory?

# Part 2: one cell and one gene --------------------------------------------

cell_id <- "BC0001-1"
cell_counts <- counts[, cell_id, drop = FALSE]

print(sum(cell_counts))
print(sum(cell_counts > 0))
named_cell_counts <- setNames(as.numeric(cell_counts), rownames(cell_counts))
print(sort(named_cell_counts, decreasing = TRUE)[1:10])

gene_counts <- counts["NKG7", ]
print(summary(as.numeric(gene_counts)))
print(mean(gene_counts > 0))

# Interpretation prompt:
# What does the cell-oriented summary reveal that the gene-oriented summary
# does not, and vice versa?

# Part 3: derive Seurat metadata manually ----------------------------------

manual_library_size <- Matrix::colSums(counts)
manual_detected_features <- Matrix::colSums(counts > 0)

stopifnot(all(manual_library_size == obj$nCount_RNA))
stopifnot(all(manual_detected_features == obj$nFeature_RNA))

qc_summary <- data.frame(
  cell_id = Cells(obj),
  library_size = manual_library_size,
  detected_features = manual_detected_features,
  truth_qc_class = obj$truth_qc_class
)

print(aggregate(
  cbind(library_size, detected_features) ~ truth_qc_class,
  qc_summary,
  median
))

# Part 4: controlled depth pair --------------------------------------------

toy <- read.csv(
  "data/toy_counts.csv",
  row.names = 1,
  check.names = FALSE
)
toy <- as.matrix(toy)

toy_proportions <- sweep(toy, 2, colSums(toy), "/")

print(toy)
print(colSums(toy))
print(round(toy_proportions, 3))

# Interpretation prompts:
# 1. Which values differ between Cell_A and Cell_B?
# 2. Which proportions differ?
# 3. List at least three mechanisms that could produce a library-size change.

# Part 5: change only the sampling probability -----------------------------

candidate <- which(
  obj$truth_qc_class == "high_quality" &
    obj$nCount_RNA > median(obj$nCount_RNA)
)[1]
original <- as.integer(counts[, candidate])

probabilities <- c(0.8, 0.5, 0.35, 0.2, 0.1)
sampling_results <- lapply(
  probabilities,
  function(probability) {
    thinned <- rbinom(
      length(original),
      size = original,
      prob = probability
    )

    data.frame(
      probability = probability,
      library_size = sum(thinned),
      detected_genes = sum(thinned > 0),
      newly_zero = sum(original > 0 & thinned == 0)
    )
  }
)
sampling_results <- do.call(rbind, sampling_results)
print(sampling_results)

library_plot <- ggplot(
  sampling_results,
  aes(probability, library_size)
) +
  geom_line() +
  geom_point() +
  labs(x = "UMI retention probability", y = "Observed library size")

zero_plot <- ggplot(
  sampling_results,
  aes(probability, newly_zero)
) +
  geom_line() +
  geom_point() +
  labs(x = "UMI retention probability", y = "Newly observed zeros")

combined_plot <- library_plot + zero_plot
print(combined_plot)
ggsave(
  "outputs/module-01-undersampling.png",
  combined_plot,
  width = 9,
  height = 4.5,
  dpi = 150
)

# Final interpretation prompts:
# - What changed across the five simulations?
# - What was held constant?
# - Why were low-count genes most likely to disappear?
# - Could an analyst identify each new sampling zero without the original data?
# - How could the change alter a biological comparison?

message("Module 1 student lab executed successfully.")
