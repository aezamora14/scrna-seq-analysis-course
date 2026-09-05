# Execute all lab code and verify numerical invariants, including degenerate ARI.
source("scripts/course_helpers.R")
stopifnot(adjusted_rand_index(rep(1, 5), rep(2, 5)) == 1,
          adjusted_rand_index(1:5, 5:1) == 1,
          adjusted_rand_index(c(1, 1, 2, 2), c(2, 2, 1, 1)) == 1)
labs <- c("04-variable-features", "05-scaling-regression", "06-pca", "07-choosing-dimensions")
for (lab in labs) {
  pdf(tempfile(fileext = ".pdf"))
  sys.source(paste0("labs/", lab, "-lab.R"), envir = new.env(parent = globalenv()))
  dev.off()
}
message("Modules 4–7 smoke test passed.")
