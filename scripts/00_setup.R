options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

renv::restore(prompt = FALSE)

message("The course R environment matches renv.lock.")
message("Next run: Rscript scripts/check_environment.R")
