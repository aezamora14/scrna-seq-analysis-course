# Execute the ordinary R fences in each representation lesson in order.
for (lesson in sort(Sys.glob("modules/0[4-7]*/index.qmd"))) {
  lines <- readLines(lesson, warn = FALSE)
  starts <- which(lines == "```r")
  blocks <- lapply(starts, function(start) {
    end <- which(seq_along(lines) > start & lines == "```")[1]
    lines[seq.int(start + 1, end - 1)]
  })
  pdf(tempfile(fileext = ".pdf"))
  eval(parse(text = unlist(blocks)), envir = new.env(parent = globalenv()))
  dev.off()
  message("Lesson code passed: ", lesson)
}
