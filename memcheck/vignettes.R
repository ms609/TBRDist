# Code to be run with
#   R -d "valgrind --tool=memcheck --leak-check=full --error-exitcode=1" --vanilla < memcheck/vignettes.R
library("TBRDist")

vigs <- vignette(package = "TBRDist")$results
if (nrow(vigs) == 0L) {
  cat("No vignettes found.\n")
  quit(status = 0)
}

cat("Running", nrow(vigs), "vignettes\n")

for (i in seq_len(nrow(vigs))) {
  topic <- vigs[i, "Item"]
  cat("\n>>> Vignette:", topic, "\n")
  v <- vignette(topic, "TBRDist")
  src <- file.path(v$Dir, "doc", v$R)
  if (!file.exists(src)) {
    cat("No tangled R file for vignette:", topic, "\n")
    next
  }
  tryCatch(
    {
      sys.source(src, envir = globalenv())
      cat("\U2713 Success:", topic, "\n")
    },
    error = function(e) {
      cat("\U2718 Error in vignette:", topic, "\n", conditionMessage(e), "\n")
    }
  )
}
cat("\nFinished running vignettes.\n")
