#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))

# Clean up leftover download artifacts when the script exits
cleanup_tmp <- function() {
  tmp <- tempdir()
  tmp_contents <- list.files(tmp, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (length(tmp_contents)) {
    unlink(tmp_contents, recursive = TRUE, force = TRUE)
  }
}
on.exit(cleanup_tmp(), add = TRUE)

# CRAN packages handled via devcontainer r-packages feature
cran_packages <- character(0)
bioc_packages <- c(
  "airway", "SummarizedExperiment", "DESeq2", "vsn", "apeglm",
  "AnnotationDbi", "org.Hs.eg.db", "clusterProfiler", "enrichplot",
  "DOSE", "IsoformSwitchAnalyzeR"
)

installed <- rownames(installed.packages(lib.loc = .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", update = FALSE, ask = FALSE)
}

needed_bioc <- setdiff(bioc_packages, installed)
if (length(needed_bioc)) {
  BiocManager::install(needed_bioc, update = FALSE, ask = FALSE)
}

needed_cran <- setdiff(cran_packages, installed)
if (length(needed_cran)) {
  install.packages(needed_cran, update = FALSE, ask = FALSE)
}
