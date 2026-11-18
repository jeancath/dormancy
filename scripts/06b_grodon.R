#!/usr/bin/env Rscript
# gRodon2 predicted minimal doubling time per genome (mode "full" = MILC + consistency + GC3 terms)
# usage: Rscript 06b_grodon.R CDS_DIR OUT.tsv [ncores]
suppressPackageStartupMessages({library(gRodon); library(Biostrings); library(parallel)})
a <- commandArgs(TRUE); cds_dir <- a[1]; out <- a[2]; nc <- ifelse(length(a) > 2, as.integer(a[3]), 4)
files <- list.files(cds_dir, pattern = "\\.cds\\.fna\\.gz$", full.names = TRUE)
done <- if (file.exists(out)) read.delim(out)$accession else character(0)
files <- files[!sub("\\.cds\\.fna\\.gz$", "", basename(files)) %in% done]
f1 <- function(f) {
  acc <- sub("\\.cds\\.fna\\.gz$", "", basename(f))
  r <- tryCatch({
    g <- readDNAStringSet(f)
    he <- grepl(" HE$", names(g))
    if (sum(he) < 10) stop("fewer than 10 HE genes")
    p <- predictGrowth(g, he, mode = "full")
    data.frame(accession = acc, d_hours = p$d, lower = p$LowerCI, upper = p$UpperCI,
               CUBHE = p$CUBHE, ConsistencyHE = p$ConsistencyHE, CUB = p$CUB, n_HE = sum(he), err = "")
  }, error = function(e) data.frame(accession = acc, d_hours = NA, lower = NA, upper = NA, CUBHE = NA,
                                     ConsistencyHE = NA, CUB = NA, n_HE = NA, err = conditionMessage(e)))
  r
}
chunks <- split(files, ceiling(seq_along(files) / 200))
for (ch in chunks) {
  res <- do.call(rbind, mclapply(ch, f1, mc.cores = nc))
  write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE,
              append = file.exists(out), col.names = !file.exists(out))
  cat(nrow(res), "\n")
}
