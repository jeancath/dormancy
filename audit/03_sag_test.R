#!/usr/bin/env Rscript
# Audit: do single-amplified genomes (SAGs; one cell, no population consensus) show the positive
# MAG-type slope? Marine only (all 92 SAGs are marine). phylolm lambda, same covariates as E_MAGs_marine.
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table)})
tab <- as.data.frame(fread("data_small/final/analysis_table.tsv.gz")); tr <- read.tree("data_small/analysis_tree.nwk")
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v)) / sd(v))
tab$P0 <- z(tab$log_T_sd_total, tab$domain); tab$lgs <- log10(tab$genome_size)
tab$src <- ifelse(tab$ncbi_genome_category == "derived from single cell", "SAG",
                  ifelse(tab$ncbi_genome_category == "derived from metagenome", "MAG", "isolate"))
rownames(tab) <- tab$accession; m <- tab[tab$domain == "marine", ]
out <- list()
for (s in c("SAG", "MAG", "isolate")) for (resp in c("DORM_z", "log1p(n_hib)", "log1p(hib_extra)")) {
  d <- m[m$src == s, ]; ph <- keep.tip(tr, d$accession); d <- d[ph$tip.label, ]
  co <- summary(phylolm(as.formula(paste(resp, "~ P0 + lgs + checkm2_completeness")), data = d, phy = ph, model = "lambda"))$coefficients
  out[[length(out) + 1]] <- data.frame(source = s, response = resp, n = nrow(d), sd_P0 = sd(d$P0), estimate = co["P0", 1], se = co["P0", 2], p = co["P0", 4])
}
r <- do.call(rbind, out); print(r); fwrite(r, "audit/sag_test.tsv", sep = "\t")
