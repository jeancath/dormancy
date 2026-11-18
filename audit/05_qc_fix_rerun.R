#!/usr/bin/env Rscript
# Audit: re-fit H1 (direct, total) and the MAG x P0 interaction after applying the pre-registered QC
# rule correctly (exclude genomes with n_HE NA, i.e. < 10 HE genes in gRodon; code_review.md M2).
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table)})
tab <- as.data.frame(fread("data_small/final/analysis_table.tsv.gz")); tr <- read.tree("data_small/analysis_tree.nwk")
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v)) / sd(v))
tab$P0 <- z(tab$log_T_sd_total, tab$domain); tab$lgs <- log10(tab$genome_size); tab$is_mag <- as.integer(as.logical(tab$is_mag))
rownames(tab) <- tab$accession
out <- list()
for (lab in c("as_published", "qc_fixed")) {
  d <- if (lab == "qc_fixed") tab[!is.na(suppressWarnings(as.numeric(tab$n_HE))), ] else tab
  ph <- keep.tip(tr, d$accession); d <- d[ph$tip.label, ]
  for (m in list(c("H1_direct", "DORM_z ~ P0 + domain + lgs + checkm2_completeness + is_mag", "P0"),
                 c("H1_total", "DORM_z ~ P0 + domain + checkm2_completeness + is_mag", "P0"),
                 c("interaction", "DORM_z ~ P0 * is_mag + domain + lgs + checkm2_completeness", "P0:is_mag"))) {
    co <- summary(phylolm(as.formula(m[2]), data = d, phy = ph, model = "lambda"))$coefficients
    out[[length(out) + 1]] <- data.frame(set = lab, model = m[1], n = nrow(d), estimate = co[m[3], 1], se = co[m[3], 2],
      l95 = co[m[3], 1] - 1.96 * co[m[3], 2], u95 = co[m[3], 1] + 1.96 * co[m[3], 2], p = co[m[3], 4])
  }
}
r <- do.call(rbind, out); print(r, digits = 4); fwrite(r, "audit/qc_fix_rerun.tsv", sep = "\t")
