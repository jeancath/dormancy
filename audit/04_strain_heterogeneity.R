#!/usr/bin/env Rscript
# Audit extension: direct test of the "MAG = population consensus" explanation using CheckM1 strain
# heterogeneity (GTDB r232 metadata; share of contamination from closely related strains, AAI > 90 %).
# Predictions if strain mixing drives the MAG slope:
#   (1) strain heterogeneity (SH) rises with P0 in MAGs;
#   (2) the MAG P0 slope (DORM_z, hib_extra) is concentrated in MAGs with SH > 0 (or rises with SH).
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table)})
tab <- as.data.frame(fread("data_small/final/analysis_table.tsv.gz")); tr <- read.tree("data_small/analysis_tree.nwk")
sh <- as.data.frame(fread("audit/strain_heterogeneity.tsv.gz")); tab <- merge(tab, sh, by = "accession")
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v)) / sd(v))
tab$P0 <- z(tab$log_T_sd_total, tab$domain); tab$lgs <- log10(tab$genome_size)
tab$is_mag <- as.integer(as.logical(tab$is_mag)); tab$sh <- tab$checkm_strain_heterogeneity
tab$sh_any <- as.integer(tab$sh > 0); tab$cc <- tab$checkm_contamination
rownames(tab) <- tab$accession
fit <- function(d, f, term, fam = "g") {
  ph <- keep.tip(tr, d$accession); d <- d[ph$tip.label, ]
  if (length(unique(d$domain)) < 2) f <- sub("\\+ domain", "", f)
  m <- if (fam == "b") phyloglm(as.formula(f), data = d, phy = ph, method = "logistic_MPLE", btol = 30) else
    phylolm(as.formula(f), data = d, phy = ph, model = "lambda")
  co <- summary(m)$coefficients; data.frame(n = nrow(d), estimate = co[term, 1], se = co[term, 2], p = co[term, ncol(co)])
}
M <- tab[tab$is_mag == 1, ]
out <- list()
cat("MAGs with SH>0:", sum(M$sh_any), "of", nrow(M), "\n")
out[[1]] <- cbind(test = "SH_any ~ P0 (MAGs, logit)", fit(M, "sh_any ~ P0 + domain + lgs + checkm2_completeness + cc", "P0", "b"))
out[[2]] <- cbind(test = "log1p(SH) ~ P0 (MAGs)", fit(M, "log1p(sh) ~ P0 + domain + lgs + checkm2_completeness + cc", "P0"))
for (resp in c("DORM_z", "log1p(hib_extra)", "log1p(n_hib)")) {
  for (g in c(0, 1)) out[[length(out) + 1]] <- cbind(test = sprintf("%s ~ P0, MAGs SH %s", resp, ifelse(g, "> 0", "= 0")),
      fit(M[M$sh_any == g, ], paste(resp, "~ P0 + domain + lgs + checkm2_completeness"), "P0"))
  M$P0_x_sh <- M$P0 * M$sh_any
  out[[length(out) + 1]] <- cbind(test = sprintf("%s: P0 x (SH>0) interaction, MAGs", resp),
      fit(M, paste(resp, "~ P0 + sh_any + P0_x_sh + domain + lgs + checkm2_completeness"), "P0_x_sh"))
  out[[length(out) + 1]] <- cbind(test = sprintf("%s ~ P0 + SH covariate (MAGs)", resp),
      fit(M, paste(resp, "~ P0 + log1p(sh) + cc + domain + lgs + checkm2_completeness"), "P0"))
}
r <- do.call(rbind, out); print(r, digits = 3); fwrite(r, "audit/strain_heterogeneity_tests.tsv", sep = "\t")
