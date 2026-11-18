#!/usr/bin/env Rscript
# Exploratory (not pre-registered): dissect the MAG vs isolate heterogeneity of the H1 slope.
#  (1) per-module slopes (log1p counts) in MAGs and in isolates
#  (2) leave-one-site-out for the 15 largest MAG sites (influence of single studies)
#  (4) negative controls: ribosomal-protein HE genes, long RSH, coding density
#  (3) same-clade comparison: families represented by both MAGs and isolates, with family fixed effects
# usage: Rscript 18_mag_isolate.R TABLE.tsv.gz TREE.nwk OUTPREFIX
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table)})
a <- commandArgs(TRUE)
tab <- as.data.frame(fread(a[1])); tr <- read.tree(a[2]); outp <- a[3]
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
tab$P0 <- z(tab$log_T_sd_total, tab$domain)
tab$lgs <- log10(tab$genome_size); tab$is_mag <- as.integer(as.logical(tab$is_mag))
tab$site <- paste(tab$domain, round(tab$lat * 2) / 2, round(tab$lon * 2) / 2)
tab$family <- sub(".*;f__([^;]*);.*", "\\1", tab$gtdb_taxonomy)
rownames(tab) <- tab$accession
fit <- function(d, resp, covs, term = "P0") {
  if (length(unique(d$domain)) < 2) covs <- setdiff(covs, "domain")
  ph <- keep.tip(tr, d$accession); d <- d[ph$tip.label, ]
  f <- as.formula(paste(resp, "~", paste(c("P0", covs), collapse = " + ")))
  co <- summary(phylolm(f, data = d, phy = ph, model = "lambda"))$coefficients
  data.frame(n = nrow(d), beta = co[term, 1], se = co[term, 2], p = co[term, 4])
}
cv <- c("domain", "lgs", "checkm2_completeness")
out <- list()
for (g in c(1, 0)) for (resp in c("DORM_z", "log1p(spo_strict)", "log1p(n_rpf)", "log1p(n_hib)", "log1p(n_sr)", "log1p(ta_pairs)")) {
  d <- tab[tab$is_mag == g, ]
  r <- cbind(analysis = "module_by_type", genome_type = ifelse(g == 1, "MAG", "isolate"), response = resp, dropped = "", fit(d, resp, cv))
  print(r); out[[length(out) + 1]] <- r
}
mags <- tab[tab$is_mag == 1, ]
top <- names(sort(table(mags$site), decreasing = TRUE))[1:15]
for (s in top) {
  d <- mags[mags$site != s, ]
  r <- cbind(analysis = "MAG_leave_one_site_out", genome_type = "MAG", response = "DORM_z",
             dropped = paste0(s, " (n=", sum(mags$site == s), ")"), fit(d, "DORM_z", cv))
  print(r); out[[length(out) + 1]] <- r
}
d <- mags[!mags$site %in% top, ]
r <- cbind(analysis = "MAG_drop_top15_sites", genome_type = "MAG", response = "DORM_z", dropped = "top 15 sites", fit(d, "DORM_z", cv))
print(r); out[[length(out) + 1]] <- r
# families with >= 5 MAGs and >= 5 isolates: interaction P0 x is_mag within shared families
fam_ok <- names(which(tapply(tab$is_mag, tab$family, function(v) sum(v == 1) >= 5 & sum(v == 0) >= 5)))
d <- tab[tab$family %in% fam_ok, ]
d$P0_x_mag <- d$P0 * d$is_mag
r <- cbind(analysis = "shared_families_interaction", genome_type = "both", response = "DORM_z",
           dropped = paste(length(fam_ok), "families"), fit(d, "DORM_z", c(cv, "is_mag", "P0_x_mag"), term = "P0_x_mag"))
print(r); out[[length(out) + 1]] <- r
# negative controls: recovery of near-universal genes should NOT track P0 (tests a MAG gene-recovery artifact)
for (g in c(1, 0)) for (resp in c("log(n_HE)", "log1p(n_rsh_long)", "coding_density")) {
  d <- tab[tab$is_mag == g & !is.na(suppressWarnings(as.numeric(tab$n_HE))), ]
  d$n_HE <- as.numeric(d$n_HE)
  r <- tryCatch(cbind(analysis = "negative_control", genome_type = ifelse(g == 1, "MAG", "isolate"), response = resp,
                      dropped = "", fit(d, resp, cv)), error = function(e) NULL)
  print(r); out[[length(out) + 1]] <- r
}
fwrite(do.call(rbind, out), paste0(outp, "_mag_isolate.tsv"), sep = "\t")
