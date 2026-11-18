#!/usr/bin/env Rscript
# Within-species DORM_z ~ P0 x is_mag (species FE, site-clustered HC1) re-fitted at completeness
# thresholds 0/90/95/98 %. Species left with < 2 genomes are dropped at each threshold.
# Added after the MQ extension showed a completeness-driven negative MAG slope (missing single-copy
# stringent-response genes in incomplete MAGs).
# usage (repo root): Rscript replication_within_species/10_completeness_strata.R TABLE.tsv.gz OUT.tsv label
suppressPackageStartupMessages({library(data.table); library(sandwich)})
a <- commandArgs(TRUE); d <- fread(a[1]); d[, comp := checkm2_completeness]; out <- list()
for (y in c("DORM_z", "log1p(hib_extra)", "log1p(n_sr)")) for (th in c(0, 90, 95, 98)) {
  dd <- d[comp >= th]; dd <- dd[dd[, .N, by = species][N >= 2], on = "species"]; dd[, species := factor(species)]
  m <- lm(as.formula(paste(y, "~ P0 * is_mag + lgs + comp + species")), dd); V <- vcovCL(m, cluster = dd$site, type = "HC1"); b <- coef(m)
  g <- c("P0", "P0:is_mag"); est <- c(b["P0"], b["P0:is_mag"], sum(b[g]))
  se <- c(sqrt(V["P0", "P0"]), sqrt(V["P0:is_mag", "P0:is_mag"]), sqrt(sum(V[g, g])))
  out[[length(out) + 1]] <- data.table(sample = a[3], model = paste0("b_speciesFE_clsite_comp", th), response = y,
    term = c("P0 (isolates)", "P0:is_mag", "P0 (MAGs)"), estimate = est, se = se, p = 2 * pnorm(-abs(est / se)),
    n_genomes = nrow(dd), n_species = nlevels(dd$species), note = sprintf("completeness >= %d; MAG genomes = %d", th, sum(dd$is_mag)))
}
o <- rbindlist(out); fwrite(o, a[2], sep = "\t"); print(o[, .(model, response, term, round(estimate, 4), round(se, 4), n_genomes)])
