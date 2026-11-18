#!/usr/bin/env Rscript
# Control: does the SAME species subset show the MAG/isolate reversal BETWEEN species?
# Species means of single-type species (all MAG or all isolate), y ~ P0 * is_mag + lgs + comp + domain,
# HC1 SEs, and with a phylum fixed effect (crude phylogeny control). If the reversal appears between
# species here but not within, the within-species null is not a sample artefact.
# usage (repo root): Rscript replication_within_species/09_between_on_subset.R TABLE.tsv.gz OUT.tsv [label]
suppressPackageStartupMessages({library(data.table); library(sandwich); library(lmtest)})
a <- commandArgs(TRUE); lab <- if (length(a) > 2) a[3] else "main"
d <- fread(a[1]); d[, comp := checkm2_completeness]
s <- d[, .(DORM_z = mean(DORM_z), hib_extra = mean(log1p(hib_extra)), P0 = mean(P0), lgs = mean(lgs), comp = mean(comp),
           mag = mean(is_mag), domain = domain[1], phylum = phylum[1], n = .N), by = species][mag %in% c(0, 1)]
s[, is_mag := as.integer(mag)]
out <- list()
for (y in c("DORM_z", "hib_extra")) for (ph in c(FALSE, TRUE)) {
  f <- as.formula(paste(y, "~ P0 * is_mag + lgs + comp + domain", if (ph) "+ phylum" else ""))
  m <- lm(f, data = s); V <- vcovHC(m, type = "HC1"); co <- coeftest(m, vcov. = V)
  g <- c(P0 = 1, `P0:is_mag` = 1); e <- sum(coef(m)[names(g)]); se <- sqrt(drop(t(g) %*% V[names(g), names(g)] %*% g))
  mk <- function(t, e, se, p) data.table(sample = lab, model = paste0("between_species_means", if (ph) "_phylumFE" else ""),
    response = y, term = t, estimate = e, se = se, p = p, n_genomes = sum(s$n), n_species = nrow(s))
  out[[length(out) + 1]] <- rbind(mk("P0 (isolates)", co["P0", 1], co["P0", 2], co["P0", 4]),
                                  mk("P0:is_mag", co["P0:is_mag", 1], co["P0:is_mag", 2], co["P0:is_mag", 4]),
                                  mk("P0 (MAGs)", e, se, 2 * pnorm(-abs(e / se))))
}
o <- rbindlist(out); fwrite(o, a[2], sep = "\t"); print(o[, 2:9], digits = 3)
