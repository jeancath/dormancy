#!/usr/bin/env Rscript
# Exploratory (18 Nov 2026): does the MAG x P0 interaction survive adjustment for the
# taxonomic and environmental composition of the two genome sources?
# Every model is fitted twice: phylolm (Pagel's lambda, model SE) and the primary procedure P (lambda fixed at the
# model's estimate, delete-a-group jackknife over G = 20 random groups of 0.5-degree grid cells, t on 19 df).
# Composition terms are entered both as main effects and as interactions with P0, so that the MAG x P0 term is the
# difference in slopes that remains after letting the slope differ by habitat, phylum, latitude and mean temperature.
# usage: Rscript 28_interaction_composition.R TABLE.tsv.gz TREE.nwk OUT.tsv
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table)})
a <- commandArgs(TRUE); G <- 20; R <- if (length(a) > 3) as.integer(a[4]) else 1
suppressPackageStartupMessages(library(parallel))
# with R > 1 the jackknife variance is averaged over R independent random groupings (SE_avg), because the SE
# from a single grouping varies widely between groupings (scripts/30_grouping_sensitivity.R). Columns *_jackknife then
# refer to SE_avg; se_jackknife_first is the SE from the first grouping alone.
tab <- as.data.frame(fread(a[1])); tr <- read.tree(a[2])
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
tab$P0 <- z(tab$log_T_sd_total, tab$domain); tab$lgs <- log10(tab$genome_size)
tab$is_mag <- as.integer(as.logical(tab$is_mag)); rownames(tab) <- tab$accession
tab$hab <- ifelse(tab$habitat == "water_unspecified", "marine", tab$habitat)          # freshwater / soil / marine
tab$sed <- as.integer(as.logical(tab$is_sediment))
tab$Tmean_z <- z(ifelse(tab$domain == "marine", tab$sst_mean, tab$tair_mean), tab$domain)
tab$abslat_z <- as.numeric(scale(abs(tab$lat)))
both <- names(which(tapply(tab$is_mag, tab$phylum, function(v) sum(v == 1) >= 20 & sum(v == 0) >= 20)))
tab$phy5 <- ifelse(tab$phylum %in% both, tab$phylum, "other")
tab$stratum <- paste(tab$phy5, tab$hab)
tab$cell <- paste(tab$domain, floor(tab$lat / 0.5), floor(tab$lon / 0.5))
tab$MAGxP0 <- tab$is_mag * tab$P0
cat("phyla with >= 20 MAG- and >= 20 isolate-represented species:", paste(both, collapse = ", "), "\n")
base <- "P0 + is_mag + MAGxP0 + lgs + checkm2_completeness"
models <- list(
  M0_unadjusted       = list(f = paste(base, "+ domain"), sub = NULL),
  M1_habitat          = list(f = paste(base, "+ hab + hab:P0 + sed + sed:P0"), sub = NULL),
  M2_phylum           = list(f = paste(base, "+ domain + phy5 + phy5:P0"), sub = NULL),
  M3_lat_meanT        = list(f = paste(base, "+ domain + abslat_z + abslat_z:P0 + Tmean_z + Tmean_z:P0"), sub = NULL),
  M4_all              = list(f = paste(base, "+ hab + hab:P0 + sed + sed:P0 + phy5 + phy5:P0 + abslat_z + abslat_z:P0 + Tmean_z + Tmean_z:P0"), sub = NULL),
  M5_common_support   = list(f = paste(base, "+ stratum + stratum:P0 + abslat_z + abslat_z:P0 + Tmean_z + Tmean_z:P0"), sub = "support"))
# common support: phylum x habitat strata with >= 15 MAG- and >= 15 isolate-represented species
ok <- names(which(tapply(tab$is_mag, tab$stratum, function(v) sum(v == 1) >= 15 & sum(v == 0) >= 15)))
cat("common-support strata:", length(ok), "\n"); print(table(tab$stratum[tab$stratum %in% ok], tab$is_mag[tab$stratum %in% ok]))
set.seed(20261118)
out <- list()
for (m in names(models)) {
  d <- if (is.null(models[[m]]$sub)) tab else tab[tab$stratum %in% ok, ]
  ph <- keep.tip(tr, d$accession); d <- d[ph$tip.label, ]
  f <- as.formula(paste("DORM_z ~", models[[m]]$f))
  fl <- phylolm(f, data = d, phy = ph, model = "lambda"); co <- summary(fl)$coefficients
  phl <- transf.branch.lengths(ph, model = "lambda", parameters = list(lambda = fl$optpar))$tree
  cells <- unique(d$cell)
  ses <- unlist(mclapply(seq_len(R), function(r) {
    set.seed(20261118 + r); grp <- sample(rep_len(seq_len(G), length(cells))); names(grp) <- cells
    e <- sapply(seq_len(G), function(g) { keep <- grp[d$cell] != g; p2 <- keep.tip(phl, d$accession[keep])
      coef(phylolm(f, data = d[p2$tip.label, ], phy = p2, model = "BM"))["MAGxP0"] })
    sqrt((G - 1) / G * sum((e - mean(e))^2)) }, mc.cores = 4))
  jse <- sqrt(mean(ses^2))
  # slope in isolates = P0 (reference levels / covariates at 0); in MAGs = P0 + MAGxP0 (reported for M0 and M5 only,
  # where the reference levels are the whole sample or a balanced set of strata)
  r <- data.frame(model = m, n = nrow(d), n_mag = sum(d$is_mag), lambda = fl$optpar,
                  n_groupings = R, se_jackknife_first = ses[1],
                  interaction = co["MAGxP0", 1], se_model = co["MAGxP0", 2], p_model = co["MAGxP0", 4],
                  ci95_lo_model = co["MAGxP0", 1] - 1.96 * co["MAGxP0", 2], ci95_hi_model = co["MAGxP0", 1] + 1.96 * co["MAGxP0", 2],
                  se_jackknife = jse[1], p_jackknife = 2 * pt(abs(co["MAGxP0", 1] / jse[1]), G - 1, lower.tail = FALSE),
                  ci95_lo_jackknife = co["MAGxP0", 1] - qt(0.975, G - 1) * jse[1],
                  ci95_hi_jackknife = co["MAGxP0", 1] + qt(0.975, G - 1) * jse[1],
                  formula = models[[m]]$f)
  print(r[, 1:12]); out[[m]] <- r
}
fwrite(do.call(rbind, out), a[3], sep = "\t")
