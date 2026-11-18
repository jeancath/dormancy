#!/usr/bin/env Rscript
# Exploratory: MAG x P0 interaction and per-genome-type slopes in a Gaussian MCMCglmm with
# phylogenetic (animal, sparse all-node inverse) AND site (0.5-degree cell) random effects.
# This is the model recommended for isolation-site predictors with many taxa per site.
# usage: Rscript 21_mcmc_site_interaction.R TABLE.tsv.gz TREE.nwk OUTPREFIX [nitt=43000]
suppressPackageStartupMessages({library(ape); library(MCMCglmm); library(data.table)})
a <- commandArgs(TRUE)
nitt <- if (length(a) > 3) as.integer(a[4]) else 43000
tab <- as.data.frame(fread(a[1])); tr <- read.tree(a[2]); outp <- a[3]
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
tab$P0 <- z(tab$log_T_sd_total, tab$domain)
tab$lgs <- log10(tab$genome_size); tab$is_mag <- as.integer(as.logical(tab$is_mag))
tab$site <- factor(paste(tab$domain, round(tab$lat * 2) / 2, round(tab$lon * 2) / 2))
tab$animal <- tab$accession
set.seed(20261117)
pr <- list(R = list(V = 1, nu = 0.002),
           G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                    G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
run <- function(d, f, label, terms) {
  ph <- keep.tip(tr, d$accession); ph$node.label <- NULL
  ph$edge.length[ph$edge.length <= 0] <- 1e-8   # GTDB branch lengths kept (audit M1); no Grafen rescaling
  Ainv <- inverseA(ph, nodes = "ALL", scale = FALSE)$Ainv
  d$site <- droplevels(d$site)
  m <- MCMCglmm(f, random = ~animal + site, ginverse = list(animal = Ainv), data = d, prior = pr,
                nitt = nitt, burnin = 3000, thin = 40, verbose = FALSE)
  s <- summary(m)$solutions; vc <- colMeans(m$VCV)
  do.call(rbind, lapply(terms, function(t) data.frame(model = label, term = t, n = nrow(d), n_sites = nlevels(d$site),
    post_mean = s[t, 1], l95 = s[t, 2], u95 = s[t, 3], pMCMC = s[t, 5], ess = effectiveSize(m$Sol[, t]),
    var_animal = vc["animal"], var_site = vc["site"], var_resid = vc["units"])))
}
out <- list()
out[[1]] <- run(tab, DORM_z ~ P0 * is_mag + domain + lgs + checkm2_completeness, "all_interaction", c("P0", "P0:is_mag"))
print(out[[1]]); fwrite(do.call(rbind, out), paste0(outp, "_mcmc_site_interaction.tsv"), sep = "\t")
out[[2]] <- run(tab[tab$is_mag == 1, ], DORM_z ~ P0 + domain + lgs + checkm2_completeness, "MAGs", "P0")
print(out[[2]]); fwrite(do.call(rbind, out), paste0(outp, "_mcmc_site_interaction.tsv"), sep = "\t")
out[[3]] <- run(tab[tab$is_mag == 0, ], DORM_z ~ P0 + domain + lgs + checkm2_completeness, "isolates", "P0")
print(out[[3]]); fwrite(do.call(rbind, out), paste0(outp, "_mcmc_site_interaction.tsv"), sep = "\t")
