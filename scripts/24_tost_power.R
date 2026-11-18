#!/usr/bin/env Rscript
# TOST operating characteristics for H1 (added 18 Nov 2026).
# For true slopes beta in {-0.015, 0, 0.005, 0.010, 0.015, 0.020}, simulate DORM_z on the confirmatory tree:
#   y = beta * P0 + phylogenetic error (Pagel lambda and sigma2 from the confirmatory H1 phylolm fit)
#       + site effect (0.5-degree cell; SD from the MCMCglmm phylogeny + site model, GTDB branch lengths),
# fit phylolm(y ~ P0 + domain, lambda) and apply TOST against +-0.015 (90 % CI inside the bounds).
# Reports P(declare equivalence) [power when |beta| < 0.015; false-equivalence rate at |beta| = 0.015],
# P(beta > 0, p < 0.05), and the mean / SD of the estimate.
# usage: Rscript 24_tost_power.R TABLE.tsv.gz TREE.nwk OUT.tsv [nsim=200] [cores=3]
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table); library(parallel)})
a <- commandArgs(TRUE)
nsim <- if (length(a) > 3) as.integer(a[4]) else 200
cores <- if (length(a) > 4) as.integer(a[5]) else 3
set.seed(20261118)
tab <- as.data.frame(fread(a[1])); tr <- read.tree(a[2])
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
tab$P0 <- z(tab$log_T_sd_total, tab$domain); tab$lgs <- log10(tab$genome_size)
tab$is_mag <- as.integer(as.logical(tab$is_mag))
tab$site <- paste(tab$domain, round(tab$lat * 2) / 2, round(tab$lon * 2) / 2)
rownames(tab) <- tab$accession
ph <- keep.tip(tr, tab$accession); tab <- tab[ph$tip.label, ]
fit <- phylolm(DORM_z ~ P0 + domain + lgs + checkm2_completeness + is_mag, data = tab, phy = ph, model = "lambda")
lam <- fit$optpar; s2 <- fit$sigma2
site_sd <- sqrt(fread("results/siteMCMC_mcmcglmm_site.tsv")$var_site[1])
cat("lambda", lam, "sigma2", s2, "site SD", site_sd, "n", nrow(tab), "sites", length(unique(tab$site)), "\n")
ph_l <- transf.branch.lengths(ph, model = "lambda", parameters = list(lambda = lam))$tree
sites <- unique(tab$site); sidx <- match(tab$site, sites)
SESOI <- 0.015
one <- function(i, b) {
  e <- rTrait(n = 1, ph_l, model = "BM", parameters = list(sigma2 = s2))[tab$accession]
  u <- rnorm(length(sites), 0, site_sd)[sidx]
  d <- data.frame(y = b * tab$P0 + e + u, P0 = tab$P0, domain = tab$domain, row.names = tab$accession)
  co <- summary(phylolm(y ~ P0 + domain, data = d, phy = ph, model = "lambda"))$coefficients["P0", ]
  est <- co[1]; se <- co[2]
  c(est = est, se = se, p = co[4], equiv = as.numeric(est - 1.645 * se > -SESOI & est + 1.645 * se < SESOI))
}
betas <- c(-0.015, 0, 0.005, 0.010, 0.015, 0.020)
res <- do.call(rbind, lapply(betas, function(b) {
  r <- do.call(rbind, mclapply(seq_len(nsim), one, b = b, mc.cores = cores))
  out <- data.frame(true_beta = b, nsim = nsim, P_declare_equivalence = mean(r[, 4]),
                    P_detect_positive = mean(r[, 3] < 0.05 & r[, 1] > 0), mean_est = mean(r[, 1]),
                    sd_est = sd(r[, 1]), mean_se = mean(r[, 2]))
  print(out); out
}))
res$lambda <- lam; res$sigma2 <- s2; res$site_sd <- site_sd
fwrite(res, a[3], sep = "\t")
