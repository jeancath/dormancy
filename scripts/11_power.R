#!/usr/bin/env Rscript
# Phase 4: simulation power analysis for H1 on the confirmatory tree.
# Residual structure (Pagel lambda, sigma2) is taken from the pilot fit of H1; the predictor values are
# the real confirmatory P0 values (no response data from the confirmatory set are used).
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table); library(parallel)})
set.seed(20261117)
tr <- read.tree("data_small/analysis_tree.nwk")
pil <- as.data.frame(fread("/home/user/data/hits/pilot_table.tsv.gz"))
rownames(pil) <- pil$accession
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v)) / sd(v))
a <- as.data.frame(fread("data_small/analysis_set.tsv.gz"))
a$P0 <- z(a$log_T_sd_total, a$domain)          # standardised over the whole analysis set
pil$P0 <- a$P0[match(pil$accession, a$accession)]
pil$lgs <- log10(pil$genome_size); pil$is_mag <- as.integer(as.logical(pil$is_mag))
ph <- keep.tip(tr, pil$accession); pil <- pil[ph$tip.label, ]
fit <- phylolm(DORM_z ~ P0 + domain + lgs + checkm2_completeness + is_mag, data = pil, phy = ph, model = "lambda")
lam <- fit$optpar; s2 <- fit$sigma2
cat("pilot: beta =", coef(fit)["P0"], " se =", summary(fit)$coefficients["P0", 2], " lambda =", lam, " sigma2 =", s2,
    " resid SD of DORM_z (pilot) =", sd(residuals(fit)), " SD DORM_z =", sd(pil$DORM_z), "\n")
conf <- a[a$pilot == 0 & a$accession %in% tr$tip.label, ]
ct <- keep.tip(tr, conf$accession); conf <- conf[match(ct$tip.label, conf$accession), ]
rownames(conf) <- conf$accession
ct_l <- transf.branch.lengths(ct, model = "lambda", parameters = list(lambda = lam))$tree
betas <- c(0, 0.005, 0.01, 0.015, 0.02, 0.03)
nsim <- 100
one <- function(i, b) {
  e <- rTrait(n = 1, ct_l, model = "BM", parameters = list(sigma2 = s2))
  d <- data.frame(y = b * conf$P0 + e[conf$accession], P0 = conf$P0, domain = conf$domain, row.names = conf$accession)
  f <- phylolm(y ~ P0 + domain, data = d, phy = ct, model = "lambda")
  co <- summary(f)$coefficients["P0", ]
  c(est = co[1], p = co[4])
}
res <- do.call(rbind, lapply(betas, function(b) {
  r <- do.call(rbind, mclapply(seq_len(nsim), one, b = b, mc.cores = 3))
  data.frame(beta = b, n = nrow(conf), power = mean(r[, 2] < 0.05 & r[, 1] > 0), mean_est = mean(r[, 1]),
             sd_est = sd(r[, 1]))
}))
print(res)
fwrite(res, "results/power_simulation.tsv", sep = "\t")
