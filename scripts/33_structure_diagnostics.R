#!/usr/bin/env Rscript
# Added 18 Nov 2026: compare the dependence structure of the real data with the simulated data of DGM B
# (scripts/31) and DGM D (scripts/31c) on several metrics, not only the between-cell share of residual variance.
# Every data set (real, 20 draws of B, 20 draws of D, slope 0) is analysed with the same model,
#   y ~ P0 + domain + lgs + checkm2_completeness + is_mag (Pagel's lambda),
# and the following are computed on its residuals:
#   lambda_hat; between-cell share and between-order share of residual variance (one-way ANOVA estimator);
#   spatial autocorrelation between DIFFERENT 0.5-degree cells: Pearson correlation of cell-mean residuals over cell
#   pairs whose centres are 0-1, 1-3 and 3-10 degrees apart (the jackknife assumes cells independent);
#   correlation of residuals of species in the SAME order but DIFFERENT cells (cross-cell phylogenetic dependence);
#   residual skewness and excess kurtosis.
# usage: Rscript 33_structure_diagnostics.R TABLE.tsv.gz TREE.nwk OUT.tsv [ndraws=20] [cores=4]
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table); library(parallel)})
a <- commandArgs(TRUE); ND <- if (length(a) > 3) as.integer(a[4]) else 20; cores <- if (length(a) > 4) as.integer(a[5]) else 4
tab <- as.data.frame(fread(a[1])); tr <- read.tree(a[2])
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
tab$P0 <- z(tab$log_T_sd_total, tab$domain); tab$lgs <- log10(tab$genome_size)
tab$is_mag <- as.integer(as.logical(tab$is_mag)); rownames(tab) <- tab$accession
ph <- keep.tip(tr, tab$accession); tab <- tab[ph$tip.label, ]
ph$edge.length[ph$edge.length <= 0] <- 1e-8
tab$cell <- paste(tab$domain, floor(tab$lat / 0.5), floor(tab$lon / 0.5))
tab$site <- paste(tab$domain, round(tab$lat * 2) / 2, round(tab$lon * 2) / 2)
tab$order <- sub(".*;o__([^;]*);.*", "\\1", tab$gtdb_taxonomy)
cells <- unique(tab$cell); cidx <- match(tab$cell, cells); sites <- unique(tab$site); sidx <- match(tab$site, sites)
oidx <- match(tab$order, unique(tab$order)); nord <- length(unique(tab$order))
form <- y ~ P0 + domain + lgs + checkm2_completeness + is_mag
share <- function(r, g) { n_i <- tapply(r, g, length); m_i <- tapply(r, g, mean); k <- length(n_i); N <- length(r)
  msb <- sum(n_i * (m_i - mean(r))^2) / (k - 1); msw <- sum((r - m_i[g])^2) / (N - k); n0 <- (N - sum(n_i^2) / N) / (k - 1)
  sb <- max(0, (msb - msw) / n0); sb / (sb + msw) }
# cell centres and pair distances (degrees, equirectangular with cos(lat) correction)
cc <- aggregate(cbind(lat, lon) ~ cell, data = tab, FUN = mean); cc <- cc[match(cells, cc$cell), ]
pr <- t(combn(length(cells), 2))
dist <- sqrt((cc$lat[pr[, 1]] - cc$lat[pr[, 2]])^2 + ((cc$lon[pr[, 1]] - cc$lon[pr[, 2]]) * cos((cc$lat[pr[, 1]] + cc$lat[pr[, 2]]) / 2 * pi / 180))^2)
bins <- list(`0-1` = which(dist <= 1), `1-3` = which(dist > 1 & dist <= 3), `3-10` = which(dist > 3 & dist <= 10))
# same order, different cell: sample 200,000 random pairs
set.seed(1); byo <- split(seq_len(nrow(tab)), tab$order)
sp <- do.call(rbind, lapply(byo[lengths(byo) >= 2], function(ix) { m <- min(length(ix) * 3, 2000); cbind(sample(ix, m, TRUE), sample(ix, m, TRUE)) }))
sp <- sp[sp[, 1] != sp[, 2] & tab$cell[sp[, 1]] != tab$cell[sp[, 2]], ]
metrics <- function(y, label, draw) {
  d <- tab; d$y <- y
  f <- phylolm(form, data = d, phy = ph, model = "lambda"); r <- unname(f$residuals)
  cm <- tapply(r, factor(tab$cell, levels = cells), mean)
  ac <- sapply(bins, function(ix) cor(cm[pr[ix, 1]], cm[pr[ix, 2]]))
  m <- mean(r); s <- sd(r)
  data.frame(data = label, draw = draw, lambda_hat = f$optpar, between_cell_share = share(r, tab$cell),
             between_order_share = share(r, tab$order), autocor_cells_0_1deg = ac[1], autocor_cells_1_3deg = ac[2],
             autocor_cells_3_10deg = ac[3], cor_same_order_diff_cell = cor(r[sp[, 1]], r[sp[, 2]]),
             skewness = mean((r - m)^3) / s^3, excess_kurtosis = mean((r - m)^4) / s^4 - 3)
}
out <- list(metrics(tab$DORM_z, "real", 0))
print(out[[1]])
# DGM B
mc <- fread("results/siteMCMC_mcmcglmm_site.tsv")
out <- c(out, mclapply(seq_len(ND), function(i) { set.seed(9e6 + i)
  e <- rTrait(n = 1, ph, model = "BM", parameters = list(sigma2 = mc$var_animal))[tab$accession] +
       rnorm(length(sites), 0, sqrt(mc$var_site))[sidx] + rnorm(nrow(tab), 0, sqrt(mc$var_resid))
  metrics(e, "DGM B", i) }, mc.cores = cores))
# DGM D (slope 0)
full <- phylolm(DORM_z ~ P0 + domain + lgs + checkm2_completeness + is_mag, data = tab, phy = ph, model = "lambda")
lam <- full$optpar; s2 <- full$sigma2; r <- full$residuals
C <- vcv(ph)[tab$accession, tab$accession]; dg <- diag(C)
eps <- s2 * (1 - lam) * dg * solve(s2 * (lam * C + diag((1 - lam) * dg)), r); u <- r - eps; rm(C); gc()
fit0 <- full$fitted.values - coef(full)["P0"] * tab$P0
out <- c(out, mclapply(seq_len(ND), function(i) { set.seed(9.5e6 + i)
  wo <- sample(c(-1, 1), nord, TRUE)[oidx]; wc <- sample(c(-1, 1), length(cells), TRUE)[cidx]
  metrics(fit0 + wo * u + wc * eps, "DGM D", i) }, mc.cores = cores))
res <- rbindlist(out); fwrite(res, a[3], sep = "\t")
sm <- res[, lapply(.SD, function(v) sprintf("%.3f [%.3f, %.3f]", mean(v), min(v), max(v))), by = data, .SDcols = 3:11]
fwrite(sm, sub("\\.tsv$", "_summary.tsv", a[3]), sep = "\t"); print(t(sm))
