#!/usr/bin/env Rscript
# NOTE: this naive wild cluster bootstrap destroyed the phylogenetic correlation between cells (lambda-hat about
# 0.04 instead of 0.94) and was stopped after 80 runs; its output is results/calibration_wild_cluster_naive_*. It is
# replaced by scripts/31c_calibration_real_structure.R (DGM D).
# Added 18 Nov 2026: calibration of procedure P under a data-generating model that keeps the REAL within-cell
# dependence of the residuals ("DGM C", wild cluster bootstrap at the boundary).
# Why: the parametric DGMs of scripts/26 and 31 reproduce lambda-hat and the SD of the slope, but not the clustering of
# residuals within 0.5-degree cells (between-cell share of residual variance about 0.13-0.17 in simulations, 0.50 in the
# real data for the same model form). Here:
#   1. impose the boundary slope b on the real data: fit (DORM_z - b * P0) ~ domain + lgs + checkm2_completeness + is_mag
#      (Pagel's lambda), giving fitted values f and residuals r;
#   2. each simulated response is y* = b * P0 + f + w_c * r, with an independent Rademacher sign w_c (+1/-1) per
#      0.5-degree cell, so every cell keeps its real residual pattern and cells are independent;
#   3. procedure P is applied with lambda re-estimated, with the original single-grouping SE and with the SE averaged over
#      K = 5 random groupings (as in scripts/31_calibration_lambda_free.R).
# 500 simulations at each boundary; restart-safe.
# usage: Rscript 31b_calibration_wild_cluster.R TABLE.tsv.gz TREE.nwk OUT_PREFIX [cores=4] -> OUT_PREFIX_sims.tsv, _summary.tsv
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table); library(parallel)})
a <- commandArgs(TRUE); cores <- if (length(a) > 3) as.integer(a[4]) else 4
G <- 20; SESOI <- 0.015; NSIM <- 500; K <- 5
tab <- as.data.frame(fread(a[1])); tr <- read.tree(a[2])
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
tab$P0 <- z(tab$log_T_sd_total, tab$domain); tab$lgs <- log10(tab$genome_size)
tab$is_mag <- as.integer(as.logical(tab$is_mag)); rownames(tab) <- tab$accession
ph <- keep.tip(tr, tab$accession); tab <- tab[ph$tip.label, ]
tab$cell <- paste(tab$domain, floor(tab$lat / 0.5), floor(tab$lon / 0.5))
cells <- unique(tab$cell); cidx <- match(tab$cell, cells)
form <- y ~ P0 + domain + lgs + checkm2_completeness + is_mag
between_share <- function(r, cell) { n_i <- tapply(r, cell, length); m_i <- tapply(r, cell, mean); k <- length(n_i); N <- length(r)
  msb <- sum(n_i * (m_i - mean(r))^2) / (k - 1); msw <- sum((r - m_i[cell])^2) / (N - k); n0 <- (N - sum(n_i^2) / N) / (k - 1)
  sb <- max(0, (msb - msw) / n0); sb / (sb + msw) }
base <- lapply(c(-0.015, 0.015), function(b) {
  d <- tab; d$y0 <- d$DORM_z - b * d$P0
  f0 <- phylolm(y0 ~ domain + lgs + checkm2_completeness + is_mag, data = d, phy = ph, model = "lambda")
  list(b = b, fit = unname(f0$fitted.values), res = unname(f0$residuals))
})
procP <- function(d, grp, extra) {
  fl <- phylolm(form, data = d, phy = ph, model = "lambda"); lam <- fl$optpar; b <- coef(fl)["P0"]
  pl <- transf.branch.lengths(ph, model = "lambda", parameters = list(lambda = lam))$tree
  jse <- function(gr) { est <- sapply(seq_len(G), function(g) { keep <- gr[d$cell] != g; p2 <- keep.tip(pl, d$accession[keep])
    coef(phylolm(form, data = d[p2$tip.label, ], phy = p2, model = "BM"))["P0"] }); sqrt((G - 1) / G * sum((est - mean(est))^2)) }
  q <- qt(0.95, G - 1); dec <- function(s) as.numeric(b - q * s > -SESOI & b + q * s < SESOI)
  ses <- c(jse(grp), sapply(seq_len(extra), function(i) { g2 <- sample(rep_len(seq_len(G), length(cells))); names(g2) <- cells; jse(g2) }))
  c(lambda_hat = lam, estimate = unname(b), se = ses[1], equiv = dec(ses[1]), se_avg = sqrt(mean(ses^2)), equiv_avg = dec(sqrt(mean(ses^2))),
    between_cell_share = between_share(unname(fl$residuals), d$cell))
}
one <- function(seed, k) {
  set.seed(seed); B <- base[[k]]; w <- sample(c(-1, 1), length(cells), replace = TRUE)[cidx]
  d <- tab; d$y <- B$b * tab$P0 + B$fit + w * B$res
  grp <- sample(rep_len(seq_len(G), length(cells))); names(grp) <- cells
  c(dgm = 3, true_beta = B$b, seed = seed, procP(d, grp, K - 1))
}
simf <- paste0(a[3], "_sims.tsv")
done <- if (file.exists(simf)) fread(simf) else data.table(true_beta = numeric(0))
for (k in 1:2) {
  b <- base[[k]]$b; have <- sum(abs(done$true_beta - b) < 1e-12)
  while (have < NSIM) {
    n <- min(40, NSIM - have); seeds <- 5e6 + 1e5 * k + have + seq_len(n)
    r <- as.data.frame(do.call(rbind, mclapply(seeds, one, k = k, mc.cores = cores)))
    fwrite(r, simf, sep = "\t", append = file.exists(simf)); have <- have + n
    cat("DGM C beta", b, "done", have, format(Sys.time()), "\n")
  }
}
s <- fread(simf)
summ <- s[, .(nsim = .N, n_equiv = sum(equiv), P_declare_equivalence = mean(equiv), n_equiv_avg = sum(equiv_avg),
              P_equiv_se_avg = mean(equiv_avg), mean_lambda_hat = mean(lambda_hat), sd_estimate = sd(estimate), mean_estimate = mean(estimate),
              mean_se = mean(se), median_se_avg = median(se_avg), mean_between_cell_share = mean(between_cell_share)), by = true_beta]
summ[, `:=`(mc_lo95 = qbeta(0.025, n_equiv, nsim - n_equiv + 1), mc_hi95 = qbeta(0.975, n_equiv + 1, nsim - n_equiv),
            avg_mc_lo95 = qbeta(0.025, n_equiv_avg, nsim - n_equiv_avg + 1), avg_mc_hi95 = qbeta(0.975, n_equiv_avg + 1, nsim - n_equiv_avg))]
fwrite(summ, paste0(a[3], "_summary.tsv"), sep = "\t"); print(summ)
