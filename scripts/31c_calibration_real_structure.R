#!/usr/bin/env Rscript
# Added 18 Nov 2026: calibration of procedure P under "DGM D", a data-generating model built from the REAL
# residuals that keeps both their phylogenetic structure and their clustering within 0.5-degree cells.
# Background: a naive wild cluster bootstrap (sign flip of whole residuals per cell; scripts/31b, results/
# calibration_wild_cluster_naive_*) reproduced the within-cell clustering but destroyed the phylogenetic correlation
# between cells (lambda-hat about 0.04 instead of 0.94), so it was not informative and was stopped after 80 runs.
# Decomposition of the real lambda-model residuals r (scratch check, reproduced below): r = u + eps with
# eps = sigma2 (1 - lambda) D V^-1 r the non-phylogenetic part (BLUP) and u the phylogenetic part. The excess
# within-cell clustering of the real data is in u (between-cell share 0.45), not in eps (0.06).
# DGM D: y* = b * P0 + (fitted values without the P0 term) + w_order * u + w_cell * eps, with independent Rademacher
# signs per GTDB order (639 orders) and per 0.5-degree cell. Pilot draws: lambda-hat 0.96-0.98, between-cell share
# 0.35-0.38 (real data 0.94 and 0.41).
# Why: the parametric DGMs of scripts/26 and 31 reproduce lambda-hat and the SD of the slope, but not the clustering of
# residuals within 0.5-degree cells (between-cell share of residual variance about 0.13-0.17 in simulations, 0.50 in the
# real data for the same model form). Here:
#   procedure P is applied with lambda re-estimated, with the original single-grouping SE and with the SE averaged over
#      K = 5 random groupings (as in scripts/31_calibration_lambda_free.R).
# 500 simulations at each boundary; restart-safe.
# usage: Rscript 31c_calibration_real_structure.R TABLE.tsv.gz TREE.nwk OUT_PREFIX [cores=4] -> OUT_PREFIX_sims.tsv, _summary.tsv
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
for (lv in "o") tab$order <- sub(".*;o__([^;]*);.*", "\\1", tab$gtdb_taxonomy)
oidx <- match(tab$order, unique(tab$order)); nord <- length(unique(tab$order))
full <- phylolm(DORM_z ~ P0 + domain + lgs + checkm2_completeness + is_mag, data = tab, phy = ph, model = "lambda")
lam <- full$optpar; s2 <- full$sigma2; r <- full$residuals
C <- vcv(ph)[tab$accession, tab$accession]; dg <- diag(C)
eps <- s2 * (1 - lam) * dg * solve(s2 * (lam * C + diag((1 - lam) * dg)), r); u <- r - eps; rm(C)
fit0 <- full$fitted.values - coef(full)["P0"] * tab$P0
cat("real residual between-cell share: total", between_share(r, tab$cell), "phylogenetic part", between_share(u, tab$cell),
    "non-phylogenetic part", between_share(eps, tab$cell), "\n")
fwrite(data.frame(part = c("total", "phylogenetic (u)", "non-phylogenetic (eps)"), var = c(var(r), var(u), var(eps)),
                  between_cell_share = c(between_share(r, tab$cell), between_share(u, tab$cell), between_share(eps, tab$cell))),
       paste0(a[3], "_real_decomposition.tsv"), sep = "\t")
base <- lapply(c(-0.015, 0.015), function(b) list(b = b))
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
  set.seed(seed); B <- base[[k]]
  wo <- sample(c(-1, 1), nord, replace = TRUE)[oidx]; wc <- sample(c(-1, 1), length(cells), replace = TRUE)[cidx]
  d <- tab; d$y <- B$b * tab$P0 + fit0 + wo * u + wc * eps
  grp <- sample(rep_len(seq_len(G), length(cells))); names(grp) <- cells
  c(dgm = 4, true_beta = B$b, seed = seed, procP(d, grp, K - 1))
}
simf <- paste0(a[3], "_sims.tsv")
done <- if (file.exists(simf)) fread(simf) else data.table(true_beta = numeric(0))
for (k in 1:2) {
  b <- base[[k]]$b; have <- sum(abs(done$true_beta - b) < 1e-12)
  while (have < NSIM) {
    n <- min(40, NSIM - have); seeds <- 7e6 + 1e5 * k + have + seq_len(n)
    r <- as.data.frame(do.call(rbind, mclapply(seeds, one, k = k, mc.cores = cores)))
    fwrite(r, simf, sep = "\t", append = file.exists(simf)); have <- have + n
    cat("DGM D beta", b, "done", have, format(Sys.time()), "\n")
  }
}
s <- fread(simf)
summ <- s[, .(nsim = .N, n_equiv = sum(equiv), P_declare_equivalence = mean(equiv), n_equiv_avg = sum(equiv_avg),
              P_equiv_se_avg = mean(equiv_avg), mean_lambda_hat = mean(lambda_hat), sd_estimate = sd(estimate), mean_estimate = mean(estimate),
              mean_se = mean(se), median_se_avg = median(se_avg), mean_between_cell_share = mean(between_cell_share)), by = true_beta]
summ[, `:=`(mc_lo95 = qbeta(0.025, n_equiv, nsim - n_equiv + 1), mc_hi95 = qbeta(0.975, n_equiv + 1, nsim - n_equiv),
            avg_mc_lo95 = qbeta(0.025, n_equiv_avg, nsim - n_equiv_avg + 1), avg_mc_hi95 = qbeta(0.975, n_equiv_avg + 1, nsim - n_equiv_avg))]
fwrite(summ, paste0(a[3], "_summary.tsv"), sep = "\t"); print(summ)
