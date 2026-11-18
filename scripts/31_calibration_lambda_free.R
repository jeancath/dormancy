#!/usr/bin/env Rscript
# Added 18 Nov 2026: calibration of procedure P with lambda RE-ESTIMATED in every simulated data set, exactly as
# P is applied to the real data (lambda-hat from the data, then fixed for the jackknife), and with two data-generating
# models (DGMs):
#   A  "original": Brownian motion on the lambda-hat-transformed tree with sigma2 from the real phylolm fit, plus an
#      independent effect per 0.5-degree site with the SD from MCMCglmm (as in scripts/26_tost_calibration.R).
#   B  "single fitted model": the variance components of the ONE model fitted to the real data with phylogeny, site and
#      residual together (MCMCglmm, results/siteMCMC_mcmcglmm_site.tsv): phylogenetic effect ~ BM on the GTDB tree
#      (branch lengths unscaled, as in the fit) with variance var_animal, site effect with var_site, residual var_resid.
#      Nothing is counted twice.
# Recorded per simulation, to compare the variance structure with the real data: lambda-hat, the slope, the jackknife SE,
# and the share of the GLS residual variance that lies between 0.5-degree cells (one-way ANOVA estimator).
# 500 simulations per DGM at each boundary (-0.015, +0.015); restart-safe (appends chunks of 40).
# For DGM B each simulation also evaluates a DE-RANDOMISED variant of P: the jackknife variance averaged over K = 5
# independent random groupings, SE_avg = sqrt(mean(SE_k^2)), still referred to t on 19 df. (It was found that on the
# real data the SE from a single grouping varies widely between groupings; the averaged SE removes most of that.)
# usage: Rscript 31_calibration_lambda_free.R TABLE.tsv.gz TREE.nwk OUT_PREFIX [cores=4]
#   -> OUT_PREFIX_sims.tsv, OUT_PREFIX_real.tsv (same descriptors on the real data), OUT_PREFIX_summary.tsv
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table); library(parallel)})
a <- commandArgs(TRUE); cores <- if (length(a) > 3) as.integer(a[4]) else 4
G <- 20; SESOI <- 0.015; NSIM <- 500; K <- 5
tab <- as.data.frame(fread(a[1])); tr <- read.tree(a[2])
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
tab$P0 <- z(tab$log_T_sd_total, tab$domain); tab$lgs <- log10(tab$genome_size)
tab$is_mag <- as.integer(as.logical(tab$is_mag)); rownames(tab) <- tab$accession
ph <- keep.tip(tr, tab$accession); tab <- tab[ph$tip.label, ]
ph$edge.length[ph$edge.length <= 0] <- 1e-8
tab$cell <- paste(tab$domain, floor(tab$lat / 0.5), floor(tab$lon / 0.5))
tab$site <- paste(tab$domain, round(tab$lat * 2) / 2, round(tab$lon * 2) / 2)
cells <- unique(tab$cell); sites <- unique(tab$site); sidx <- match(tab$site, sites)
mc <- fread("results/siteMCMC_mcmcglmm_site.tsv")
full <- phylolm(DORM_z ~ P0 + domain + lgs + checkm2_completeness + is_mag, data = tab, phy = ph, model = "lambda")
lam0 <- full$optpar; s2 <- full$sigma2
phl0 <- transf.branch.lengths(ph, model = "lambda", parameters = list(lambda = lam0))$tree
between_share <- function(r, cell) {  # one-way ANOVA estimator of the between-cell variance share of residuals
  n_i <- tapply(r, cell, length); m_i <- tapply(r, cell, mean); k <- length(n_i); N <- length(r)
  msb <- sum(n_i * (m_i - mean(r))^2) / (k - 1); msw <- sum((r - m_i[cell])^2) / (N - k)
  n0 <- (N - sum(n_i^2) / N) / (k - 1); sb <- max(0, (msb - msw) / n0); sb / (sb + msw)
}
procP <- function(d, form, grp, extra = 0) {   # lambda estimated on d, then fixed for the jackknife (procedure P)
  fl <- phylolm(form, data = d, phy = ph, model = "lambda"); lam <- fl$optpar; b <- coef(fl)["P0"]
  pl <- transf.branch.lengths(ph, model = "lambda", parameters = list(lambda = lam))$tree
  jse <- function(gr) { est <- sapply(seq_len(G), function(g) { keep <- gr[d$cell] != g; p2 <- keep.tip(pl, d$accession[keep])
    coef(phylolm(form, data = d[p2$tip.label, ], phy = p2, model = "BM"))["P0"] }); sqrt((G - 1) / G * sum((est - mean(est))^2)) }
  se <- jse(grp); q <- qt(0.95, G - 1)
  dec <- function(s) as.numeric(b - q * s > -SESOI & b + q * s < SESOI)
  out <- c(lambda_hat = lam, estimate = unname(b), se = se,
    p_tost = max(pt((b + SESOI) / se, G - 1, lower.tail = FALSE), pt((SESOI - b) / se, G - 1, lower.tail = FALSE)),
    equiv = dec(se), between_cell_share = between_share(unname(fl$residuals), d$cell), se_avg = NA, equiv_avg = NA)
  if (extra > 0) {
    ses <- c(se, sapply(seq_len(extra), function(i) { g2 <- sample(rep_len(seq_len(G), length(cells))); names(g2) <- cells; jse(g2) }))
    out["se_avg"] <- sqrt(mean(ses^2)); out["equiv_avg"] <- dec(out["se_avg"])
  }
  out
}
realf <- paste0(a[3], "_real.tsv")
if (!file.exists(realf)) {
  set.seed(20261118); grp <- sample(rep_len(seq_len(G), length(cells))); names(grp) <- cells
  r <- procP(tab, DORM_z ~ P0 + domain + lgs + checkm2_completeness + is_mag, grp)
  fwrite(data.frame(t(r), sigma2_lambda_model = s2), realf, sep = "\t"); print(r)
}
one <- function(seed, b, dgm) {
  set.seed(seed)
  e <- if (dgm == "A") rTrait(n = 1, phl0, model = "BM", parameters = list(sigma2 = s2))[tab$accession] + rnorm(length(sites), 0, sqrt(mc$var_site))[sidx]
       else rTrait(n = 1, ph, model = "BM", parameters = list(sigma2 = mc$var_animal))[tab$accession] +
            rnorm(length(sites), 0, sqrt(mc$var_site))[sidx] + rnorm(nrow(tab), 0, sqrt(mc$var_resid))
  d <- data.frame(accession = tab$accession, y = b * tab$P0 + e, P0 = tab$P0, domain = tab$domain, cell = tab$cell, row.names = tab$accession)
  grp <- sample(rep_len(seq_len(G), length(cells))); names(grp) <- cells
  c(dgm = ifelse(dgm == "A", 1, 2), true_beta = b, seed = seed, procP(d, y ~ P0 + domain, grp, extra = ifelse(dgm == "B", K - 1, 0)))
}
simf <- paste0(a[3], "_sims.tsv")
done <- if (file.exists(simf)) fread(simf) else data.table(dgm = numeric(0), true_beta = numeric(0))
plan <- expand.grid(b = c(-0.015, 0.015), dgm = c("B", "A"), stringsAsFactors = FALSE)
for (k in seq_len(nrow(plan))) {
  b <- plan$b[k]; dg <- plan$dgm[k]; dn <- ifelse(dg == "A", 1, 2)
  have <- sum(done$dgm == dn & abs(done$true_beta - b) < 1e-12)
  while (have < NSIM) {
    n <- min(40, NSIM - have); seeds <- 3e6 + 1e5 * k + have + seq_len(n)
    r <- as.data.frame(do.call(rbind, mclapply(seeds, one, b = b, dgm = dg, mc.cores = cores)))
    fwrite(r, simf, sep = "\t", append = file.exists(simf)); have <- have + n
    cat("DGM", dg, "beta", b, "done", have, format(Sys.time()), "\n")
  }
}
s <- fread(simf)
summ <- s[, .(nsim = .N, P_declare_equivalence = mean(equiv), n_equiv = sum(equiv), mean_lambda_hat = mean(lambda_hat),
              sd_estimate = sd(estimate), mean_se = mean(se), median_se = median(se), q05_se = quantile(se, 0.05), q95_se = quantile(se, 0.95),
              mean_between_cell_share = mean(between_cell_share),
              P_equiv_se_avg = mean(equiv_avg), n_equiv_avg = sum(equiv_avg), median_se_avg = median(se_avg)), by = .(dgm, true_beta)]
summ[, `:=`(mc_lo95 = qbeta(0.025, n_equiv, nsim - n_equiv + 1), mc_hi95 = qbeta(0.975, n_equiv + 1, nsim - n_equiv))]
summ[!is.na(n_equiv_avg), `:=`(avg_mc_lo95 = qbeta(0.025, n_equiv_avg, nsim - n_equiv_avg + 1), avg_mc_hi95 = qbeta(0.975, n_equiv_avg + 1, nsim - n_equiv_avg))]
summ[, dgm := ifelse(dgm == 1, "A original", "B single fitted model")]
fwrite(summ, paste0(a[3], "_summary.tsv"), sep = "\t"); print(summ)
