#!/usr/bin/env Rscript
# Calibration of the PRIMARY cluster-robust equivalence procedure (added 18 Nov 2026).
#
# Primary procedure P (applied identically to the real data and to every simulated data set):
#   1. lambda-hat from phylolm(model = "lambda") on the full real data (fixed thereafter);
#   2. slope from phylolm(model = "BM") on the lambda-hat-transformed tree (identical point estimate);
#   3. SE from a delete-a-group jackknife over G = 20 random groups of 0.5-degree grid cells (floor-aligned),
#      refitting with lambda fixed;
#   4. TOST against +-0.015 with a t distribution on G - 1 = 19 df (equivalent iff the 90 % t-CI lies inside).
# Data-generating model for the simulations: real tree and real P0 values; phylogenetic error (Brownian motion on
# the lambda-hat tree, sigma2 from the real fit); site effect per 0.5-degree cell with the SD estimated by the
# MCMCglmm phylogeny + site model; true slopes {-0.015, +0.015} (Type I error at the boundaries) and {0, 0.005, 0.010}.
# Results are appended per chunk (restart-safe).
# usage: Rscript 26_tost_calibration.R TABLE.tsv.gz TREE.nwk OUT_PREFIX [cores=3]
#   -> OUT_PREFIX_real.tsv (procedure P on the real data), OUT_PREFIX_sims.tsv (one row per simulation),
#      OUT_PREFIX_summary.tsv (rates with Monte Carlo 95 % CIs; the committed summary was recomputed in Python with
#      Clopper-Pearson CIs, see results/REPORT.md section 3.10)
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table); library(parallel)})
a <- commandArgs(TRUE); cores <- if (length(a) > 3) as.integer(a[4]) else 3
G <- 20; SESOI <- 0.015
plan <- data.frame(beta = c(-0.015, 0.015, 0, 0.005, 0.010), nsim = c(5000, 5000, 2000, 1000, 1000))
tab <- as.data.frame(fread(a[1])); tr <- read.tree(a[2])
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
tab$P0 <- z(tab$log_T_sd_total, tab$domain); tab$lgs <- log10(tab$genome_size)
tab$is_mag <- as.integer(as.logical(tab$is_mag)); rownames(tab) <- tab$accession
ph <- keep.tip(tr, tab$accession); tab <- tab[ph$tip.label, ]
tab$cell <- paste(tab$domain, floor(tab$lat / 0.5), floor(tab$lon / 0.5))          # jackknife blocks
tab$site <- paste(tab$domain, round(tab$lat * 2) / 2, round(tab$lon * 2) / 2)     # simulated site effects
full <- phylolm(DORM_z ~ P0 + domain + lgs + checkm2_completeness + is_mag, data = tab, phy = ph, model = "lambda")
lam <- full$optpar; s2 <- full$sigma2
phl <- transf.branch.lengths(ph, model = "lambda", parameters = list(lambda = lam))$tree
site_sd <- sqrt(fread("results/siteMCMC_mcmcglmm_site.tsv")$var_site[1])
cells <- unique(tab$cell)

procP <- function(d, form, grp) {
  b <- coef(phylolm(form, data = d, phy = phl, model = "BM"))["P0"]
  est <- sapply(seq_len(G), function(g) {
    keep <- grp[d$cell] != g; p2 <- keep.tip(phl, d$accession[keep])
    coef(phylolm(form, data = d[p2$tip.label, ], phy = p2, model = "BM"))["P0"]
  })
  se <- sqrt((G - 1) / G * sum((est - mean(est))^2)); q <- qt(0.95, G - 1)
  c(beta = unname(b), se = se, lo90 = unname(b) - q * se, hi90 = unname(b) + q * se,
    p_tost = max(pt((b + SESOI) / se, G - 1, lower.tail = FALSE), pt((SESOI - b) / se, G - 1, lower.tail = FALSE)),
    equiv = as.numeric(b - q * se > -SESOI & b + q * se < SESOI))
}

# procedure P on the real data (primary, direct and total effect; MAGs; isolates)
if (!file.exists(paste0(a[3], "_real.tsv"))) {
  set.seed(20261118); grp <- sample(rep_len(seq_len(G), length(cells))); names(grp) <- cells
  forms <- list(H1_direct = DORM_z ~ P0 + domain + lgs + checkm2_completeness + is_mag,
                H1_total = DORM_z ~ P0 + domain + checkm2_completeness + is_mag)
  real <- do.call(rbind, lapply(names(forms), function(n) data.frame(set = n, n = nrow(tab), t(procP(tab, forms[[n]], grp)))))
  for (gt in c(1, 0)) {
    d <- tab[tab$is_mag == gt, ]; keep_cells <- unique(d$cell)
    g2 <- grp[keep_cells]
    # refit lambda-free BM on the subset tree
    phl_bak <- phl; phl <<- keep.tip(phl_bak, d$accession)
    real <- rbind(real, data.frame(set = ifelse(gt == 1, "MAGs", "isolates"), n = nrow(d),
                                   t(procP(d[phl$tip.label, ], DORM_z ~ P0 + domain + lgs + checkm2_completeness, g2))))
    phl <<- phl_bak
  }
  real$lambda <- lam; real$G <- G
  fwrite(real, paste0(a[3], "_real.tsv"), sep = "\t"); print(real)
}

# simulations (restart-safe)
simf <- paste0(a[3], "_sims.tsv")
done <- if (file.exists(simf)) setnames(fread(simf), 1, "true_beta") else data.frame(true_beta = numeric(0))
sites <- unique(tab$site); sidx <- match(tab$site, sites)
form_sim <- y ~ P0 + domain
one <- function(seed, b) {
  set.seed(seed)
  e <- rTrait(n = 1, phl, model = "BM", parameters = list(sigma2 = s2))[tab$accession]
  d <- data.frame(accession = tab$accession, y = b * tab$P0 + e + rnorm(length(sites), 0, site_sd)[sidx],
                  P0 = tab$P0, domain = tab$domain, cell = tab$cell, row.names = tab$accession)
  grp <- sample(rep_len(seq_len(G), length(cells))); names(grp) <- cells
  c(true_beta = b, seed = seed, procP(d, form_sim, grp))
}
for (k in seq_len(nrow(plan))) {
  b <- plan$beta[k]; have <- sum(abs(done$true_beta - b) < 1e-12, na.rm = TRUE)
  while (have < plan$nsim[k]) {
    n <- min(200, plan$nsim[k] - have)
    seeds <- 1e6 * (k + 1) + have + seq_len(n)
    r <- do.call(rbind, mclapply(seeds, one, b = b, mc.cores = cores))
    r <- as.data.frame(r)   # columns: true_beta, seed, beta (estimate), se, lo90, hi90, p_tost, equiv
    fwrite(r, simf, sep = "\t", append = file.exists(simf))
    have <- have + n; cat("beta", b, "done", have, "\n")
  }
}
s <- fread(simf); setnames(s, 1:3, c("true_beta", "seed", "estimate"))   # older runs wrote two "beta" headers
summ <- s[, .(nsim = .N, P_declare_equivalence = mean(equiv), mean_estimate = mean(estimate), sd_estimate = sd(estimate),
              mean_jackknife_se = mean(se)), by = true_beta]
summ[, mc_se := sqrt(P_declare_equivalence * (1 - P_declare_equivalence) / nsim)]
summ[, `:=`(mc_lo95 = pmax(0, P_declare_equivalence - 1.96 * mc_se), mc_hi95 = pmin(1, P_declare_equivalence + 1.96 * mc_se))]
fwrite(summ, paste0(a[3], "_summary.tsv"), sep = "\t"); print(summ)
