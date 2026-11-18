#!/usr/bin/env Rscript
# Exploratory (NOT pre-registered; see results/deviations.md): robustness of H1 to pseudoreplication of
# isolation sites. Many species (especially MAGs) share one sampling site and therefore one predictor value,
# so phylolm SEs treat them as independent environmental replicates.
#   (1) delete-a-group jackknife over sites (G groups of whole sites) -> cluster-robust SE
#   (2) site-capped resampling: at most K species per site, R replicates -> distribution of beta
#   (3) MCMCglmm Gaussian with phylogenetic (animal) + site random effects
# usage: Rscript 16_cluster_robust.R TABLE.tsv.gz TREE.nwk OUTPREFIX [G=100] [R=50] [K=3] [mcmc]
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table)})
a <- commandArgs(TRUE)
G <- if (length(a) > 3) as.integer(a[4]) else 100
R <- if (length(a) > 4) as.integer(a[5]) else 50
K <- if (length(a) > 5) as.integer(a[6]) else 3
do_mcmc <- length(a) > 6 && a[7] == "mcmc"
tab <- as.data.frame(fread(a[1])); tr <- read.tree(a[2]); outp <- a[3]
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
tab$P0 <- z(tab$log_T_sd_total, tab$domain)
tab$lgs <- log10(tab$genome_size); tab$is_mag <- as.integer(as.logical(tab$is_mag))
# a site = 0.5-degree cell (the land climate grid; marine 0.25-degree cells nest inside it)
tab$site <- paste(tab$domain, round(tab$lat * 2) / 2, round(tab$lon * 2) / 2)
tab <- tab[complete.cases(tab[, c("DORM_z", "P0", "lgs", "checkm2_completeness")]), ]
rownames(tab) <- tab$accession
cat("species", nrow(tab), "sites", length(unique(tab$site)), "\n")

fitb <- function(d, covs) {
  ph <- keep.tip(tr, d$accession); d <- d[ph$tip.label, ]
  if (length(unique(d$domain)) < 2) covs <- setdiff(covs, "domain")
  if (length(unique(d$is_mag)) < 2) covs <- setdiff(covs, "is_mag")
  f <- as.formula(paste("DORM_z ~", paste(c("P0", covs), collapse = " + ")))
  co <- summary(phylolm(f, data = d, phy = ph, model = "lambda"))$coefficients
  co["P0", 1:2]
}
sets <- list(
  H1_direct = list(filter = function(d) d, covs = c("domain", "lgs", "checkm2_completeness", "is_mag")),
  H1_total  = list(filter = function(d) d, covs = c("domain", "checkm2_completeness", "is_mag")),
  MAGs      = list(filter = function(d) d[d$is_mag == 1, ], covs = c("domain", "lgs", "checkm2_completeness")),
  isolates  = list(filter = function(d) d[d$is_mag == 0, ], covs = c("domain", "lgs", "checkm2_completeness")))
set.seed(20261117)
out <- list(); boot <- list()
for (nm in names(sets)) {
  d0 <- sets[[nm]]$filter(tab); covs <- sets[[nm]]$covs
  full <- fitb(d0, covs)
  sites <- unique(d0$site); grp <- sample(rep_len(seq_len(G), length(sites))); names(grp) <- sites
  jk <- sapply(seq_len(G), function(g) fitb(d0[grp[d0$site] != g, ], covs)[1])
  se_jk <- sqrt((G - 1) / G * sum((jk - mean(jk))^2))
  capb <- sapply(seq_len(R), function(r) {
    idx <- unlist(lapply(split(seq_len(nrow(d0)), d0$site), function(i) if (length(i) > K) sample(i, K) else i))
    fitb(d0[idx, ], covs)[1]
  })
  out[[nm]] <- data.frame(set = nm, n = nrow(d0), n_sites = length(sites), beta = full[1], se_phylolm = full[2],
                          se_site_jackknife = se_jk, p_site_jackknife = 2 * pnorm(-abs(full[1] / se_jk)),
                          ci_lo_jk = full[1] - 1.96 * se_jk, ci_hi_jk = full[1] + 1.96 * se_jk,
                          cap_K = K, cap_n_mean = NA, cap_beta_median = median(capb),
                          cap_beta_q025 = quantile(capb, 0.025), cap_beta_q975 = quantile(capb, 0.975),
                          cap_frac_positive = mean(capb > 0))
  boot[[nm]] <- data.frame(set = nm, rep = seq_len(R), beta = capb)
  print(out[[nm]])
  fwrite(do.call(rbind, out), paste0(outp, "_cluster_robust.tsv"), sep = "\t")
}
fwrite(do.call(rbind, boot), paste0(outp, "_site_capped_betas.tsv"), sep = "\t")

if (do_mcmc) {
  suppressPackageStartupMessages(library(MCMCglmm))
  ph <- keep.tip(tr, tab$accession); ph$node.label <- NULL
  ph$edge.length[ph$edge.length <= 0] <- 1e-8   # GTDB branch lengths kept (audit M1); no Grafen rescaling
  Ainv <- inverseA(ph, nodes = "ALL", scale = FALSE)$Ainv
  tab$animal <- tab$accession; tab$site <- factor(tab$site)
  pr <- list(R = list(V = 1, nu = 0.002),
             G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                      G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
  m <- MCMCglmm(DORM_z ~ P0 + domain + lgs + checkm2_completeness + is_mag, random = ~animal + site,
                ginverse = list(animal = Ainv), data = tab, prior = pr,
                nitt = 33000, burnin = 3000, thin = 30, verbose = FALSE)
  s <- summary(m)$solutions["P0", ]
  vc <- colMeans(m$VCV)
  res <- data.frame(model = "MCMCglmm_DORM_z_animal_site", n = nrow(tab), post_mean = s[1], l95 = s[2], u95 = s[3],
                    pMCMC = s[5], var_animal = vc["animal"], var_site = vc["site"], var_resid = vc["units"],
                    ess_P0 = effectiveSize(m$Sol[, "P0"]))
  print(res)
  fwrite(res, paste0(outp, "_mcmcglmm_site.tsv"), sep = "\t")
}
