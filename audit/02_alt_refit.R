#!/usr/bin/env Rscript
# audit task 4: re-fit H1 (direct) and the MAG x P0 interaction (plus MAG-only / isolate-only
# slopes) with methods different from the original phylolm(lambda):
#  A  phylolm OU (OUfixedRoot) on the GTDB branch lengths
#  B  lmer with nested taxonomic random intercepts (phylum/class/order/family/genus) + crossed
#     site (0.5-degree cell) random intercept: phylogeny AND site non-independence together
#  C  MCMCglmm with phylogenetic (animal, GTDB branch lengths, NOT Grafen) + site random effects
# Data: data_small/final/analysis_table.tsv.gz (pilot == 0), tree data_small/analysis_tree.nwk.
# usage (repo root): Rscript audit/02_alt_refit.R OUTDIR [A,B,C] [nitt]
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table); library(lme4); library(lmerTest)})
a <- commandArgs(TRUE); outd <- a[1]; which <- if (length(a) > 1) strsplit(a[2], ",")[[1]] else c("A", "B", "C")
nitt <- if (length(a) > 2) as.integer(a[3]) else 33000
tab <- as.data.frame(fread("data_small/final/analysis_table.tsv.gz")); tr <- read.tree("data_small/analysis_tree.nwk")
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
tab$P0 <- z(tab$log_T_sd_total, tab$domain)
tab$lgs <- log10(tab$genome_size); tab$is_mag <- as.integer(as.logical(tab$is_mag))
tab$site <- factor(paste(tab$domain, round(tab$lat * 2) / 2, round(tab$lon * 2) / 2))
rk <- function(r) sub(paste0(".*;", r, "__([^;]*).*"), "\\1", tab$gtdb_taxonomy)
for (r in c("p", "c", "o", "f", "g")) tab[[paste0("tx_", r)]] <- factor(rk(r))
rownames(tab) <- tab$accession
cat("n", nrow(tab), "sites", nlevels(tab$site), "genera", nlevels(tab$tx_g), "\n")
specs <- list(
  H1_direct   = list(d = function(x) x, f = "DORM_z ~ P0 + domain + lgs + checkm2_completeness + is_mag", term = "P0"),
  interaction = list(d = function(x) x, f = "DORM_z ~ P0 * is_mag + domain + lgs + checkm2_completeness", term = "P0:is_mag"),
  MAGs        = list(d = function(x) x[x$is_mag == 1, ], f = "DORM_z ~ P0 + domain + lgs + checkm2_completeness", term = "P0"),
  isolates    = list(d = function(x) x[x$is_mag == 0, ], f = "DORM_z ~ P0 + domain + lgs + checkm2_completeness", term = "P0"))
res <- list(); dir.create(outd, showWarnings = FALSE, recursive = TRUE)
save <- function() fwrite(do.call(rbind, res), file.path(outd, "alt_refit.tsv"), sep = "\t")
rowf <- function(method, nm, d, est, se, p, lo, hi, extra = "") data.frame(method = method, model = nm, n = nrow(d),
  estimate = est, se = se, p = p, l95 = lo, u95 = hi, extra = extra)
for (nm in names(specs)) {
  s <- specs[[nm]]; d <- s$d(tab); ph <- keep.tip(tr, d$accession); d <- d[ph$tip.label, ]
  if ("A" %in% which) {
    t0 <- Sys.time()
    fit <- tryCatch(phylolm(as.formula(s$f), data = d, phy = ph, model = "OUfixedRoot"), error = function(e) e)
    if (!inherits(fit, "error")) {
      co <- summary(fit)$coefficients
      res[[length(res) + 1]] <- rowf("A_phylolm_OU", nm, d, co[s$term, 1], co[s$term, 2], co[s$term, 4],
        co[s$term, 1] - 1.96 * co[s$term, 2], co[s$term, 1] + 1.96 * co[s$term, 2],
        sprintf("alpha=%.4g; AIC_OU=%.1f; AIC_lambda=%.1f", fit$optpar, AIC(fit),
                AIC(phylolm(as.formula(s$f), data = d, phy = ph, model = "lambda"))))
    } else res[[length(res) + 1]] <- rowf("A_phylolm_OU", nm, d, NA, NA, NA, NA, NA, conditionMessage(fit))
    print(tail(res, 1)[[1]]); save(); cat("A", nm, format(Sys.time() - t0), "\n")
  }
  if ("B" %in% which) {
    fB <- as.formula(paste(s$f, "+ (1 | tx_p/tx_c/tx_o/tx_f/tx_g) + (1 | site)"))
    m <- lmer(fB, data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
    co <- summary(m)$coefficients; vc <- as.data.frame(VarCorr(m))
    res[[length(res) + 1]] <- rowf("B_lmer_taxonomy_site", nm, d, co[s$term, 1], co[s$term, 2], co[s$term, 5],
      co[s$term, 1] - 1.96 * co[s$term, 2], co[s$term, 1] + 1.96 * co[s$term, 2],
      paste(sprintf("%s=%.4f", vc$grp, vc$vcov), collapse = "; "))
    print(tail(res, 1)[[1]]); save()
    # same model without the site term: isolates the effect of modelling site
    m0 <- lmer(as.formula(paste(s$f, "+ (1 | tx_p/tx_c/tx_o/tx_f/tx_g)")), data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
    co <- summary(m0)$coefficients
    res[[length(res) + 1]] <- rowf("B0_lmer_taxonomy_only", nm, d, co[s$term, 1], co[s$term, 2], co[s$term, 5],
      co[s$term, 1] - 1.96 * co[s$term, 2], co[s$term, 1] + 1.96 * co[s$term, 2], "")
    print(tail(res, 1)[[1]]); save()
  }
}
if ("C" %in% which) {
  suppressPackageStartupMessages(library(MCMCglmm))
  set.seed(20261117)
  pr <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                                   G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
  for (nm in c("interaction", "H1_direct")) {
    s <- specs[[nm]]; d <- s$d(tab); ph <- keep.tip(tr, d$accession); ph$node.label <- NULL
    ph$edge.length[ph$edge.length <= 0] <- 1e-8
    Ainv <- inverseA(ph, nodes = "ALL", scale = FALSE)$Ainv   # GTDB branch lengths kept
    d$animal <- d$accession; d$site <- droplevels(d$site); t0 <- Sys.time()
    m <- MCMCglmm(as.formula(s$f), random = ~animal + site, ginverse = list(animal = Ainv), data = d, prior = pr,
                  nitt = nitt, burnin = 3000, thin = 30, verbose = FALSE)
    sm <- summary(m)$solutions; vc <- colMeans(m$VCV)
    res[[length(res) + 1]] <- rowf("C_MCMCglmm_phylo_site_GTDBbl", nm, d, sm[s$term, 1], sd(m$Sol[, s$term]), sm[s$term, 5],
      sm[s$term, 2], sm[s$term, 3], sprintf("ESS=%.0f; var_animal=%.4f; var_site=%.4f; var_resid=%.4f; nitt=%d; %s",
      effectiveSize(m$Sol[, s$term]), vc["animal"], vc["site"], vc["units"], nitt, format(Sys.time() - t0)))
    print(tail(res, 1)[[1]]); save()
  }
}
