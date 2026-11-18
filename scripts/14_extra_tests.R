#!/usr/bin/env Rscript
# Pre-registered extras on the confirmatory set:
#  (1) H6 joint model DORM_z ~ seasonal + residual (+covariates) and Wald test of beta_resid - beta_seasonal
#  (2) Sensitivity (x): MCMCglmm Poisson per module with a phylogenetic random effect (if time allows)
# usage: Rscript 14_extra_tests.R TABLE.tsv.gz TREE.nwk OUTPREFIX [mcmc]
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table)})
a <- commandArgs(TRUE)
tab <- as.data.frame(fread(a[1])); tr <- read.tree(a[2]); outp <- a[3]
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
tab$seas <- z(tab$log_T_sd_seasonal, tab$domain); tab$resid <- z(tab$log_T_sd_resid, tab$domain)
tab$P0 <- z(tab$log_T_sd_total, tab$domain)
tab$lgs <- log10(tab$genome_size); tab$is_mag <- as.integer(as.logical(tab$is_mag))
rownames(tab) <- tab$accession
ph <- keep.tip(tr, tab$accession); tab <- tab[ph$tip.label, ]
f <- phylolm(DORM_z ~ seas + resid + domain + lgs + checkm2_completeness + is_mag, data = tab, phy = ph, model = "lambda")
V <- vcov(f); b <- coef(f)
dlt <- b["resid"] - b["seas"]; se <- sqrt(V["resid", "resid"] + V["seas", "seas"] - 2 * V["resid", "seas"])
out <- data.frame(test = "H6_wald_resid_minus_seasonal", n = nrow(tab), b_seasonal = b["seas"], b_resid = b["resid"],
                  diff = dlt, se_diff = se, z = dlt / se, p = 2 * pnorm(-abs(dlt / se)),
                  cor_predictors = cor(tab$seas, tab$resid))
print(out)
fwrite(out, paste0(outp, "_H6_wald.tsv"), sep = "\t")

if (length(a) > 3 && a[4] == "mcmc") {
  suppressPackageStartupMessages(library(MCMCglmm))
  set.seed(1)
  ph2 <- ph; ph2$node.label <- NULL
  ph2$edge.length[ph2$edge.length <= 0] <- 1e-8   # GTDB branch lengths kept (audit M1); no Grafen rescaling
  Ainv <- inverseA(ph2, nodes = "ALL", scale = FALSE)$Ainv
  tab$animal <- tab$accession
  pr <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
  res <- list()
  for (resp in c("n_rpf", "n_hib", "n_sr")) {
    fm <- as.formula(paste(resp, "~ P0 + domain + lgs + checkm2_completeness + is_mag"))
    m <- MCMCglmm(fm, random = ~animal, ginverse = list(animal = Ainv), family = "poisson", data = tab,
                  prior = pr, nitt = 23000, burnin = 3000, thin = 20, verbose = FALSE)
    s <- summary(m)$solutions["P0", ]
    res[[resp]] <- data.frame(response = resp, post_mean = s[1], l95 = s[2], u95 = s[3], pMCMC = s[5])
    print(res[[resp]])
  }
  fwrite(do.call(rbind, res), paste0(outp, "_mcmcglmm.tsv"), sep = "\t")
}
