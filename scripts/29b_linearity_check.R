#!/usr/bin/env Rscript
# numerical check that the procedure-P slope error and jackknife SE scale linearly with the noise (k = 1 vs 0.5).
# usage (repo root): Rscript scripts/29b_linearity_check.R  -> prints ratios (expected 0.5 and 0.5)
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table)})
tab <- as.data.frame(fread("data_small/final/analysis_table.tsv.gz")); tr <- read.tree("data_small/analysis_tree.nwk")
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
tab$P0 <- ave(tab$log_T_sd_total, tab$domain, FUN = function(v) (v - mean(v)) / sd(v)); rownames(tab) <- tab$accession
ph <- keep.tip(tr, tab$accession); tab <- tab[ph$tip.label, ]
phl <- transf.branch.lengths(ph, model = "lambda", parameters = list(lambda = 0.936))$tree
tab$cell <- paste(tab$domain, floor(tab$lat / 0.5), floor(tab$lon / 0.5)); cells <- unique(tab$cell)
set.seed(1); e <- rTrait(n = 1, phl, model = "BM", parameters = list(sigma2 = 0.01))[tab$accession]
grp <- sample(rep_len(1:20, length(cells))); names(grp) <- cells
f <- function(y) { d <- data.frame(accession = tab$accession, y = y, P0 = tab$P0, domain = tab$domain, cell = tab$cell, row.names = tab$accession)
  b <- coef(phylolm(y ~ P0 + domain, data = d, phy = phl, model = "BM"))["P0"]
  est <- sapply(1:20, function(g) { keep <- grp[d$cell] != g; p2 <- keep.tip(phl, d$accession[keep])
    coef(phylolm(y ~ P0 + domain, data = d[p2$tip.label, ], phy = p2, model = "BM"))["P0"] })
  c(b = unname(b), se = sqrt(19 / 20 * sum((est - mean(est))^2))) }
r1 <- f(-0.015 * tab$P0 + e); r2 <- f(-0.015 * tab$P0 + 0.5 * e)
cat("k=1: dev", r1["b"] + 0.015, "se", r1["se"], " k=0.5: dev", r2["b"] + 0.015, "se", r2["se"], " ratios", (r2["b"] + 0.015) / (r1["b"] + 0.015), r2["se"] / r1["se"], "\n")
