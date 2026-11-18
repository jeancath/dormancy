#!/usr/bin/env Rscript
# Added 18 Nov 2026: how much does procedure P on the REAL data depend on (a) the random assignment of 0.5-degree
# cells to the 20 jackknife groups, (b) the groups being geographically scattered rather than contiguous, and (c) lambda
# being held fixed?  Procedure P (scripts/26_tost_calibration.R) assigns the 1,642 floor-aligned 0.5-degree cells at
# random to 20 groups (about 82 cells each, scattered across the globe); the deletion unit is a whole group of cells.
# Variants, all with the same point estimate and the same TOST margin (+-0.015, 90 % t-interval):
#   random20      200 independent random assignments of cells to 20 groups (lambda fixed)
#   one_cell      delete-one-cell jackknife over every cell (no grouping randomness; df = cells - 1; lambda fixed)
#   contiguous20  20 spatially contiguous groups from k-means on cell centroids (unit-sphere coordinates), 20 k-means
#                 seeds (lambda fixed)
#   lambda_free   the original grouping (seed 20261118), lambda re-estimated in every jackknife replicate
# usage: Rscript 30_grouping_sensitivity.R TABLE.tsv.gz TREE.nwk OUT_PREFIX [cores=4]
#   -> OUT_PREFIX_runs.tsv (one row per variant x replicate), OUT_PREFIX_summary.tsv
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table); library(parallel)})
a <- commandArgs(TRUE); cores <- if (length(a) > 3) as.integer(a[4]) else 4
G <- 20; SESOI <- 0.015
tab <- as.data.frame(fread(a[1])); tr <- read.tree(a[2])
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
tab$P0 <- z(tab$log_T_sd_total, tab$domain); tab$lgs <- log10(tab$genome_size)
tab$is_mag <- as.integer(as.logical(tab$is_mag)); rownames(tab) <- tab$accession
ph <- keep.tip(tr, tab$accession); tab <- tab[ph$tip.label, ]
tab$cell <- paste(tab$domain, floor(tab$lat / 0.5), floor(tab$lon / 0.5))
full <- phylolm(DORM_z ~ P0 + domain + lgs + checkm2_completeness + is_mag, data = tab, phy = ph, model = "lambda")
lam <- full$optpar
phl <- transf.branch.lengths(ph, model = "lambda", parameters = list(lambda = lam))$tree
cells <- unique(tab$cell)
cc <- aggregate(cbind(lat, lon) ~ cell, data = tab, FUN = mean)
xyz <- with(cc, cbind(cos(lat * pi / 180) * cos(lon * pi / 180), cos(lat * pi / 180) * sin(lon * pi / 180), sin(lat * pi / 180)))
rownames(xyz) <- cc$cell
sets <- list(H1_direct = list(d = tab, f = DORM_z ~ P0 + domain + lgs + checkm2_completeness + is_mag),
             MAGs = list(d = tab[tab$is_mag == 1, ], f = DORM_z ~ P0 + domain + lgs + checkm2_completeness),
             isolates = list(d = tab[tab$is_mag == 0, ], f = DORM_z ~ P0 + domain + lgs + checkm2_completeness))
fitBM <- function(d, f, tree) coef(phylolm(f, data = d[tree$tip.label, ], phy = tree, model = "BM"))["P0"]
summ1 <- function(set, variant, rep, b, est, nG) {
  se <- sqrt((nG - 1) / nG * sum((est - mean(est))^2)); q <- qt(0.95, nG - 1)
  data.frame(set = set, variant = variant, rep = rep, G = nG, beta = b, se = se, lo90 = b - q * se, hi90 = b + q * se,
             lo95 = b - qt(0.975, nG - 1) * se, hi95 = b + qt(0.975, nG - 1) * se,
             p_tost = max(pt((b + SESOI) / se, nG - 1, lower.tail = FALSE), pt((SESOI - b) / se, nG - 1, lower.tail = FALSE)),
             equiv = as.integer(b - q * se > -SESOI & b + q * se < SESOI),
             p_gt0 = pt(b / se, nG - 1, lower.tail = FALSE), p_gt_sesoi = pt((b - SESOI) / se, nG - 1, lower.tail = FALSE))
}
out <- list()
for (s in names(sets)) {
  d <- sets[[s]]$d; f <- sets[[s]]$f; tp <- keep.tip(phl, d$accession); d <- d[tp$tip.label, ]
  b <- fitBM(d, f, tp); dc <- unique(d$cell)
  jk <- function(grp) sapply(sort(unique(grp)), function(g) { keep <- grp[d$cell] != g; fitBM(d[keep, ], f, keep.tip(tp, d$accession[keep])) })
  # (a) 200 random groupings (the first uses the original seed for the full-data set)
  r <- mclapply(1:200, function(i) {
    set.seed(if (i == 1) 20261118 else 1000 + i)
    grp <- sample(rep_len(seq_len(G), length(cells))); names(grp) <- cells
    grp <- grp[dc]; summ1(s, "random20", i, b, jk(grp), length(unique(grp)))
  }, mc.cores = cores)
  out[[length(out) + 1]] <- do.call(rbind, r); cat(s, "random20 done\n")
  # (b) delete-one-cell
  est <- unlist(mclapply(dc, function(cl) { keep <- d$cell != cl; fitBM(d[keep, ], f, keep.tip(tp, d$accession[keep])) }, mc.cores = cores))
  out[[length(out) + 1]] <- summ1(s, "one_cell", 1, b, est, length(dc)); cat(s, "one_cell done\n")
  # (c) contiguous k-means groups
  r <- mclapply(1:20, function(i) {
    set.seed(i); km <- kmeans(xyz[dc, , drop = FALSE], centers = G, nstart = 5, iter.max = 100)
    summ1(s, "contiguous20", i, b, jk(km$cluster), G)
  }, mc.cores = cores)
  out[[length(out) + 1]] <- do.call(rbind, r); cat(s, "contiguous20 done\n")
  # (d) lambda re-estimated in every replicate (original grouping), full data set only
  if (s == "H1_direct") {
    set.seed(20261118); grp <- sample(rep_len(seq_len(G), length(cells))); names(grp) <- cells
    est <- unlist(mclapply(seq_len(G), function(g) { keep <- grp[d$cell] != g; dd <- d[keep, ]; p2 <- keep.tip(ph, dd$accession)
      coef(phylolm(f, data = dd[p2$tip.label, ], phy = p2, model = "lambda"))["P0"] }, mc.cores = cores))
    out[[length(out) + 1]] <- summ1(s, "lambda_free", 1, coef(full)["P0"], est, G); cat(s, "lambda_free done\n")
  }
  fwrite(do.call(rbind, out), paste0(a[3], "_runs.tsv"), sep = "\t")
}
runs <- rbindlist(out)
sm <- runs[, .(n_runs = .N, beta = beta[1], se_median = median(se), se_q05 = quantile(se, 0.05), se_q95 = quantile(se, 0.95),
               p_tost_median = median(p_tost), share_equivalent = mean(equiv), p_gt_sesoi_median = median(p_gt_sesoi)),
           by = .(set, variant)]
fwrite(sm, paste0(a[3], "_summary.tsv"), sep = "\t"); print(sm)
