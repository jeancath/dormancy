#!/usr/bin/env Rscript
# Spatial-block sensitivity (added 18 Nov 2026): does the conclusion depend on defining a
# "site" as a 0.5-degree cell? For block sizes of 0.5, 1, 2 and 5 degrees:
#   (a) delete-a-group jackknife over blocks (G = 100 groups) for H1 direct, MAGs and isolates -> SE, TOST
#   (b) one random species per block, R resamples -> distribution (quantiles) of the H1 slope.
# usage: Rscript 25_spatial_blocks.R TABLE.tsv.gz TREE.nwk OUT_PREFIX [R=100]
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table)})
a <- commandArgs(TRUE); R <- if (length(a) > 3) as.integer(a[4]) else 100
set.seed(20261118)
tab <- as.data.frame(fread(a[1])); tr <- read.tree(a[2])
tab <- tab[tab$pilot == 0 & tab$accession %in% tr$tip.label, ]
z <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
tab$P0 <- z(tab$log_T_sd_total, tab$domain); tab$lgs <- log10(tab$genome_size)
tab$is_mag <- as.integer(as.logical(tab$is_mag)); rownames(tab) <- tab$accession
fitb <- function(d, covs) {
  ph <- keep.tip(tr, d$accession); d <- d[ph$tip.label, ]
  if (length(unique(d$domain)) < 2) covs <- setdiff(covs, "domain")
  co <- summary(phylolm(as.formula(paste("DORM_z ~", paste(c("P0", covs), collapse = " + "))), data = d, phy = ph,
                        model = "lambda"))$coefficients
  co["P0", 1:2]
}
sets <- list(H1_direct = list(f = function(d) d, c = c("domain", "lgs", "checkm2_completeness", "is_mag")),
             MAGs = list(f = function(d) d[d$is_mag == 1, ], c = c("domain", "lgs", "checkm2_completeness")),
             isolates = list(f = function(d) d[d$is_mag == 0, ], c = c("domain", "lgs", "checkm2_completeness")))
jk <- list(); dist <- list()
for (bs in c(0.5, 1, 2, 5)) {
  tab$block <- paste(tab$domain, floor(tab$lat / bs), floor(tab$lon / bs))
  for (nm in names(sets)) {
    d0 <- sets[[nm]]$f(tab); cv <- sets[[nm]]$c; full <- fitb(d0, cv)
    bl <- unique(d0$block); G <- min(100, length(bl)); grp <- sample(rep_len(seq_len(G), length(bl))); names(grp) <- bl
    est <- sapply(seq_len(G), function(g) fitb(d0[grp[d0$block] != g, ], cv)[1])
    se <- sqrt((G - 1) / G * sum((est - mean(est))^2))
    p_tost <- max(pnorm((full[1] + 0.015) / se, lower.tail = FALSE), pnorm((0.015 - full[1]) / se, lower.tail = FALSE))
    jk[[length(jk) + 1]] <- data.frame(block_deg = bs, set = nm, n = nrow(d0), n_blocks = length(bl), beta = full[1],
                                       se_phylolm = full[2], se_block_jackknife = se, ci90_lo = full[1] - 1.645 * se,
                                       ci90_hi = full[1] + 1.645 * se, p_TOST = p_tost)
    print(tail(jk, 1)[[1]])
  }
  if (bs >= 1) {
    b <- sapply(seq_len(R), function(r) {
      idx <- unlist(lapply(split(seq_len(nrow(tab)), tab$block), function(i) if (length(i) > 1) sample(i, 1) else i))
      fitb(tab[idx, ], sets$H1_direct$c)[1]
    })
    dist[[length(dist) + 1]] <- data.frame(block_deg = bs, n_blocks = length(unique(tab$block)), R = R,
      median = median(b), q025 = quantile(b, 0.025), q25 = quantile(b, 0.25), q75 = quantile(b, 0.75),
      q975 = quantile(b, 0.975), frac_positive = mean(b > 0), frac_inside_SESOI = mean(abs(b) < 0.015))
    print(tail(dist, 1)[[1]])
  }
  fwrite(do.call(rbind, jk), paste0(a[3], "_jackknife.tsv"), sep = "\t")
  if (length(dist)) fwrite(do.call(rbind, dist), paste0(a[3], "_one_per_block.tsv"), sep = "\t")
}
