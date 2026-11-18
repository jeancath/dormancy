#!/usr/bin/env Rscript
# Prune GTDB r232 bac120 species tree to the species in the analysis set; tips relabelled to the
# analysis genome accession (one genome per species). Branch lengths: GTDB substitutions/site.
suppressPackageStartupMessages(library(ape))
a <- read.delim(gzfile("data_small/analysis_set.tsv.gz"))
tr <- read.tree("/home/user/data/gtdb/bac120.tree")
cat("tree tips:", Ntip(tr), "\n")
m <- match(a$gtdb_genome_representative, tr$tip.label)
cat("analysis species with rep on tree:", sum(!is.na(m)), "of", nrow(a), "\n")
keep <- a[!is.na(m), ]
pt <- keep.tip(tr, keep$gtdb_genome_representative)
pt$tip.label <- keep$accession[match(pt$tip.label, keep$gtdb_genome_representative)]
pt$node.label <- NULL
# zero-length branches break some GLS fits; give them a tiny length
pt$edge.length[pt$edge.length <= 0] <- 1e-6
if (!is.binary(pt)) pt <- multi2di(pt, random = FALSE)
pt$edge.length[pt$edge.length <= 0] <- 1e-6
write.tree(pt, "data_small/analysis_tree.nwk")
cat("pruned tips:", Ntip(pt), " ultrametric:", is.ultrametric(pt), "\n")
