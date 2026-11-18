#!/usr/bin/env Rscript
# step 4: within-species models (EXPLORATORY replication of an exploratory finding).
# usage (repo root): Rscript replication_within_species/04_models.R TABLE.tsv.gz OUTPREFIX [label]
suppressPackageStartupMessages({library(data.table); library(lme4); library(lmerTest); library(sandwich); library(lmtest)})
a <- commandArgs(TRUE); outp <- a[2]; lab <- if (length(a) > 2) a[3] else "main"
d <- as.data.frame(fread(a[1]))
d$is_mag <- as.integer(d$is_mag); d$species <- factor(d$species); d$site <- factor(d$site)
d$comp <- d$checkm2_completeness
# Mundlak split: within-species deviation and species mean of P0
d$P0_sp <- ave(d$P0, d$species); d$P0_w <- d$P0 - d$P0_sp
d$hib_systems <- as.numeric(d$hib_systems)
resp <- c(DORM_z = "DORM_z", log_n_hib = "log1p(n_hib)", hib_extra = "log1p(hib_extra)",
          log_hpf_short = "log1p(hpf_short)", log_n_balon = "log1p(n_balon)", log_n_sr = "log1p(n_sr)",
          log_n_sas = "log1p(n_sas)", log_spo = "log1p(spo_strict)", log_rpf = "log1p(n_rpf)",
          log_ta = "log1p(ta_pairs)")
res <- list()
add <- function(model, response, term, est, se, p, n, nsp, note = "") {
  res[[length(res) + 1]] <<- data.frame(sample = lab, model = model, response = response, term = term, estimate = est,
                                         se = se, p = p, n_genomes = n, n_species = nsp, note = note)
}
cv <- "lgs + comp + is_mag"
# descriptive: within-species SD of P0
sdw <- sd(d$P0_w); sdb <- sd(tapply(d$P0, d$species, mean))
add("descriptive", "P0", "within_species_SD", sdw, NA, NA, nrow(d), nlevels(d$species),
    sprintf("between-species SD of species-mean P0 = %.3f; species with >1 site = %d", sdb, sum(tapply(d$site, d$species, function(v) length(unique(v))) > 1)))
cr <- function(fit, cl) coeftest(fit, vcov. = vcovCL(fit, cluster = cl, type = "HC1"))
for (rn in names(resp)) {
  y <- resp[[rn]]
  # (a) random-slope mixed model, Mundlak within/between P0
  f1 <- as.formula(sprintf("%s ~ P0_w + P0_sp + %s + (1 + P0_w | species)", y, cv))
  m1 <- tryCatch(lmer(f1, data = d, control = lmerControl(optimizer = "bobyqa")), error = function(e) NULL)
  if (!is.null(m1)) { co <- summary(m1)$coefficients
    add("a_lmer_randslope", rn, "P0_within", co["P0_w", 1], co["P0_w", 2], co["P0_w", 5], nrow(d), nlevels(d$species),
        if (isSingular(m1)) "singular fit" else "")
    add("a_lmer_randslope", rn, "P0_between", co["P0_sp", 1], co["P0_sp", 2], co["P0_sp", 5], nrow(d), nlevels(d$species)) }
  # (a) as specified in the replication plan: raw P0 with (1 + P0 | species)
  f1b <- as.formula(sprintf("%s ~ P0 + %s + (1 + P0 | species)", y, cv))
  m1b <- tryCatch(lmer(f1b, data = d, control = lmerControl(optimizer = "bobyqa")), error = function(e) NULL)
  if (!is.null(m1b)) { co <- summary(m1b)$coefficients
    add("a_lmer_spec", rn, "P0", co["P0", 1], co["P0", 2], co["P0", 5], nrow(d), nlevels(d$species), if (isSingular(m1b)) "singular fit" else "") }
  # (a) species fixed effects; SEs clustered by species and (e) by site
  f2 <- as.formula(sprintf("%s ~ P0 + %s + species", y, cv)); m2 <- lm(f2, data = d)
  for (cl in c("species", "site")) { co <- cr(m2, d[[cl]])
    add(paste0("a_speciesFE_cl", cl), rn, "P0", co["P0", 1], co["P0", 2], co["P0", 4], nrow(d), nlevels(d$species)) }
  # (b) within-species MAG/isolate reversal: P0 x is_mag with species FE
  f3 <- as.formula(sprintf("%s ~ P0 * is_mag + lgs + comp + species", y)); m3 <- lm(f3, data = d)
  for (cl in c("species", "site")) { co <- cr(m3, d[[cl]])
    add(paste0("b_speciesFE_cl", cl), rn, "P0:is_mag", co["P0:is_mag", 1], co["P0:is_mag", 2], co["P0:is_mag", 4], nrow(d), nlevels(d$species))
    add(paste0("b_speciesFE_cl", cl), rn, "P0 (isolates)", co["P0", 1], co["P0", 2], co["P0", 4], nrow(d), nlevels(d$species)) }
  # slope in MAGs = P0 + P0:is_mag (site-clustered)
  V <- vcovCL(m3, cluster = d$site, type = "HC1"); b <- coef(m3)
  g <- c(P0 = 1, `P0:is_mag` = 1); est <- sum(b[names(g)]); se <- sqrt(drop(t(g) %*% V[names(g), names(g)] %*% g))
  add("b_speciesFE_clsite", rn, "P0 (MAGs)", est, se, 2 * pnorm(-abs(est / se)), nrow(d), nlevels(d$species))
  # (b) random-slope version of the interaction, with (e) site random intercept
  f4 <- as.formula(sprintf("%s ~ P0_w * is_mag + P0_sp + lgs + comp + (1 + P0_w | species) + (1 | site)", y))
  m4 <- tryCatch(lmer(f4, data = d, control = lmerControl(optimizer = "bobyqa")), error = function(e) NULL)
  if (!is.null(m4)) { co <- summary(m4)$coefficients
    add("b_e_lmer_randslope_site", rn, "P0_within:is_mag", co["P0_w:is_mag", 1], co["P0_w:is_mag", 2], co["P0_w:is_mag", 5], nrow(d), nlevels(d$species), if (isSingular(m4)) "singular fit" else "")
    add("b_e_lmer_randslope_site", rn, "P0_within (isolates)", co["P0_w", 1], co["P0_w", 2], co["P0_w", 5], nrow(d), nlevels(d$species)) }
  # (e) site random intercept added to the (a) model
  f5 <- as.formula(sprintf("%s ~ P0_w + P0_sp + %s + (1 + P0_w | species) + (1 | site)", y, cv))
  m5 <- tryCatch(lmer(f5, data = d, control = lmerControl(optimizer = "bobyqa")), error = function(e) NULL)
  if (!is.null(m5)) { co <- summary(m5)$coefficients
    add("e_lmer_randslope_site", rn, "P0_within", co["P0_w", 1], co["P0_w", 2], co["P0_w", 5], nrow(d), nlevels(d$species), if (isSingular(m5)) "singular fit" else "") }
  # genome-type-specific within-species slopes (species FE within each type)
  for (g in c(1, 0)) { dd <- d[d$is_mag == g, ]; dd$species <- droplevels(dd$species)
    dd <- dd[dd$species %in% names(which(table(dd$species) >= 2)), ]; dd$species <- droplevels(dd$species)
    if (nrow(dd) < 30) next
    m6 <- lm(as.formula(sprintf("%s ~ P0 + lgs + comp + species", y)), data = dd); co <- cr(m6, dd$site)
    add("c_speciesFE_within_type_clsite", rn, ifelse(g == 1, "P0 (MAG-only species sets)", "P0 (isolate-only species sets)"),
        co["P0", 1], co["P0", 2], co["P0", 4], nrow(dd), nlevels(dd$species)) }
  cat(rn, "done\n")
}
# (d) mixed species: paired MAG vs isolate comparison
mx <- d[d$mixed_species == 1, ]; mx$species <- droplevels(mx$species)
if (nrow(mx) > 0) {
agg <- do.call(rbind, lapply(split(mx, mx$species), function(s) {
  M <- s[s$is_mag == 1, ]; I <- s[s$is_mag == 0, ]
  data.frame(species = s$species[1], domain = s$domain[1], nM = nrow(M), nI = nrow(I),
             dP0 = mean(M$P0) - mean(I$P0), dDORM = mean(M$DORM_z) - mean(I$DORM_z),
             dhibx = mean(M$hib_extra) - mean(I$hib_extra), dlhibx = mean(log1p(M$hib_extra)) - mean(log1p(I$hib_extra)),
             dhib = mean(M$n_hib) - mean(I$n_hib), dcomp = mean(M$comp) - mean(I$comp), dlgs = mean(M$lgs) - mean(I$lgs),
             dsr = mean(M$n_sr) - mean(I$n_sr), dta = mean(log1p(M$ta_pairs)) - mean(log1p(I$ta_pairs)))
}))
fwrite(agg, paste0(outp, "_mixed_pairs.tsv"), sep = "\t")
ns <- nrow(agg); nmx <- nrow(mx)
for (v in c("dDORM", "dlhibx", "dhibx", "dhib", "dsr", "dta")) {
  tt <- t.test(agg[[v]]); add("d_paired_mean_MAG_minus_isolate", v, "mean difference", unname(tt$estimate), unname(tt$stderr), tt$p.value, nmx, ns,
                              sprintf("sign test: %d of %d species positive", sum(agg[[v]] > 0), sum(agg[[v]] != 0)))
  f <- lm(as.formula(paste(v, "~ dP0 + dcomp + dlgs")), data = agg); co <- coeftest(f, vcov. = vcovHC(f, type = "HC1"))
  add("d_paired_diff_on_dP0", v, "dP0", co["dP0", 1], co["dP0", 2], co["dP0", 4], nmx, ns, "Δ(MAG−isolate) ~ ΔP0 + Δcompleteness + Δlog genome size, HC1 SE")
}
# genome-level version within mixed species: species FE, is_mag main effect = within-species MAG-isolate gap
for (rn in c("DORM_z", "hib_extra", "log_n_hib", "log_hpf_short", "log_n_balon", "log_n_sr")) {
  y <- resp[[rn]]
  yv <- eval(parse(text = y), mx); wv <- yv - ave(yv, mx$species)
  if (max(abs(wv)) < 1e-9) {  # response constant within every mixed species: not estimable
    for (tm in c("is_mag (MAG-isolate gap at P0=0)", "P0:is_mag")) add("d_mixed_speciesFE_clsite", rn, tm, NA, NA, NA, nmx, ns, "no within-species variation")
    next
  }
  m7 <- lm(as.formula(sprintf("%s ~ is_mag * P0 + lgs + comp + species", y)), data = mx); co <- cr(m7, mx$site)
  add("d_mixed_speciesFE_clsite", rn, "is_mag (MAG-isolate gap at P0=0)", co["is_mag", 1], co["is_mag", 2], co["is_mag", 4], nmx, ns)
  add("d_mixed_speciesFE_clsite", rn, "P0:is_mag", co["is_mag:P0", 1], co["is_mag:P0", 2], co["is_mag:P0", 4], nmx, ns)
}
}
out <- do.call(rbind, res); fwrite(out, paste0(outp, "_results.tsv"), sep = "\t"); print(out[, c(2:7)], digits = 3)
