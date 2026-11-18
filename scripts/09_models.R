#!/usr/bin/env Rscript
# Phylogenetic regressions for the Dormancy Rule.
# usage: Rscript 09_models.R TABLE.tsv.gz TREE.nwk SUBSET OUTPREFIX [CONFIG.tsv]
#   CONFIG defaults to the pre-registered scripts/models_config.tsv; exploratory_config.tsv holds non-pre-registered models
#   SUBSET: pilot | confirmatory | all
# Continuous/count responses: phylolm(model = "lambda") on log1p(count) (Gaussian, Pagel's lambda).
# Binary responses: phyloglm(method = "logistic_MPLE") (Ives & Garland 2010).
# Every model: response ~ predictor + domain + log10(genome_size) + checkm2_completeness + is_mag
suppressPackageStartupMessages({library(ape); library(phylolm); library(data.table)})
a <- commandArgs(TRUE)
tab <- fread(a[1]); tr <- read.tree(a[2]); subset <- a[3]; outp <- a[4]
if (subset == "pilot") tab <- tab[pilot == 1]
if (subset == "confirmatory") tab <- tab[pilot == 0]
tab <- tab[accession %in% tr$tip.label]
tab[, is_mag := as.integer(as.logical(is_mag))]
tab[, lgs := log10(genome_size)]
tab <- as.data.frame(tab); rownames(tab) <- tab$accession

# predictors are standardised within domain (z) so marine and land scales are comparable
zwithin <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
pred_raw <- c(T_logsd_total = "log_T_sd_total", T_logsd_seasonal = "log_T_sd_seasonal",
              T_logsd_resid = "log_T_sd_resid", T_ac1 = "T_ac1_resid", T_colwellP = "T_colwell_P",
              T_colwellM = "T_colwell_M", T_colwellC = "T_colwell_C",
              T_logsd_total_spmed = "log_T_sd_total_spmed", T_mean_z = "T_mean", abslat_z = "abs_lat", T_fracseas = "T_frac_seasonal", T_logsd_interann = "log_T_sd_interann",
              chl_logsd_total = "log_chl_sd_total", chl_logsd_resid = "log_chl_sd_resid",
              prec_logcv_total = "log_prec_cv_total", prec_logcv_resid = "log_prec_cv_resid", hao_lognb = "log_hao_nb", snb = "snb_genus")
tab$is_marine <- as.integer(tab$domain == "marine")
for (p in names(pred_raw)) if (pred_raw[p] %in% names(tab)) tab[[p]] <- zwithin(tab[[pred_raw[p]]], tab$domain)

cfg <- fread(if (length(a) > 4) a[5] else "scripts/models_config.tsv")  # columns: id, response, family, predictor, covars, filter
res <- list()
for (k in seq_len(nrow(cfg))) {
  r <- cfg[k]
  d <- tab
  if (nzchar(r$filter) && r$filter != "none") d <- d[with(d, eval(parse(text = r$filter))), , drop = FALSE]
  covs <- if (nzchar(r$covars) && r$covars != "none") strsplit(r$covars, "\\+")[[1]] else character(0)
  if (length(unique(d$domain)) < 2) covs <- setdiff(covs, "domain")
  vars <- unique(c(r$response, r$predictor, covs))
  vars <- vars[vars != "domain" | length(unique(d$domain)) > 1]
  d <- d[complete.cases(d[, intersect(vars, names(d)), drop = FALSE]), , drop = FALSE]
  if (nrow(d) < 30) next
  ph <- keep.tip(tr, d$accession)
  d <- d[ph$tip.label, ]
  f <- as.formula(paste(if (r$family == "gaussian_log1p") paste0("log1p(", r$response, ")") else r$response,
                        "~", paste(c(r$predictor, covs), collapse = " + ")))
  fit <- tryCatch({
    if (r$family == "binary") {
      if (length(unique(d[[r$response]])) < 2) stop("response invariant")
      phyloglm(f, data = d, phy = ph, method = "logistic_MPLE", btol = 30)
    } else phylolm(f, data = d, phy = ph, model = "lambda")
  }, error = function(e) e)
  if (inherits(fit, "error")) {
    res[[k]] <- data.frame(id = r$id, n = nrow(d), estimate = NA, se = NA, p = NA, lambda_or_alpha = NA,
                           n_pos = NA, note = conditionMessage(fit)); next
  }
  co <- summary(fit)$coefficients
  est <- co[r$predictor, 1]; se <- co[r$predictor, 2]; p <- co[r$predictor, ncol(co)]
  par <- if (r$family == "binary") fit$alpha else fit$optpar
  res[[k]] <- data.frame(id = r$id, n = nrow(d), estimate = est, se = se, p = p, lambda_or_alpha = par,
                         n_pos = if (r$family == "binary") sum(d[[r$response]]) else sum(d[[r$response]] > 0),
                         note = "")
  cat(sprintf("%-28s n=%5d  b=%+.4f  se=%.4f  p=%.3g  par=%.3f\n", r$id, nrow(d), est, se, p, par))
}
out <- do.call(rbind, res)
out$p_bh <- NA
fam <- cfg$fdr_family[match(out$id, cfg$id)]
for (g in unique(fam[!is.na(fam) & fam != "none"])) {
  i <- which(fam == g); out$p_bh[i] <- p.adjust(out$p[i], "BH")
}
fwrite(out, paste0(outp, "_models.tsv"), sep = "\t")
