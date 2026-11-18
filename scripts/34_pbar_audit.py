#!/usr/bin/env python3
"""Added 18 Nov 2026: audit of the final procedure P-bar.

Definition. Let b be the GLS slope of P0 (lambda fixed at the full-data estimate). For a random assignment k of the
C = 1,642 occupied 0.5-degree cells to G = 20 groups, let b_(k,g) be the slope with all species of group g removed and
b_(k,.) their mean. The single-grouping jackknife variance is
    V_k = (G - 1) / G * sum_g (b_(k,g) - b_(k,.))^2 .
P-bar averages VARIANCES over K independent assignments, then takes the root:
    SE_bar = sqrt( (1/K) * sum_k V_k ),    K = 200 on the real data (5 in the calibration simulations),
and declares equivalence at margin M = 0.015 iff  b - t_{0.95,19} SE_bar > -M  and  b + t_{0.95,19} SE_bar < M
(p_TOST = max of the two one-sided t_19 p-values). Averaging SEs instead, (1/K) sum_k sqrt(V_k), gives a slightly
smaller value (Jensen); both are reported.

Outputs (repo root):
  results/pbar_se_distribution.tsv      SE over the 200 real-data groupings: quantiles, both averages, decisions
  results/pbar_calibration_by_boundary.tsv  false-equivalence rate per DGM, per boundary, per K, Clopper-Pearson CIs;
                                        K = 1, K = 5 and an 'oracle' SE (the simulation SD of the slope) bracket K = 200
  results/pbar_estimator_family.tsv     the decision under every SE estimator considered; the intersection-union rule
  results/figures/fig_pbar_se_distribution.png
usage: python3 scripts/34_pbar_audit.py"""
import numpy as np, pandas as pd
from scipy.stats import t as T, beta as B
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
M, G = 0.015, 20; q = T.ppf(0.95, G - 1)
def cp(k, n): return (B.ppf(.025, k, n - k + 1) if k else 0.0, B.ppf(.975, k + 1, n - k) if k < n else 1.0)
def ptost(b, se, df): return max(T.sf((b + M) / se, df), T.sf((M - b) / se, df))
# ---- 1. SE distribution on the real data
runs = pd.read_csv('results/grouping_sensitivity_runs.tsv', sep='\t')
rows = []
for s, g in runs[runs.variant == 'random20'].groupby('set', sort=False):
    b = g.beta.iloc[0]; se = g.se.to_numpy()
    sev, ses = np.sqrt(np.mean(se ** 2)), np.mean(se)
    rows.append(dict(set=s, n_groupings=len(se), beta=b, se_min=se.min(), se_q05=np.quantile(se, .05), se_q25=np.quantile(se, .25),
                     se_median=np.median(se), se_q75=np.quantile(se, .75), se_q95=np.quantile(se, .95), se_max=se.max(),
                     cv_of_se=se.std(ddof=1) / se.mean(), se_original_grouping=se[0], rank_of_original=int((se < se[0]).sum()) + 1,
                     SE_bar_mean_of_variances=sev, SE_mean_of_SEs=ses,
                     p_tost_SE_bar=ptost(b, sev, 19), p_tost_mean_SE=ptost(b, ses, 19),
                     share_single_groupings_equivalent=float(np.mean((b - q * se > -M) & (b + q * se < M)))))
d1 = pd.DataFrame(rows); d1.to_csv('results/pbar_se_distribution.tsv', sep='\t', index=False, float_format='%.6g')
print(d1.round(5).T.to_string())
h = runs[(runs.variant == 'random20') & (runs.set == 'H1_direct')]
fig, ax = plt.subplots(figsize=(6.2, 3.4), dpi=150)
ax.hist(h.se, bins=30, color='#7f8fa6', edgecolor='white')
for x, lab, c in [(h.se.iloc[0], 'original grouping', '#c0392b'), (d1.SE_bar_mean_of_variances.iloc[0], 'SE-bar (200)', '#1f6f8b'),
                  (runs[(runs.set == 'H1_direct') & (runs.variant == 'one_cell')].se.iloc[0], 'delete-one-cell', '#2d8a4e')]:
    ax.axvline(x, color=c, lw=1.8, label=lab)
ax.axvline((M - h.beta.iloc[0]) / q, color='k', ls='--', lw=1, label='SE above which\nequivalence fails')
ax.set_xlabel('jackknife SE of the pooled slope (one random grouping of cells)'); ax.set_ylabel('groupings')
ax.legend(fontsize=7, frameon=False); fig.tight_layout(); fig.savefig('results/figures/fig_pbar_se_distribution.png'); plt.close(fig)
# ---- 2. calibration by boundary
cal = []
old = pd.read_csv('results/tost_calibration_sims.tsv', sep='\t'); old['dgm'] = 'A original, lambda fixed (first calibration)'
lf = pd.read_csv('results/calibration_lambda_free_sims.tsv', sep='\t'); lf['dgm'] = lf.dgm.map({1: 'A original, lambda re-estimated', 2: 'B single fitted model, lambda re-estimated'})
rs = pd.read_csv('results/calibration_real_structure_sims.tsv', sep='\t'); rs['dgm'] = 'D real residual structure, lambda re-estimated'
for df in (old, lf, rs):
    for (dg, tb), g in df[df.true_beta.abs() > 0.0149].groupby(['dgm', 'true_beta']):
        n = len(g); e = g.estimate.to_numpy(); sd = (e - tb).std(ddof=1)
        dec = {'K = 1 (one grouping)': g.equiv.to_numpy()}
        if 'equiv_avg' in g and g.equiv_avg.notna().all(): dec['K = 5'] = g.equiv_avg.to_numpy()
        dec['oracle SE (simulation SD of the slope; limit of large K)'] = ((e - q * sd > -M) & (e + q * sd < M)).astype(int)
        for kk, v in dec.items():
            k = int(np.sum(v)); lo, hi = cp(k, n)
            cal.append(dict(dgm=dg, boundary=tb, K=kk, nsim=n, n_equivalent=k, rate=k / n, mc_lo95=lo, mc_hi95=hi,
                            sd_slope=sd, mean_se_one_grouping=g.se.mean()))
d2 = pd.DataFrame(cal); d2.to_csv('results/pbar_calibration_by_boundary.tsv', sep='\t', index=False, float_format='%.5g')
pd.set_option('display.width', 220); print(d2[['dgm', 'boundary', 'K', 'nsim', 'rate', 'mc_lo95', 'mc_hi95']].round(4).to_string())
# ---- 3. estimator family and the intersection-union rule (pooled slope)
ft = pd.read_csv('results/genome_source_final_table.tsv', sep='\t')
pool = ft[ft.quantity == 'pooled slope (H1 direct)'].copy()
sb = pd.read_csv('results/spatial_blocks_jackknife.tsv', sep='\t')
fam = [dict(estimator=r.procedure, role=r.role, se=r.se, p_tost=r.p_tost, equivalent=bool(r.equivalent_at_0015)) for _, r in pool.iterrows()]
for _, r in sb[sb.set == 'H1_direct'].iterrows():
    fam.append(dict(estimator=f'block jackknife, {r.block_deg:g}-degree grid blocks ({int(r.n_blocks)} blocks in 100 groups; scripts/25)', role='sensitivity (robustness)',
                    se=r.se_block_jackknife, p_tost=r.p_TOST, equivalent=bool(r.ci90_lo > -M and r.ci90_hi < M)))
gs1 = runs[(runs.variant == 'random20') & (runs.set == 'H1_direct')]
fam.append(dict(estimator='single random groupings (200 realisations of P; not an estimator choice)', role='randomisation of P',
                se=float(gs1.se.median()), p_tost=float(gs1.p_tost.max()), equivalent=bool(gs1.equiv.all())))
d3 = pd.DataFrame(fam)
est = d3[~d3.role.str.startswith(('superseded', 'randomisation'))]
iu = dict(estimator='INTERSECTION-UNION: equivalence only if every non-superseded estimator above declares it', role='decision rule',
          se=np.nan, p_tost=est.p_tost.max(), equivalent=bool(est.equivalent.all()))
d3 = pd.concat([d3, pd.DataFrame([iu])], ignore_index=True)
d3.to_csv('results/pbar_estimator_family.tsv', sep='\t', index=False, float_format='%.5g')
print(d3.to_string())
