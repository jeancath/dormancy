#!/usr/bin/env python3
"""Formal equivalence tests (TOST, two one-sided tests; Schuirmann 1987; Lakens 2017) of the H1-type slopes
against the pre-registered SESOI of +-0.015 DORM_z per SD of P0.
- For Wald-type estimates (phylolm, lme4, species fixed effects): TOST p = max(p_lower, p_upper) with
  z = (beta + 0.015)/se and (0.015 - beta)/se; equivalent at alpha = 0.05 if the 90 % CI lies inside +-SESOI.
- For MCMCglmm posteriors: only the 95 % HPD is stored, so equivalence is declared if the 95 % HPD lies inside
  +-SESOI. This is conservative relative to the 90 % interval used by TOST.
usage (repo root): python3 scripts/22_equivalence_tests.py -> results/equivalence_tests.tsv"""
import pandas as pd
from scipy.stats import norm

SESOI = 0.015
rows = []
def wald(label, source, b, se, n=None):
    z_lo, z_hi = (b + SESOI) / se, (SESOI - b) / se
    p_tost = max(norm.sf(z_lo), norm.sf(z_hi))
    lo90, hi90 = b - 1.645 * se, b + 1.645 * se
    rows.append(dict(model=label, source=source, n=n, estimate=b, se=se, ci90_lo=lo90, ci90_hi=hi90,
                     p_nonzero=2 * norm.sf(abs(b / se)), p_TOST=p_tost,
                     equivalent=bool(p_tost < 0.05), method='TOST (Wald, normal)'))
def hpd(label, source, b, lo95, hi95, n=None):
    rows.append(dict(model=label, source=source, n=n, estimate=b, se=None, ci90_lo=None, ci90_hi=None,
                     ci95_lo=lo95, ci95_hi=hi95, p_nonzero=None, p_TOST=None,
                     equivalent=bool(lo95 > -SESOI and hi95 < SESOI), method='95 % HPD inside +-SESOI (conservative)'))

c = pd.read_csv('results/confirmatory_models.tsv', sep='\t').set_index('id')
for i in ['H1_DORM_P0', 'H1_DORM_P0_total', 'H6_seasonal', 'H6_resid', 'S_isolates', 'S_comp95', 'S_spmed', 'S_marine',
          'S_land', 'S_noBacCPR', 'S_growth', 'S_mobile', 'S_Tmean', 'X_chl_total', 'X_prec_cv', 'X_hao_nb', 'X_snb']:
    wald(i, 'results/confirmatory_models.tsv', c.at[i, 'estimate'], c.at[i, 'se'], int(c.at[i, 'n']))
q = pd.read_csv('audit/qc_fix_rerun.tsv', sep='\t')
for _, r in q[(q.set == 'qc_fixed') & q.model.str.startswith('H1')].iterrows():
    wald(f'{r.model} (QC-corrected)', 'audit/qc_fix_rerun.tsv', r.estimate, r.se, int(r.n))
s = pd.read_csv('results/siteMCMC_mcmcglmm_site.tsv', sep='\t').iloc[0]
hpd('H1 pooled, MCMCglmm phylogeny + site (GTDB branch lengths)', 'results/siteMCMC_mcmcglmm_site.tsv',
    s.post_mean, s.l95, s.u95, int(s.n))
w = pd.read_csv('replication_within_species/results.tsv', sep='\t')
w = w[(w['sample'] == 'hq') & (w.model == 'a_speciesFE_clsite') & (w.response == 'DORM_z') & (w.term == 'P0')].iloc[0]
wald('Within-species pooled slope (species FE, site-clustered)', 'replication_within_species/results.tsv',
     w.estimate, w.se, int(w.n_genomes))
k = pd.read_csv('results/siteK1_cluster_robust.tsv', sep='\t').set_index('set')
for st in ['H1_direct', 'H1_total']:
    wald(f'{st} with site-grouped jackknife SE (most conservative)', 'results/siteK1_cluster_robust.tsv',
         k.at[st, 'beta'], k.at[st, 'se_site_jackknife'], int(k.at[st, 'n']))
out = pd.DataFrame(rows)
out.to_csv('results/equivalence_tests.tsv', sep='\t', index=False, float_format='%.6g')
print(out[['model', 'n', 'estimate', 'ci90_lo', 'ci90_hi', 'p_TOST', 'equivalent']].round(4).to_string())
