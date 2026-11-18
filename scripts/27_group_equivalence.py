#!/usr/bin/env python3
"""Genome-source-specific slopes, one-sided tests against 0 and against the SESOI , formal interaction tests, within-group equivalence, and equivalence-margin
sensitivity (added 18 Nov 2026). All rows are assembled from committed result files.
usage (repo root): python3 scripts/27_group_equivalence.py -> results/group_equivalence.tsv, results/margin_sensitivity.tsv"""
import numpy as np, pandas as pd
from scipy.stats import norm, t as tdist

def row(target, method, source, b, se=None, lo95=None, hi95=None, df=None, n=None):
    if se is not None:
        q90 = tdist.ppf(0.95, df) if df else norm.ppf(0.95); q95 = tdist.ppf(0.975, df) if df else norm.ppf(0.975)
        lo90, hi90, lo95, hi95 = b - q90 * se, b + q90 * se, b - q95 * se, b + q95 * se
    else:
        lo90 = hi90 = np.nan
    out = dict(target=target, method=method, source=source, n=n, estimate=b, se=se, ci90_lo=lo90, ci90_hi=hi90,
               ci95_lo=lo95, ci95_hi=hi95)
    # one-sided tests: is the slope above 0, and above the SESOI (H0: beta <= 0.015)?
    if se is not None:
        sf = (lambda z: tdist.sf(z, df)) if df else norm.sf
        out['p_two_sided_0'] = 2 * sf(abs(b) / se); out['p_one_sided_gt0'] = sf(b / se)
        out['p_one_sided_gt_0.015'] = sf((b - 0.015) / se); out['p_one_sided_lt_-0.015'] = sf((-0.015 - b) / se)
    for m in (0.010, 0.015, 0.020):
        if se is not None:
            sf = (lambda z: tdist.sf(z, df)) if df else norm.sf
            p = max(sf((b + m) / se), sf((m - b) / se)); out[f'p_TOST_{m:.3f}'] = p
            out[f'equiv_{m:.3f}'] = bool(lo90 > -m and hi90 < m)
        else:
            out[f'p_TOST_{m:.3f}'] = np.nan; out[f'equiv_{m:.3f}'] = bool(lo95 > -m and hi95 < m)
    return out

R = []
c = pd.read_csv('results/confirmatory_models.tsv', sep='\t').set_index('id')
e = pd.read_csv('results/exploratory_models.tsv', sep='\t').set_index('id')
k1 = pd.read_csv('results/siteK1_cluster_robust.tsv', sep='\t').set_index('set')
P = pd.read_csv('results/tost_calibration_real.tsv', sep='\t').set_index('set')
mc_int = pd.read_csv('results/exploratory_mcmc_site_interaction.tsv', sep='\t')
mc_pool = pd.read_csv('results/siteMCMC_mcmcglmm_site.tsv', sep='\t').iloc[0]
ws = pd.read_csv('replication_within_species/results.tsv', sep='\t'); ws = ws[ws['sample'] == 'hq']
def wsrow(model, term):
    x = ws[(ws.model == model) & (ws.response == 'DORM_z') & (ws.term == term)].iloc[0]; return x.estimate, x.se, int(x.n_genomes)

for tgt, cid, kset in (('pooled (H1 direct)', 'H1_DORM_P0', 'H1_direct'),):
    R.append(row(tgt, 'phylolm lambda, model SE', 'confirmatory_models', c.at[cid, 'estimate'], c.at[cid, 'se'], n=int(c.at[cid, 'n'])))
    R.append(row(tgt, 'site jackknife, G=100 (0.5-deg centred cells)', 'siteK1_cluster_robust', k1.at[kset, 'beta'], k1.at[kset, 'se_site_jackknife'], n=int(k1.at[kset, 'n'])))
    R.append(row(tgt, 'PRIMARY procedure P: lambda fixed, jackknife G=20 (0.5-deg grid cells), t19', 'tost_calibration_real', P.at[kset, 'beta'], P.at[kset, 'se'], df=19, n=int(P.at[kset, 'n'])))
    R.append(row(tgt, 'MCMCglmm phylogeny + site (95% HPD)', 'siteMCMC_mcmcglmm_site', mc_pool.post_mean, lo95=mc_pool.l95, hi95=mc_pool.u95, n=int(mc_pool.n)))
    b, s, n = wsrow('a_speciesFE_clsite', 'P0'); R.append(row('pooled, WITHIN species', 'species fixed effects, site-clustered SE', 'replication_within_species', b, s, n=n))
for grp, eid, kset, mcm in (('MAG-represented species', 'E_MAGs_only', 'MAGs', 'MAGs'), ('isolate-represented species', None, 'isolates', 'isolates')):
    if eid:
        R.append(row(grp, 'phylolm lambda, model SE', 'exploratory_models', e.at[eid, 'estimate'], e.at[eid, 'se'], n=int(e.at[eid, 'n'])))
    else:
        R.append(row(grp, 'phylolm lambda, model SE', 'confirmatory_models (S_isolates)', c.at['S_isolates', 'estimate'], c.at['S_isolates', 'se'], n=int(c.at['S_isolates', 'n'])))
    R.append(row(grp, 'site jackknife, G=100 (0.5-deg centred cells)', 'siteK1_cluster_robust', k1.at[kset, 'beta'], k1.at[kset, 'se_site_jackknife'], n=int(k1.at[kset, 'n'])))
    R.append(row(grp, 'procedure P: lambda fixed, jackknife G=20, t19', 'tost_calibration_real', P.at[kset, 'beta'], P.at[kset, 'se'], df=19, n=int(P.at[kset, 'n'])))
    x = mc_int[(mc_int.model == mcm) & (mc_int.term == 'P0')].iloc[0]
    R.append(row(grp, 'MCMCglmm phylogeny + site (95% HPD)', 'exploratory_mcmc_site_interaction', x.post_mean, lo95=x.l95, hi95=x.u95, n=int(x.n)))
    term = 'P0 (MAGs)' if kset == 'MAGs' else 'P0 (isolates)'
    b, s, n = wsrow('b_speciesFE_clsite', term); R.append(row(grp.replace('species', 'genomes') + ', WITHIN species', 'species fixed effects, site-clustered SE', 'replication_within_species', b, s, n=n))
# interaction (difference between genome sources)
R.append(row('MAG x P0 interaction (difference in slopes)', 'phylolm lambda, model SE', 'exploratory_models', e.at['E_mag_x_P0', 'estimate'], e.at['E_mag_x_P0', 'se'], n=int(e.at['E_mag_x_P0', 'n'])))
x = mc_int[(mc_int.model == 'all_interaction') & (mc_int.term == 'P0:is_mag')].iloc[0]
R.append(row('MAG x P0 interaction (difference in slopes)', 'MCMCglmm phylogeny + site (95% HPD)', 'exploratory_mcmc_site_interaction', x.post_mean, lo95=x.l95, hi95=x.u95, n=int(x.n)))
b, s, n = wsrow('b_speciesFE_clsite', 'P0:is_mag'); R.append(row('MAG x P0 interaction, WITHIN species', 'species fixed effects, site-clustered SE', 'replication_within_species', b, s, n=n))
out = pd.DataFrame(R)
out.to_csv('results/group_equivalence.tsv', sep='\t', index=False, float_format='%.5g')
pd.set_option('display.width', 250)
print(out[['target', 'method', 'estimate', 'ci90_lo', 'ci90_hi', 'ci95_lo', 'ci95_hi', 'equiv_0.010', 'equiv_0.015', 'equiv_0.020']].round(4).to_string())
