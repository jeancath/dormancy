#!/usr/bin/env python3
"""Added 18 Nov 2026: ONE table for the pooled, MAG, isolate and interaction slopes, linking every number to its
procedure and source file, and a reconciliation of numbers reported in earlier rounds that are now superseded.

Final procedure ("P-bar"): procedure P (lambda fixed at the full-data estimate, GLS point estimate, delete-a-group
jackknife over 20 random groups of 0.5-degree cells, t on 19 df) with the jackknife VARIANCE AVERAGED over 200 random
groupings (SE = sqrt(mean SE_k^2)). Reason: on the real data the SE from one grouping varies widely between groupings
(scripts/30_grouping_sensitivity.R); the originally reported grouping gave the 4th-smallest SE of 200 for the pooled
slope. Calibration of the averaged-SE version: scripts/31 and 31b. For the interaction, the average is over 50 groupings
(scripts/28_interaction_composition.R).
usage (repo root): python3 scripts/32_final_genome_source_table.py
  -> results/genome_source_final_table.tsv, results/superseded_numbers.tsv"""
import numpy as np, pandas as pd
from scipy.stats import t as T, norm
M = 0.015
runs = pd.read_csv('results/grouping_sensitivity_runs.tsv', sep='\t')
g = pd.read_csv('results/group_equivalence.tsv', sep='\t')
ic = pd.read_csv('results/interaction_composition.tsv', sep='\t').set_index('model')
rows = []
def add(quantity, procedure, role, b, se, df, source, lo95=None, hi95=None, note=''):
    r = dict(quantity=quantity, procedure=procedure, role=role, estimate=b, se=se, df=df)
    if se is not None and not np.isnan(se):
        d = (lambda x: T.sf(x, df)) if df else norm.sf
        q90 = T.ppf(.95, df) if df else norm.ppf(.95); q95 = T.ppf(.975, df) if df else norm.ppf(.975)
        r.update(ci90_lo=b - q90 * se, ci90_hi=b + q90 * se, ci95_lo=b - q95 * se, ci95_hi=b + q95 * se,
                 p_two_sided=2 * d(abs(b) / se), p_one_sided_gt0=d(b / se), p_one_sided_lt0=d(-b / se),
                 p_one_sided_gt_sesoi=d((b - M) / se), p_one_sided_lt_minus_sesoi=d((-M - b) / se),
                 p_tost=max(d((b + M) / se), d((M - b) / se)), equivalent_at_0015=bool(b - q90 * se > -M and b + q90 * se < M))
    else:
        r.update(ci95_lo=lo95, ci95_hi=hi95, equivalent_at_0015=bool(lo95 > -M and hi95 < M))
    r.update(source=source, note=note); rows.append(r)
lab = {'H1_direct': 'pooled slope (H1 direct)', 'MAGs': 'MAG-represented species: slope', 'isolates': 'isolate-represented species: slope'}
for s in ('H1_direct', 'MAGs', 'isolates'):
    x = runs[(runs.set == s) & (runs.variant == 'random20')]; b = x.beta.iloc[0]
    add(lab[s], 'P-bar: P with jackknife variance averaged over 200 random groupings, t19', 'FINAL', b,
        float(np.sqrt((x.se ** 2).mean())), 19, 'results/grouping_sensitivity_runs.tsv')
    o = runs[(runs.set == s) & (runs.variant == 'one_cell')].iloc[0]
    add(lab[s], 'delete-one-cell jackknife (no grouping randomness), t on cells-1 df', 'sensitivity', b, o.se, int(o.G) - 1,
        'results/grouping_sensitivity_runs.tsv')
    c = runs[(runs.set == s) & (runs.variant == 'contiguous20')]
    add(lab[s], '20 spatially contiguous groups (k-means), variance averaged over 20 k-means seeds, t19', 'sensitivity', b,
        float(np.sqrt((c.se ** 2).mean())), 19, 'results/grouping_sensitivity_runs.tsv')
    first = x[x.rep == 1].iloc[0]
    add(lab[s], 'P with the single original grouping (seed 20261118), t19', 'superseded (earlier report)', b, first.se, 19,
        'results/tost_calibration_real.tsv', note=f'SE rank {int((x.se < first.se).sum()) + 1} of 200 groupings')
    if s == 'H1_direct':
        lf = runs[(runs.set == s) & (runs.variant == 'lambda_free')].iloc[0]
        add(lab[s], 'P with lambda re-estimated in every jackknife replicate (original grouping), t19', 'sensitivity', lf.beta, lf.se, 19,
            'results/grouping_sensitivity_runs.tsv')
# model-based and other estimators from results/group_equivalence.tsv
tg = {'pooled (H1 direct)': lab['H1_direct'], 'MAG-represented species': lab['MAGs'], 'isolate-represented species': lab['isolates'],
      'MAG x P0 interaction (difference in slopes)': 'MAG x P0 interaction (between species)'}
for _, r in g.iterrows():
    if r.target not in tg or r.method.startswith(('PRIMARY', 'procedure P')): continue
    if r.method.startswith('MCMCglmm'):
        add(tg[r.target], r.method, 'sensitivity', r.estimate, None, None, 'results/' + r.source + '.tsv', r.ci95_lo, r.ci95_hi)
    else:
        add(tg[r.target], r.method, 'sensitivity', r.estimate, r.se, None, 'results/' + r.source + '.tsv')
m0 = ic.loc['M0_unadjusted']
add('MAG x P0 interaction (between species)', f'P-bar: jackknife variance averaged over {int(m0.n_groupings)} random groupings, t19', 'FINAL',
    m0.interaction, m0.se_jackknife, 19, 'results/interaction_composition.tsv')
m5 = ic.loc['M5_common_support']
add('MAG x P0 interaction, common support (18 phylum x habitat strata)', f'P-bar ({int(m5.n_groupings)} groupings), t19', 'sensitivity',
    m5.interaction, m5.se_jackknife, 19, 'results/interaction_composition.tsv')
for _, r in g[g.target.str.contains('WITHIN')].iterrows():
    add(r.target, r.method, 'FINAL (within-species analysis)', r.estimate, r.se, None, 'replication_within_species/results.tsv')
out = pd.DataFrame(rows)
out.to_csv('results/genome_source_final_table.tsv', sep='\t', index=False, float_format='%.8g')
pd.set_option('display.width', 250)
print(out[['quantity', 'procedure', 'role', 'estimate', 'se', 'ci90_lo', 'ci90_hi', 'ci95_lo', 'ci95_hi', 'p_two_sided',
           'p_one_sided_gt0', 'p_one_sided_gt_sesoi', 'p_tost']].round(4).to_string())

# ---- reconciliation of numbers reported in earlier rounds (old values as reported; new values from the table above)
def fin(q, proc_start='P-bar'):
    return out[(out.quantity == q) & out.procedure.str.startswith(proc_start)].iloc[0]
P_, Mg, Is, In = fin(lab['H1_direct']), fin(lab['MAGs']), fin(lab['isolates']), fin('MAG x P0 interaction (between species)')
f4 = lambda x: f'{x:.4f}'; f3 = lambda x: f'{x:.3f}'
S = [
 ('pooled slope: SE, 90 % CI, p_TOST', 'SE 0.0022; [0.0026, 0.0103]; p_TOST = 0.0005', 'earlier drafts of the paper and report',
  f'SE {f4(P_.se)}; [{f4(P_.ci90_lo)}, {f4(P_.ci90_hi)}]; p_TOST = {f3(P_.p_tost)}',
  'the single original random grouping gave the 4th-smallest SE of 200 (scripts/30); now averaged over 200 groupings'),
 ('MAG slope in the abstract', '+0.018, 90 % CI 0.008-0.028, 95 % CI 0.006-0.029', 'an earlier abstract',
  f'+{f3(Mg.estimate)}, 90 % CI {f3(Mg.ci90_lo)}-{f3(Mg.ci90_hi)}, 95 % CI {f3(Mg.ci95_lo)}-{f3(Mg.ci95_hi)}',
  'the abstract mixed three estimators (phylolm estimate, G=100 jackknife 90 % CI, single-grouping P 95 % CI); now one procedure'),
 ('MAG slope tests', 'p = 0.005 two-sided; p = 0.003 for beta > 0; p = 0.33 for beta > 0.015', 'an earlier draft',
  f'p = {f3(Mg.p_two_sided)} two-sided; p = {f3(Mg.p_one_sided_gt0)} for beta > 0; p = {Mg.p_one_sided_gt_sesoi:.2f} for beta > 0.015',
  '0.005 and 0.003 were the two-sided and one-sided p of the same single-grouping test; replaced by the averaged-SE test'),
 ('isolate slope', '90 % CI [-0.016, -0.001]; p_TOST = 0.090; p = 0.055 two-sided', 'an earlier draft',
  f'90 % CI [{f3(Is.ci90_lo)}, {f3(Is.ci90_hi)}]; p_TOST = {f3(Is.p_tost)}; p = {f3(Is.p_two_sided)} two-sided', 'single grouping replaced by the averaged SE'),
 ('MAG x P0 interaction, cluster-robust', 'single-grouping jackknife CIs in the composition table', 'an earlier draft',
  f'+{f3(In.estimate)}, 95 % CI [{f3(In.ci95_lo)}, {f3(In.ci95_hi)}] (50 groupings)', 'single grouping replaced by the averaged SE'),
 ('power at +-0.010 at the "real-data noise level"', '69-70 %', 'an earlier draft and deviations item 21',
  'withdrawn; 37 % (simulations, whose slope SD 0.0035 equals the typical real-data SE)',
  'the premise (real-data SE 0.0022) came from the unusual grouping; the simulations were not noisier than the data'),
]
pd.DataFrame(S, columns=['quantity', 'old_value', 'where_reported', 'final_value', 'reason']).to_csv(
    'results/superseded_numbers.tsv', sep='\t', index=False)
print(pd.DataFrame(S).iloc[:, [0, 1, 3]].to_string())
