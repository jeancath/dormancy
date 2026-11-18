#!/usr/bin/env python3
"""Apply the pre-registered decision rules to the confirmatory model table and write a Markdown summary.
usage: 17_decisions.py CONF_MODELS.tsv EXPL_MODELS.tsv H6_WALD.tsv CLUSTER_ROBUST.tsv OUT.md
Rules (prereg/PREREGISTRATION.md section 5), SESOI = 0.015 DORM_z units per SD of P0:
  supported            beta > 0, p < 0.05, and the 95% CI reaches or exceeds +SESOI
  practically null     95% CI entirely inside (-SESOI, +SESOI)
  contradicted         beta < 0 and p < 0.05
  inconclusive         anything else
"""
import sys, pandas as pd, numpy as np

SESOI = 0.015
conf, expl, wald, clus, out = sys.argv[1:6]
m = pd.read_csv(conf, sep='\t')
m['lo'] = m.estimate - 1.96 * m.se
m['hi'] = m.estimate + 1.96 * m.se


def rule(r):
    if pd.isna(r.estimate):
        return 'not estimable'
    if r.estimate < 0 and r.p < 0.05:
        return 'contradicted'
    if r.lo > -SESOI and r.hi < SESOI:
        return 'practically null' + (' (statistically non-zero)' if r.p < 0.05 else '')
    if r.estimate > 0 and r.p < 0.05 and r.hi >= SESOI:
        return 'supported'
    return 'inconclusive'


def fmt(df, decide=False):
    rows = ['| model | n | β | 95 % CI | p | p (BH) | λ/α |' + (' decision |' if decide else ''),
            '|---|---|---|---|---|---|---|' + ('---|' if decide else '')]
    for _, r in df.iterrows():
        if pd.isna(r.estimate):
            rows.append(f'| {r.id} | {r.n} | — | — | — | — | — |' + (' not estimable |' if decide else ''))
            continue
        pbh = '' if pd.isna(r.get('p_bh', np.nan)) else f'{r.p_bh:.3g}'
        rows.append(f'| {r.id} | {int(r.n):,} | {r.estimate:+.4f} | [{r.lo:+.4f}, {r.hi:+.4f}] | {r.p:.3g} | {pbh} | '
                    f'{r.lambda_or_alpha:.2f} |' + (f' {rule(r)} |' if decide else ''))
    return '\n'.join(rows)


L = ['# Confirmatory results: decisions under the pre-registered rules', '',
     f'The SESOI is ±{SESOI} DORM_z units per 1 SD of P0. A practically null result is one whose 95 % CI '
     'lies inside ±SESOI. The β values for count modules are on the log1p scale. Binary models are '
     'phyloglm logit slopes.', '']
groups = [('Primary (H1, direct and total effect) and the genome-size path', ['H1_DORM_P0', 'H1_DORM_P0_total', 'PATH_gsize_P0'], True),
          ('Modules (H1a–d, BH within the family)', ['H1a_sporulator', 'H1a_spo_strict', 'H1b_rpf', 'H1b_rpf_actino', 'H1c_hib', 'H1d_sr'], False),
          ('Predictable vs unpredictable variability (H6, BH within the family)', [i for i in m.id if i.startswith('H6')], False),
          ('Marine TA and streamlining (H3, H4)', ['H3_TA_marine', 'H3_TA_marine_mob', 'H3_TA_P0', 'H4_streamlined'], False),
          ('16S copy-number replication (H5)', ['H5_ssu', 'PATH_ssu_P0'], False),
          ('Pre-registered sensitivity analyses of H1', [i for i in m.id if i.startswith('S_')], True),
          ('Secondary predictors', [i for i in m.id if i.startswith('X_')], False)]
for title, ids, dec in groups:
    sub = m[m.id.isin(ids)].set_index('id').loc[[i for i in ids if i in set(m.id)]].reset_index()
    L += [f'## {title}', '', fmt(sub, dec), '']
w = pd.read_csv(wald, sep='\t')
L += ['## H6 Wald test (joint model, β residual − β seasonal)', '',
      w.round(4).to_markdown(index=False), '']
# H3 rule: supported if the marine coefficient is NOT significantly negative
h3 = m.set_index('id').loc['H3_TA_marine']
L += [f"**H3 decision:** the marine coefficient is {h3.estimate:+.3f} (p = {h3.p:.2g}). "
      + ('It is significantly negative, so **H3 is not supported**: after genome-size correction, marine genomes still carry fewer paired TA loci.'
         if (h3.estimate < 0 and h3.p < 0.05) else 'It is not significantly negative, so **H3 is supported**.'), '']
h5 = m.set_index('id').loc['H5_ssu']
L += [f"**H5 decision (replication):** β = {h5.estimate:+.3f}, p = {h5.p:.2g}. The replication "
      + ('succeeded.' if (h5.estimate > 0 and h5.p < 0.05) else 'failed (β > 0 and p < 0.05 were both required).'), '']
e = pd.read_csv(expl, sep='\t')
e['lo'] = e.estimate - 1.96 * e.se; e['hi'] = e.estimate + 1.96 * e.se
L += ['# Exploratory analyses (not pre-registered)', '', fmt(e), '']
c = pd.read_csv(clus, sep='\t')
L += ['## Site pseudoreplication (exploratory)', '', c.round(4).to_markdown(index=False), '']
open(out, 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
