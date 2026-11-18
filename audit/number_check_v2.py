#!/usr/bin/env python3
"""Re-check (orchestrator, 18 Nov 2026) of the numbers added or changed in the post-audit revision of
paper/academic/paper.md, paper/student/student_report.md and results/REPORT.md.
Each expected string is computed from the committed source files, then searched for in each document.
usage (repo root): python3 audit/number_check_v2.py -> audit/number_check_v2.tsv"""
import pandas as pd

import os
# Works in the private research repository and in the public release (paper/manuscript.md; no student version).
_paths = {'paper': 'paper/academic/paper.md' if os.path.exists('paper/academic/paper.md') else 'paper/manuscript.md',
          'student': 'paper/student/student_report.md', 'report': 'results/REPORT.md'}
docs = {k: open(v).read() for k, v in _paths.items() if os.path.exists(v)}
f = lambda x, d=3: (f'{x:+.{d}f}').replace('-', '−')
ci = lambda lo, hi, d=3: f'[{f(lo, d)}, {f(hi, d)}]'
rows = []
def chk(claim, expected, where):
    for w in [w for w in where if w in docs]:
        rows.append(dict(claim=claim, document=w, expected=expected, status='found' if expected in docs[w] else 'NOT FOUND'))

# within-species replication (HQ, species FE, site-clustered)
r = pd.read_csv('replication_within_species/results.tsv', sep='\t'); r = r[r['sample'] == 'hq']
def ws(model, resp, term):
    x = r[(r.model == model) & (r.response == resp) & (r.term == term)].iloc[0]
    return x.estimate, x.estimate - 1.96 * x.se, x.estimate + 1.96 * x.se
for (m, resp, term, lab) in [('b_speciesFE_clsite', 'DORM_z', 'P0:is_mag', 'within interaction'),
                             ('b_speciesFE_clsite', 'DORM_z', 'P0 (MAGs)', 'within MAG slope'),
                             ('b_speciesFE_clsite_comp95', 'DORM_z', 'P0 (MAGs)', 'within MAG slope >=95%'),
                             ('b_speciesFE_clsite_comp95', 'DORM_z', 'P0:is_mag', 'within interaction >=95%')]:
    e, lo, hi = ws(m, resp, term)
    chk(lab, f'{f(e)} {ci(lo, hi)}', ['paper', 'report'])
e, lo, hi = ws('a_speciesFE_clsite', 'DORM_z', 'P0'); chk('within pooled', f'{f(e)} {ci(lo, hi)}', ['paper', 'report'])
# MCMCglmm on GTDB branch lengths
m = pd.read_csv('results/confirmatory_mcmcglmm.tsv', sep='\t').set_index('response')
for k, lab in [('n_rpf', 'Rpf'), ('n_hib', 'hibernation'), ('n_sr', 'stringent response')]:
    chk(f'MCMC Poisson {lab}', f'{f(m.at[k, "post_mean"])} {ci(m.at[k, "l95"], m.at[k, "u95"])}', ['paper', 'report'])
s = pd.read_csv('results/siteMCMC_mcmcglmm_site.tsv', sep='\t').iloc[0]
chk('MCMC pooled phylo+site', f'{f(s.post_mean, 4)} {ci(s.l95, s.u95, 4)}', ['paper', 'report'])
i = pd.read_csv('results/exploratory_mcmc_site_interaction.tsv', sep='\t')
for model, term, lab in [('all_interaction', 'P0:is_mag', 'MCMC interaction'), ('MAGs', 'P0', 'MCMC MAGs'), ('isolates', 'P0', 'MCMC isolates')]:
    x = i[(i.model == model) & (i.term == term)]
    if len(x):
        x = x.iloc[0]; chk(lab, f'{f(x.post_mean)} {ci(x.l95, x.u95)}', ['paper', 'report'])
# audit-derived numbers
st = pd.read_csv('audit/strain_heterogeneity_tests.tsv', sep='\t').set_index('test')
chk('MAG slope, zero strain heterogeneity', f(st.at['DORM_z ~ P0, MAGs SH = 0', 'estimate']), ['paper'])
chk('MAG slope, SH > 0', f(st.at['DORM_z ~ P0, MAGs SH > 0', 'estimate']), ['paper'])
# corrected rounding / wording items
chk('MAG CI lower bound', '+0.018 [0.011, 0.025]', ['paper'])
chk('MAG CI (report)', '+0.018 [+0.011, +0.025]', ['report'])
chk('top-10 site share', '42 %', ['paper', 'report'])
chk('Balon <5 °C', '21 %', ['paper', 'student', 'report'])
chk('Balon >20 °C', '8 %', ['paper', 'student', 'report'])
chk('S_marine upper limit', '+0.0155', ['paper', 'report'])

# additions
eq = pd.read_csv('results/equivalence_tests.tsv', sep='\t').set_index('model')
chk('TOST H1 direct', f"| {eq.at['H1_DORM_P0', 'p_TOST']:.4f} |", ['paper'])
chk('TOST jackknife', f"| {eq.at['H1_direct with site-grouped jackknife SE (most conservative)', 'p_TOST']:.3f} |", ['paper'])
rv = pd.read_csv('results/review_models.tsv', sep='\t').set_index('id')
import math
for i in ['R_TA_comp95', 'R_TA_comp98', 'R_TA_comp98_cont1']:
    chk(f'TA {i}', f"{f(rv.at[i, 'estimate'])} | {round(100 * (1 - math.exp(rv.at[i, 'estimate'])))} %", ['paper'])
cd = pd.read_csv('results/collection_date_models.tsv', sep='\t').set_index('id')
e = cd.loc['D_H1_in_window']; chk('H1 in window', f"{f(e.estimate, 4)} [{f(e.estimate - 1.96 * e.se, 4)}, {f(e.estimate + 1.96 * e.se, 4)}]", ['paper', 'report'])

# additions: calibration of procedure P, procedure P on real data, collection-date tally, genome-source table
cal = pd.read_csv('results/tost_calibration_summary.tsv', sep='\t').set_index('true_beta')
pc = lambda b: f"{100 * cal.at[b, 'P_declare_equivalence']:.1f} % ({100 * cal.at[b, 'mc_lo95_clopper_pearson']:.1f}–{100 * cal.at[b, 'mc_hi95_clopper_pearson']:.1f} %"
for b in (-0.015, 0.015, 0.0, 0.005, 0.010):
    chk(f'calibration P(equiv) at {b}', pc(b), ['paper'])
sims = pd.read_csv('results/tost_calibration_sims.tsv', sep='\t', header=None, skiprows=1)
bd = sims[sims[0].abs() > 0.0149]; chk('calibration boundary combined', f"{100 * bd[7].mean():.2f} %", ['paper', 'report'])
P = pd.read_csv('results/tost_calibration_real.tsv', sep='\t').set_index('set')
dt = pd.read_csv('data_small/final/collection_dates.tsv', sep='\t')
an = pd.read_csv('data_small/final/analysis_table.tsv.gz', sep='\t', usecols=['accession', 'pilot'])
dt = dt[dt.accession.isin(an.loc[an.pilot == 0, 'accession'])]
chk('dated genomes', f"{dt.year_start.notna().sum():,}", ['paper', 'report'])
chk('in-window share of all', f"{100 * dt.in_window.astype(bool).mean():.1f} %", ['paper', 'report'])
chk('in-window share of dated', f"{100 * dt.in_window.astype(bool).sum() / dt.year_start.notna().sum():.1f} %", ['paper', 'report'])
g = pd.read_csv('results/group_equivalence.tsv', sep='\t')
x = g[(g.target == 'MAG x P0 interaction (difference in slopes)') & g.method.str.startswith('phylolm')].iloc[0]
chk('interaction phylolm', f"| {x.estimate:+.3f} | [{x.ci90_lo:.3f}, {x.ci90_hi:.3f}] | [{x.ci95_lo:.3f}, {x.ci95_hi:.3f}] |", ['paper'])

# additions: noise-scale power, one-sided MAG tests, composition-adjusted interaction, provenance counts
ic = pd.read_csv('results/interaction_composition.tsv', sep='\t').set_index('model')
for m in ic.index:
    chk(f'composition {m}', f"| +{ic.at[m, 'interaction']:.3f} |", ['paper', 'report'])
    chk(f'composition {m} P CI', f"[{ic.at[m, 'ci95_lo_jackknife']:.3f}, {ic.at[m, 'ci95_hi_jackknife']:.3f}]", ['paper', 'report'])
pv = pd.read_csv('data_small/final/biosample_retrieval_log.tsv', sep='\t')
dr = pd.read_csv('data_small/final/biosample_discarded_returns.tsv', sep='\t')
chk('provenance discarded', f"{len(dr)} misresolved", ['report'])
chk('provenance efetch', f"{int(pv.step1_own_record_returned.sum()):,} records fetched by accession", ['report'])

# additions: final procedure P-bar (results/genome_source_final_table.tsv), grouping sensitivity, scope, BioSample flow
ft = pd.read_csv('results/genome_source_final_table.tsv', sep='\t'); FN = ft[ft.role == 'FINAL'].set_index('quantity')
pb = FN.loc['pooled slope (H1 direct)']
chk('P-bar pooled 90% CI', f"{pb.ci90_lo:.4f}–{pb.ci90_hi:.4f}", ['paper'])
chk('P-bar pooled 90% CI (report)', f"[{pb.ci90_lo:.4f}, {pb.ci90_hi:.4f}]", ['report'])
chk('P-bar pooled p_TOST', f"p_TOST = {pb.p_tost:.3f}", ['paper', 'report'])
chk('P-bar pooled p vs 0', f"p = {pb.p_two_sided:.2f}", ['paper', 'report'])
mg = FN.loc['MAG-represented species: slope']
chk('P-bar MAG row', f"| {mg.estimate:+.3f} | [{mg.ci90_lo:.3f}, {mg.ci90_hi:.3f}] | [{mg.ci95_lo:.3f}, {mg.ci95_hi:.3f}] |", ['paper'])
chk('P-bar MAG p two-sided', f"p = {mg.p_two_sided:.3f}", ['paper', 'report'])
chk('P-bar MAG p > 0', f"p = {mg.p_one_sided_gt0:.3f}", ['paper', 'report'])
chk('P-bar MAG p > SESOI', f"p = {mg.p_one_sided_gt_sesoi:.2f}", ['paper', 'report'])
iso = FN.loc['isolate-represented species: slope']
chk('P-bar isolate p_TOST', f"p_TOST = {iso.p_tost:.3f}", ['paper', 'report'])
gs = pd.read_csv('results/grouping_sensitivity_runs.tsv', sep='\t'); rr = gs[(gs.set == 'H1_direct') & (gs.variant == 'random20')]
chk('groupings declaring equivalence', f"{int(rr.equiv.sum())} of {len(rr)}", ['paper'])
chk('rank of original grouping', f"{int((rr.se < rr[rr.rep == 1].se.iloc[0]).sum()) + 1}th-smallest SE of {len(rr)}", ['paper'])
oc = gs[(gs.set == 'H1_direct') & (gs.variant == 'one_cell')].iloc[0]
chk('delete-one-cell p_TOST', f"p_TOST = {oc.p_tost:.3f}", ['paper'])
sc = pd.read_csv('replication_within_species/inference_scope.tsv', sep='\t').set_index('quantity')['value'].astype(float)
chk('scope: species with variation', f"{int(sc['species with within-species P0 variation'])} of 1,473", ['report'])
bt = pd.read_csv('replication_within_species/inference_scope_by_type.tsv', sep='\t').set_index('species_type')
al = bt.loc['all']
chk('scope: coverage table, all row', f"| **All** | **{int(al.species):,}** | **{int(al.MAG_genomes):,} / {int(al.isolate_genomes):,}** | **{int(al.species_with_P0_variation)}** | **{int(al.species_with_P0_variation_among_MAGs)}** | **{int(al.species_with_P0_variation_among_isolates)}** | **{int(al.species_with_P0_variation_in_both_sources)}** |", ['paper'])
for k, lab in (('MAG-only', 'MAG-only'), ('isolate-only', 'Isolate-only'), ('mixed', 'Both sources')):
    r = bt.loc[k]
    chk(f'scope: coverage row {k}', f"| {lab} | {int(r.species):,} | {int(r.MAG_genomes):,} / {int(r.isolate_genomes):,} | {int(r.species_with_P0_variation)} | {int(r.species_with_P0_variation_among_MAGs)} | {int(r.species_with_P0_variation_among_isolates)} | {int(r.species_with_P0_variation_in_both_sources)} | {100 * r.share_of_within_species_P0_variation:.0f} % |", ['paper'])
# calibration by boundary and estimator family
cb = pd.read_csv('results/pbar_calibration_by_boundary.tsv', sep='\t')
for _, r in cb.iterrows():
    chk(f'calibration {r.dgm} {r.boundary} {r.K[:6]}', f"{100 * r.rate:.1f} ({100 * r.mc_lo95:.1f}–{100 * r.mc_hi95:.1f})", ['paper'])
ef = pd.read_csv('results/pbar_estimator_family.tsv', sep='\t')
chk('largest sensitivity p_TOST', f"largest p_TOST among them is {ef.iloc[-1].p_tost:.3f}", ['paper'])
sd = pd.read_csv('results/pbar_se_distribution.tsv', sep='\t').set_index('set').loc['H1_direct']
chk('SE range over groupings', f"{sd.se_min:.4f}–{sd.se_max:.4f}", ['paper'])
chk('scope: MAG-identifying species (report)', f"{int(sc['species with MAG genomes at >= 2 P0 values'])} species (1,214 MAGs)", ['report'])
fl = pd.read_csv('results/biosample_retrieval_flow.tsv', sep='\t').set_index('step')['count']
chk('BioSample returned', f"| {fl['step 1: records returned']:,} |", ['paper', 'report'])
chk('BioSample nothing returned', f"| {fl['  of which nothing was returned for them']} |", ['paper'])
for w, txt in docs.items():
    for bad in ['held-out pilot (−0.25', 'replicated in the held-out pilot', 'about 17 % below 5 °C', 'about 5 % above',
                'metagenomic observation of Lennon & Jones', '94 % of the pattern', 'before any confirmatory genome was annotated',
                'Nine of ten', '{{MC_', 'reversal vanishes', 'disappeared. The interaction', 'not an effect of genome source', 'The difference vanished', 'The simulated data were noisier', 'power is 99 % at ±0.005', 'p_TOST = 0.0005: the']:
        rows.append(dict(claim=f'withdrawn wording absent: "{bad}"', document=w, expected='absent',
                         status='absent' if bad not in txt else 'STILL PRESENT'))
out = pd.DataFrame(rows); out.to_csv('audit/number_check_v2.tsv', sep='\t', index=False)
print(out.to_string()); print('\nproblems:', ((out.status == 'NOT FOUND') | (out.status == 'STILL PRESENT')).sum())
