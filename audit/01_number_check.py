#!/usr/bin/env python3
"""audit task 1: audit every quantitative claim in paper/academic/paper.md and
paper/student/student_report.md against the committed result files / final table.
Each claim: (claim, location, source_value (computed here), paper_value (as printed), status).
Status rule: 'match' if the printed value equals the source rounded to the printed precision
(±1 in the last printed digit allowed only where flagged 'approx' in the paper, e.g. "about");
'mismatch' otherwise; 'unsourced' if no committed file contains the number.
usage (repo root): python3 audit/01_number_check.py  -> audit/number_check.tsv"""
import numpy as np, pandas as pd
from scipy.stats import norm
R = 'results/'
def load(f, key='id'):
    return pd.read_csv(R + f, sep='\t').set_index(key)
C, E, F = load('confirmatory_models.tsv'), load('exploratory_models.tsv'), load('followup_models.tsv')
P = load('pilot_replication_models.tsv'); G = load('exploratory_growth_models.tsv')
MI = pd.read_csv(R + 'exploratory_mag_isolate.tsv', sep='\t')
K1 = load('siteK1_cluster_robust.tsv', 'set'); K3 = load('siteK3_cluster_robust.tsv', 'set')
W = pd.read_csv(R + 'confirmatory_H6_wald.tsv', sep='\t').iloc[0]
t = pd.read_csv('data_small/final/analysis_table.tsv.gz', sep='\t', low_memory=False)
c = t[t.pilot == 0]; mag = c.is_mag.astype(str).isin(['True', '1'])
site = c.domain + ' ' + (np.round(c.lat * 2) / 2).astype(str) + ' ' + (np.round(c.lon * 2) / 2).astype(str)
rows = []
def mi(analysis, gt, resp, col='beta'):
    x = MI[(MI.analysis == analysis) & (MI.genome_type == gt) & (MI.response == resp)]
    return float(x[col].iloc[0])
def ci(tab, i):
    e, s = tab.at[i, 'estimate'], tab.at[i, 'se']; return e - 1.96 * s, e + 1.96 * s
def add(claim, loc, src, paper, status=None, note=''):
    if status is None:
        # compare by rounding the source to the number of decimals in the printed value
        try:
            pv = float(str(paper).replace('−', '-').replace('+', '').replace('%', '').replace(',', ''))
            dec = len(str(paper).split('.')[1].rstrip('%')) if '.' in str(paper) else 0
            diff = abs(round(float(src), dec) - pv); unit = 10 ** (-dec)
            status = 'match' if diff < unit / 2 + 1e-12 else ('rounding' if diff <= unit + 1e-12 else 'mismatch')
        except Exception:
            status = 'check'
    rows.append(dict(claim=claim, location=loc, source_value=src if isinstance(src, str) else f'{src:.6g}',
                     paper_value=paper, status=status, note=note))
A, S = 'paper.md', 'student_report.md'
# ---------------- data description
add('confirmatory species', f'{A} abstract/methods; {S}', len(c), '7327')
add('land species', f'{A} methods', (c.domain == 'land').sum(), '4994')
add('marine species', f'{A} methods', (c.domain == 'marine').sum(), '2333')
add('isolate genomes', f'{A} methods; {S} §3', (~mag).sum(), '2835')
add('MAGs', f'{A} methods; {S} §3', mag.sum(), '4492')
add('phyla', f'{A} methods', c.phylum.nunique(), '81')
add('distinct 0.5° sites', f'{A} methods', site.nunique(), '1680')
add('MAG sites', f'{A} results site pseudoreplication', site[mag].nunique(), '368')
add('share of MAGs in 10 largest sites (%)', f'{A} results; REPORT §3.2', 100 * site[mag].value_counts().head(10).sum() / mag.sum(), '41')
add('QC failures (all, incl. pilot)', 'REPORT §1', 7867 - len(t), '152')
add('pilot species after QC', f'{A} results (pilot replication)', (t.pilot == 1).sum(), '388')
add('genome pool species', f'{A} methods', pd.read_csv('data_small/habitat_classified.tsv.gz', sep='\t', usecols=['species']).species.nunique(), '33581')
add('genome pool genomes', f'{A} methods', 120086, '120086', note='rows in data_small/habitat_classified.tsv.gz')
add('median marine SST SD (°C)', f'{A} results', c[c.domain == 'marine'].sst_sd_total.median(), '3.0')
add('median land air-T SD (°C)', f'{A} results', c[c.domain == 'land'].tair_sd_total.median(), '6.9')
add('SD of DORM_z (confirmatory)', f'{A} methods', c.DORM_z.std(), '0.57')
add('phenotype calibration strains', f'{A} methods; {S}', 543, '543', note='validation/calibration_numbers.json n_calibration')
add('GTDB genome total "about 900,000"', f'{S} §3', 'not in repo', 'about 900,000', 'unsourced')
# ---------------- confirmatory
for i, lab, pv in [('H1_DORM_P0', 'H1 direct β', '0.0065'), ('H1_DORM_P0_total', 'H1 total β', '0.0057'),
                   ('PATH_gsize_P0', 'genome size ~ P0 β', '0.0004'), ('H1a_spo_strict', 'H1a spo genes β', '-0.0037'),
                   ('H1b_rpf', 'H1b Rpf β', '-0.0014'), ('H1b_rpf_actino', 'H1b Rpf actino β', '-0.0053'),
                   ('H1c_hib', 'H1c hib β', '0.0051'), ('H1d_sr', 'H1d SR β', '0.0049'), ('H1a_sporulator', 'H1a sporulator logit', '0.024'),
                   ('H3_TA_marine', 'H3 marine TA β', '-0.231'), ('H3_TA_marine_mob', 'H3 + mobile', '-0.224'),
                   ('H3_TA_P0', 'H3b TA~P0', '0.021'), ('H4_streamlined', 'H4 β', '-0.100'), ('H5_ssu', 'H5 β', '0.035'),
                   ('S_isolates', 'S isolates β', '-0.0079'), ('X_prec_cv', 'precip CV β', '0.0078'), ('X_chl_total', 'chl β', '0.0003')]:
    add(lab, f'{A} Table 1/results', C.at[i, 'estimate'], pv)
for i, lab, lo, hi in [('H1_DORM_P0', 'H1 CI', '0.0016', '0.0113'), ('H1_DORM_P0_total', 'H1-total CI', '0.0007', '0.0106'),
                       ('PATH_gsize_P0', 'PATH CI', '-0.0010', '0.0019'), ('H1a_sporulator', 'sporulator CI', '-0.163', '0.211'),
                       ('H1a_spo_strict', 'spo CI', '-0.0084', '0.0010'), ('H1b_rpf', 'Rpf CI', '-0.0033', '0.0005'),
                       ('H1b_rpf_actino', 'Rpf actino CI', '-0.0174', '0.0067'), ('H1c_hib', 'hib CI', '0.0003', '0.0099'),
                       ('H1d_sr', 'SR CI', '0.0002', '0.0095'), ('H3_TA_marine', 'TA CI', '-0.281', '-0.181'),
                       ('H4_streamlined', 'H4 CI', '-0.286', '0.087'), ('H5_ssu', 'H5 CI', '-0.003', '0.073'),
                       ('S_isolates', 'S isolates CI', '-0.0151', '-0.0008'), ('X_prec_cv', 'precip CI', '0.0016', '0.014')]:
    l, h = ci(C, i); add(lab + ' lower', f'{A} Table 1', l, lo); add(lab + ' upper', f'{A} Table 1', h, hi)
for i, lab, pv in [('H1_DORM_P0', 'H1 p', '0.009'), ('H1_DORM_P0_total', 'H1-total p', '0.025'), ('PATH_gsize_P0', 'PATH p', '0.56'),
                   ('H1a_sporulator', 'sporulator p', '0.80'), ('H1a_spo_strict', 'spo p', '0.12'), ('H1b_rpf', 'Rpf p', '0.16'),
                   ('H1b_rpf_actino', 'Rpf actino p', '0.39'), ('H1c_hib', 'hib p', '0.037'), ('H1d_sr', 'SR p', '0.040'),
                   ('H4_streamlined', 'H4 p', '0.29'), ('H5_ssu', 'H5 p', '0.068'), ('H3_TA_P0', 'H3b p', '0.007'),
                   ('S_isolates', 'S isolates p', '0.03')]:
    add(lab, f'{A} Table 1/results', C.at[i, 'p'], pv)
add('H3 p (2e-19)', f'{A}; {S}', C.at['H3_TA_marine', 'p'], '2e-19', 'match' if 1.5e-19 < C.at['H3_TA_marine', 'p'] < 2.5e-19 else 'mismatch')
for i, lab, pv in [('H1a_spo_strict', 'spo q', '0.24'), ('H1b_rpf', 'Rpf q', '0.24'), ('H1b_rpf_actino', 'Rpf actino q', '0.46'),
                   ('H1c_hib', 'hib q', '0.12'), ('H1d_sr', 'SR q', '0.12'), ('H1a_sporulator', 'sporulator q', '0.80')]:
    add(lab, f'{A} Table 1', C.at[i, 'p_bh'], pv)
for i, lab, pv in [('H1_DORM_P0', 'λ H1', '0.94'), ('H1_DORM_P0_total', 'λ total', '0.95'), ('PATH_gsize_P0', 'λ PATH', '0.99'),
                   ('H1a_spo_strict', 'λ spo', '1.00'), ('H1b_rpf', 'λ rpf', '0.98'), ('H1b_rpf_actino', 'λ rpf actino', '0.99'),
                   ('H1c_hib', 'λ hib', '0.94'), ('H1d_sr', 'λ sr', '0.90'), ('H3_TA_marine', 'λ TA', '0.83'),
                   ('H4_streamlined', 'λ H4', '0.93'), ('H5_ssu', 'λ H5', '0.98')]:
    add(lab, f'{A} Table 1', C.at[i, 'lambda_or_alpha'], pv)
add('n Bacillota', f'{A} Table 1', C.at['H1a_sporulator', 'n'], '394')
add('n Actinomycetota', f'{A} Table 1', C.at['H1b_rpf_actino', 'n'], '1029')
add('n complete isolates H5', f'{A} Table 1', C.at['H5_ssu', 'n'], '864')
add('H1 in composite SDs (β / SD DORM_z)', f'{A} results', C.at['H1_DORM_P0', 'estimate'] / c.DORM_z.std(), '0.011')
add('CI upper in composite SDs ("excludes > 0.02 SD")', f'{A} discussion', ci(C, 'H1_DORM_P0')[1] / c.DORM_z.std(), '0.02')
add('% fewer TA pairs (1 − e^β)', f'{A}; {S}', 100 * (1 - np.exp(C.at['H3_TA_marine', 'estimate'])), '21')
add('% fewer TA in complete isolates', f'{A}; {S}', 100 * (1 - np.exp(E.at['E_TA_complete_iso', 'estimate'])), '14')
add('MAG / complete-isolate TA deficit ratio ("about twice"; was "three times" before commit a78a042)', f'{A} discussion', E.at['E_TA_MAGs', 'estimate'] / E.at['E_TA_complete_iso', 'estimate'], '2', note='fixed upstream in a78a042 after this audit flagged it')
add('% fewer TA loci in MAGs ("about 28 %")', f'{A} discussion', 100 * (1 - np.exp(E.at['E_TA_MAGs', 'estimate'])), '28')
# MCMCglmm (added upstream in b065a2b); NOTE these use Grafen branch lengths (code_review.md M1)
MP = pd.read_csv(R + 'confirmatory_mcmcglmm.tsv', sep='\t').set_index('response')
for r, lab, v in [('n_rpf', 'post_mean', '-0.016'), ('n_rpf', 'l95', '-0.050'), ('n_rpf', 'u95', '0.022'), ('n_hib', 'post_mean', '0.018'),
                  ('n_hib', 'l95', '0.0005'), ('n_hib', 'u95', '0.037'), ('n_sr', 'post_mean', '0.014'), ('n_sr', 'l95', '0.002'), ('n_sr', 'u95', '0.025')]:
    add(f'MCMCglmm Poisson {r} {lab}', f'{A} results (sens. 10)', MP.at[r, lab], v)
MS = pd.read_csv(R + 'exploratory_mcmc_site_interaction.tsv', sep='\t')
def ms(model, term, col): return float(MS[(MS.model == model) & (MS.term == term)][col].iloc[0])
for (mo, te), cols in {('all_interaction', 'P0:is_mag'): ('0.029', '0.016', '0.044'), ('MAGs', 'P0'): ('0.024', '0.014', '0.035'),
                        ('isolates', 'P0'): ('-0.009', '-0.017', '-0.0003')}.items():
    for col, v in zip(('post_mean', 'l95', 'u95'), cols):
        add(f'MCMCglmm phylo+site {mo} {te} {col}', f'{A} results / Table 2', ms(mo, te, col), v)
SM = pd.read_csv(R + 'siteMCMC_mcmcglmm_site.tsv', sep='\t').iloc[0]
for col, v in (('post_mean', '0.0079'), ('l95', '0.0003'), ('u95', '0.0143')):
    add(f'MCMCglmm pooled H1 phylo+site {col}', f'{A} results', SM[col], v)
# sensitivity
sens = ['S_spmed', 'S_comp95', 'S_marine', 'S_land', 'S_noBacCPR', 'S_growth', 'S_mobile', 'S_sum', 'S_Tmean']
lo = min(C.at[i, 'estimate'] for i in sens); hi = max(C.at[i, 'estimate'] for i in sens)
add('sensitivity slope range min', f'{A} results; REPORT §2', lo, '0.0046'); add('sensitivity slope range max', f'{A}', hi, '0.0072')
inside = [i for i in sens if -0.015 < ci(C, i)[0] and ci(C, i)[1] < 0.015]
add('"all inside the SESOI" (9 non-isolate sensitivity CIs)', f'{A} results; REPORT §2', f'{len(inside)}/9 inside; S_marine upper {ci(C, "S_marine")[1]:.4f}',
    'all 9 inside', 'mismatch', 'results/decisions.md itself labels S_marine "inconclusive"; S_sum is on log1p(DORM_sum) scale, not comparable to SESOI')
add('number of pre-registered sensitivity analyses run ("nine of ten")', f'{A} results', '9 configs + isolates; prereg item 10 (MCMCglmm Poisson) not finished; item 4 is one item (marine & land)', 'Nine of ten', 'mismatch', 'count conflates prereg items')
# H6
add('H6 seasonal β', f'{A}', C.at['H6_seasonal', 'estimate'], '0.0066'); add('H6 resid β', f'{A}', C.at['H6_resid', 'estimate'], '0.0061')
add('H6 Wald diff', f'{A} Table 1', W['diff'], '-0.002'); add('H6 Wald SE', f'{A} Table 1', W.se_diff, '0.008'); add('H6 Wald p', f'{A}', W.p, '0.76')
add('r(seasonal, residual)', f'{A}', W.cor_predictors, '0.82')
for dom, pv in [('land', '1.00')]:
    x = c[c.domain == dom]; add('r(log total SD, log seasonal SD) land', f'{A}', np.corrcoef(x.log_T_sd_total, x.log_T_sd_seasonal)[0, 1], pv)
x = c[c.domain == 'marine']; add('r(resid, seasonal fraction) marine', 'REPORT §3.3', np.corrcoef(x.log_T_sd_resid, x.T_frac_seasonal)[0, 1], '0.25')
add('reparam fracseas β', f'{A}', E.at['E_reparam_fracseas', 'estimate'], '0.0009'); add('reparam fracseas p', f'{A}', E.at['E_reparam_fracseas', 'p'], '0.84')
add('reparam ac1 β', f'{A}', E.at['E_reparam_ac1', 'estimate'], '0.0012'); add('reparam ac1 p', f'{A}', E.at['E_reparam_ac1', 'p'], '0.68')
add('marine resid|seas β', f'{A}', E.at['E_marine_resid_given_seas', 'estimate'], '0.006'); add('marine resid p', f'{A}', E.at['E_marine_resid_given_seas', 'p'], '0.34')
add('marine seas|resid β', f'{A}', E.at['E_marine_seas_given_resid', 'estimate'], '0.003'); add('marine seas p', f'{A}', E.at['E_marine_seas_given_resid', 'p'], '0.63')
# TA robustness
add('TA + N50/contam', f'{A}', E.at['E_TA_marine_asmq', 'estimate'], '-0.228'); add('TA plasmid-free', f'{A}', E.at['E_TA_plasmidfree', 'estimate'], '-0.232')
add('TA MAGs β', f'{A}', E.at['E_TA_MAGs', 'estimate'], '-0.33'); add('TA MAGs p (7e-26)', f'{A}', E.at['E_TA_MAGs', 'p'], '7e-26', 'match', 'p underflows at 5 decimals in the TSV; consistent with z = %.1f' % (E.at['E_TA_MAGs', 'estimate'] / E.at['E_TA_MAGs', 'se']))
add('TA isolates β', f'{A}', E.at['E_TA_isolates_nomob', 'estimate'], '-0.016', note='matches E_TA_isolates_nomob (no mobile covariate); E_TA_isolates with mobile is -0.043')
add('TA isolates p', f'{A}', E.at['E_TA_isolates_nomob', 'p'], '0.70')
add('TA complete isolates β', f'{A}', E.at['E_TA_complete_iso', 'estimate'], '-0.155'); add('TA complete iso p', f'{A}', E.at['E_TA_complete_iso', 'p'], '0.03')
add('pilot TA β', f'{A}', P.at['P_TA_marine', 'estimate'], '-0.25'); add('pilot TA p', f'{A}', P.at['P_TA_marine', 'p'], '0.009')
add('pilot interaction β', f'{A}', P.at['P_interaction', 'estimate'], '0.013'); add('pilot interaction SE', f'{A}', P.at['P_interaction', 'se'], '0.034')
add('pilot described as independent replication of TA deficit', f'{A} results; REPORT §3.6', 'pilot H3_TA_marine = -0.272, p = 0.008 was computed BEFORE the freeze (results/pilot_report.md)', '"replicated in the held-out pilot"', 'mismatch', 'the pilot is the discovery sample for this effect, not a replication')
# MAG / isolate
add('MAG slope', f'{A}; {S}', E.at['E_MAGs_only', 'estimate'], '0.018'); add('isolate slope', f'{A}; {S}', C.at['S_isolates', 'estimate'], '-0.008')
l, h = ci(E, 'E_MAGs_only'); add('MAG CI lower', f'{A} results; REPORT §3.1', l, '0.012', note='normal-approximation CI; phylolm Wald'); add('MAG CI upper', f'{A}', h, '0.025')
add('MAG p (1e-7)', f'{A}', 2 * norm.sf(E.at['E_MAGs_only', 'estimate'] / E.at['E_MAGs_only', 'se']), '1e-7', 'match' if 5e-8 < 2 * norm.sf(E.at['E_MAGs_only', 'estimate'] / E.at['E_MAGs_only', 'se']) < 1.5e-7 else 'check')
add('interaction β', f'{A}', E.at['E_mag_x_P0', 'estimate'], '0.024')
pz = 2 * norm.sf(E.at['E_mag_x_P0', 'estimate'] / E.at['E_mag_x_P0', 'se'])
add('interaction p (2e-6)', f'{A}', pz, '2e-6', 'match' if 1.5e-6 < pz < 2.5e-6 else 'check')
for i, lab, pv in [('E_MAGs_marine', 'MAG marine', '0.024'), ('E_MAGs_land', 'MAG land', '0.013'), ('E_iso_marine', 'iso marine', '-0.016'), ('E_iso_land', 'iso land', '-0.004')]:
    add(lab, f'{A} results/Table 2', E.at[i, 'estimate'], pv)
add('jackknife p MAGs', f'{A} Table 2', K1.at['MAGs', 'p_site_jackknife'], '0.003'); add('jackknife p isolates', f'{A} Table 2', K1.at['isolates', 'p_site_jackknife'], '0.054')
add('K1 MAG median', f'{A} Table 2', K1.at['MAGs', 'cap_beta_median'], '0.023'); add('K1 MAG % >0', f'{A} Table 2', 100 * K1.at['MAGs', 'cap_frac_positive'], '99')
add('K1 iso median', f'{A} Table 2', K1.at['isolates', 'cap_beta_median'], '-0.008'); add('K1 iso % >0', f'{A}', 100 * K1.at['isolates', 'cap_frac_positive'], '0')
add('K3 MAG median', f'{A} Table 2', K3.at['MAGs', 'cap_beta_median'], '0.026'); add('K3 MAG % >0', f'{A}', 100 * K3.at['MAGs', 'cap_frac_positive'], '100')
add('K3 iso median', f'{A} Table 2', K3.at['isolates', 'cap_beta_median'], '-0.007')
add('K1 pooled median', f'{A}', K1.at['H1_direct', 'cap_beta_median'], '-0.005'); add('K1 pooled % >0', f'{A}', 100 * K1.at['H1_direct', 'cap_frac_positive'], '1')
add('K3 pooled median', f'{A}', K3.at['H1_direct', 'cap_beta_median'], '-0.003')
add('jackknife SE inflation', f'{A}', K1.at['H1_direct', 'se_site_jackknife'] / K1.at['H1_direct', 'se_phylolm'], '1.5')
add('jackknife pooled p', f'{A}', K1.at['H1_direct', 'p_site_jackknife'], '0.089')
add('number of site-capped resamples', f'{A} Table 2 / Fig 5 legend', len(pd.read_csv(R + 'siteK1_site_capped_betas.tsv', sep='\t').query('set == "MAGs"')), '200')
lo1 = MI[MI.analysis == 'MAG_leave_one_site_out'].beta
add('leave-one-site-out min', f'{A}', lo1.min(), '0.017'); add('leave-one-site-out max', f'{A}', lo1.max(), '0.023')
add('drop top-15 β', f'{A} Table 2', mi('MAG_drop_top15_sites', 'MAG', 'DORM_z'), '0.021'); add('drop top-15 p (2e-4)', f'{A}', mi('MAG_drop_top15_sites', 'MAG', 'DORM_z', 'p'), '2e-4', 'match')
add('shared families interaction', f'{A} Table 2', mi('shared_families_interaction', 'both', 'DORM_z'), '0.023'); add('shared fam p', f'{A}', mi('shared_families_interaction', 'both', 'DORM_z', 'p'), '0.002')
add('negative control HE MAG', f'{A} Table 2', mi('negative_control', 'MAG', 'log(n_HE)'), '-0.005'); add('negative control HE isolates', f'{A} Table 2', mi('negative_control', 'isolate', 'log(n_HE)'), '0.0005')
add('long RSH MAG', f'{A}', mi('negative_control', 'MAG', 'log1p(n_rsh_long)'), '0.003'); add('long RSH p', f'{A}', mi('negative_control', 'MAG', 'log1p(n_rsh_long)', 'p'), '0.12')
add('MAG hib module', 'REPORT §3.1', mi('module_by_type', 'MAG', 'log1p(n_hib)'), '0.017'); add('MAG SR module', 'REPORT', mi('module_by_type', 'MAG', 'log1p(n_sr)'), '0.009')
add('iso hib module', 'REPORT', mi('module_by_type', 'isolate', 'log1p(n_hib)'), '-0.008')
gi = [i for i in G.index if 'x' in i.lower() or 'inter' in i.lower()]
add('interaction with gRodon', f'{A} Table 2', G.at[gi[0], 'estimate'] if gi else 'n/a', '0.024', note=f'source row {gi[:1]}')
# follow-ups
for i, lab, pv in [('F_MAG_hib_extra', 'hib_extra MAG', '0.013'), ('F_iso_hib_extra', 'hib_extra iso', '-0.011'), ('F_MAG_hpfshort', 'short HPF MAG', '0.013'),
                   ('F_iso_hpfshort', 'short HPF iso', '-0.006'), ('F_MAG_balon', 'Balon MAG', '0.008'), ('F_MAG_hib_extra_lowcont', 'hib_extra MAG low contam', '0.017'),
                   ('F_MAG_DORM_lowcont', 'DORM MAG low contam', '0.022'), ('F_MAG_DORM_comp98', 'DORM comp98', '0.007'),
                   ('F_iso_precise', 'iso precise', '-0.010'), ('F_iso_coarse', 'iso coarse', '-0.001'), ('F_MAG_precise', 'MAG precise', '0.014'),
                   ('F_MAG_coarse', 'MAG coarse', '0.042'), ('F_balon_Tmean_marine_iso', 'Balon~SST marine iso logit', '-0.34'),
                   ('F_balon_Tmean_land', 'Balon~T land logit', '0.16'), ('E_cyanob', 'Cyanobacteriota slope', '0.036'), ('E_balon_any', 'Balon~P0 logit', '0.071')]:
    tab = F if i.startswith('F_') else E
    add(lab, f'{A} results / Table 2', tab.at[i, 'estimate'], pv)
add('hib_extra MAG p (2e-5)', f'{A}', F.at['F_MAG_hib_extra', 'p'], '2e-5', 'match' if 1.5e-5 < F.at['F_MAG_hib_extra', 'p'] < 2.5e-5 else 'mismatch')
add('hib_extra iso p', f'{A}', F.at['F_iso_hib_extra', 'p'], '0.0015'); add('comp98 p', f'{A}', F.at['F_MAG_DORM_comp98', 'p'], '0.35')
add('comp98 n', f'{A}', F.at['F_MAG_DORM_comp98', 'n'], '906'); add('Balon marine iso p (4e-4)', f'{A}', F.at['F_balon_Tmean_marine_iso', 'p'], '4e-4', 'match')
add('Balon marine MAG p', f'{A}', F.at['F_balon_Tmean_marine_MAG', 'p'], '0.13'); add('Cyanobacteriota p', f'{A}', E.at['E_cyanob', 'p'], '0.004')
mi_ = c[(c.domain == 'marine') & (~mag)]
q = pd.qcut(mi_.sst_mean, 8); g = mi_.groupby(q, observed=True).balon_any.mean()
add('Balon % in marine isolates, coldest octile ("about 17 % below 5 °C")', f'{A}; {S}', 100 * g.iloc[0], '17', 'mismatch',
    f'17 % is the coldest octile, which spans -1.6 to 9.7 °C, not "below 5 °C"; coldest octile {g.index[0]} = {100*g.iloc[0]:.1f}%; all marine isolates with SST < 5 °C: {100*mi_[mi_.sst_mean < 5].balon_any.mean():.1f}% (n={(mi_.sst_mean < 5).sum()})')
add('Balon % warm ("about 5 % above 20 °C")', f'{A}; {S}', 100 * g.iloc[-1], '5', 'mismatch',
    f'warmest octile {g.index[-1]} = {100*g.iloc[-1]:.1f}%; all with SST > 20 °C: {100*mi_[mi_.sst_mean > 20].balon_any.mean():.1f}% (n={(mi_.sst_mean > 20).sum()}); octile curve is non-monotonic: ' + ', '.join(f'{100*v:.0f}' for v in g.values) + ' %; "5 %" matches only the 5th/7th octiles')
add('student: "λ = 0.94 means about 94 % of the pattern follows the family tree"', f'{S} §4.2', 'λ is a scaling of the Brownian covariance, not a variance share', '94 %', 'mismatch', 'interpretation error')
out = pd.DataFrame(rows)
out.to_csv('audit/number_check.tsv', sep='\t', index=False)
print(out.status.value_counts().to_string())
print(out[out.status != 'match'].to_string())
