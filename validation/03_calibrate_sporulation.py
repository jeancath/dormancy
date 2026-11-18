#!/usr/bin/env python3
"""Calibrate the genome-based sporulation calls against independent phenotype labels.

usage: python3 validation/03_calibrate_sporulation.py WORKDIR NAME
  reads validation/sporulation_labels.tsv, WORKDIR/NAME_modules.tsv, WORKDIR/NAME_hits.tsv.gz
  writes validation/calibration_metrics.tsv, validation/misclassified.tsv,
         validation/per_genome_scores.tsv, validation/fig_sporulation_roc.png,
         validation/calibration_numbers.json (numbers quoted in sporulation_calibration.md)

Scores come from scripts/08_aggregate_modules.py unchanged. The per-gene presence needed to list
missing genes is re-derived here with the same competitive-assignment rule, and the script asserts
that this re-derivation reproduces 08's spo_strict / spo_relaxed for every genome, so it is a
diagnostic view of 08's result, not a second scorer.
"""
import sys, json
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve, matthews_corrcoef
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

W, NAME = sys.argv[1], sys.argv[2]
B = 2000
rng = np.random.default_rng(20261117)
BLUE, ORANGE = '#2a78d6', '#eb6834'

lab = pd.read_csv('validation/sporulation_labels.tsv', sep='\t')
lab['acc'] = [('RS_' if a.startswith('GCF') else 'GB_') + a for a in lab.accession]
mod = pd.read_csv(f'{W}/{NAME}_modules.tsv', sep='\t').rename(columns={'accession': 'acc'}).drop_duplicates('acc')
d = lab.merge(mod, on='acc', how='left')
d['scored'] = d.spo_strict.notna()
print('labelled:', len(d), 'scored:', int(d.scored.sum()))

# ---------------------------------------------------------------- per-gene diagnostics
dic = pd.read_csv('dictionary/dormancy_dictionary.tsv', sep='\t').dropna(subset=['model_name'])
core = pd.read_csv('dictionary/sporulation_core.tsv', sep='\t')
scored_genes = core[core.in_score == 'yes'].gene.tolist()
gene_of = {}
for r in dic.itertuples():
    if r.module in ('sporulation_core', 'sporulation_extended'):
        gene_of.setdefault(r.model_name, r.gene)
spor_models = set(gene_of) - {'K06381'}
hits = pd.read_csv(f'{W}/{NAME}_hits.tsv.gz', sep='\t')
hits['prot'] = hits.accession + '|' + hits.contig.astype(str) + '|' + hits.gene_idx.astype(str)
sg = hits[hits.model.isin(spor_models)].sort_values('score', ascending=False)
best = sg.drop_duplicates('prot').set_index('prot').model.map(gene_of)
sg = sg[sg.model.map(gene_of).values == best.reindex(sg.prot).values]
sg['gene'] = sg.model.map(gene_of)
strict = sg[sg['pass'] == 1].groupby('accession').gene.apply(set)
nc = sg[sg.pass_nc == 1].groupby('accession').gene.apply(set)
thr = pd.read_csv('dictionary/thresholds.tsv', sep='\t').set_index('model_name')


def genes(acc, kind):
    s = strict.get(acc, set())
    return s | nc.get(acc, set()) if kind == 'relaxed' else s


chk = d[d.scored]
re_strict = chk.acc.map(lambda a: sum(g in genes(a, 'strict') for g in scored_genes))
re_relax = chk.acc.map(lambda a: sum(g in genes(a, 'relaxed') for g in scored_genes))
assert (re_strict.values == chk.spo_strict.values).all(), 'diagnostic re-derivation != 08 spo_strict'
assert (re_relax.values == chk.spo_relaxed.values).all(), 'diagnostic re-derivation != 08 spo_relaxed'


def sigk_fragments(acc):
    """Sub-threshold SigK-model hits on >= 2 proteins: candidate skin-interrupted sigK."""
    h = hits[(hits.accession == acc) & hits.model.isin(['spore_sigmaK', 'PRK05803.1'])]
    h = h[h.prot.map(lambda p: best.get(p) == 'sigK')]
    return h[h['pass'] == 0].prot.nunique()


def near_miss(acc, gene):
    """Best score of the gene's models in this genome vs its cutoff (only best-assigned proteins)."""
    h = sg[(sg.accession == acc) & (sg.gene == gene)]
    if not len(h):
        return ''
    t = h.iloc[0]
    cut = thr.threshold.get(t.model, np.nan)
    return f'{gene}:{t.model} {t.score:.0f}/{cut:.0f}'


# ---------------------------------------------------------------- metrics
def lineage(r):
    if r.family == 'Bacillaceae':
        return 'Bacillaceae'
    if r['class'] == 'Clostridia':
        return 'Clostridia'
    if r['class'] == 'Erysipelotrichia':
        return 'Erysipelotrichia'
    if r['class'] == 'Negativicutes':
        return 'Negativicutes'
    if r.order == 'Lactobacillales':
        return 'Lactobacillales'
    if r.family == 'Paenibacillaceae':
        return 'Paenibacillaceae'
    if r['class'] == 'Mollicutes':
        return 'Mollicutes'
    if r['class'] == 'Tissierellia':
        return 'Tissierellia'
    return 'other Bacillales'


RULES = {
    'Spo0A & >=25 (pre-registered)': lambda x: (x.spo0A == 1) & (x.spo_strict >= 25),
    'Spo0A & >=17': lambda x: (x.spo0A == 1) & (x.spo_strict >= 17),
    'Spo0A & >=20': lambda x: (x.spo0A == 1) & (x.spo_strict >= 20),
    'Spo0A & >=27': lambda x: (x.spo0A == 1) & (x.spo_strict >= 27),
    'Spo0A only': lambda x: x.spo0A == 1,
    'Spo0A & relaxed >=25': lambda x: (x.spo0A == 1) & (x.spo_relaxed >= 25),
}


def conf(y, p):
    y, p = np.asarray(y, bool), np.asarray(p, bool)
    tp, fn = int((y & p).sum()), int((y & ~p).sum())
    fp, tn = int((~y & p).sum()), int((~y & ~p).sum())
    sens = tp / (tp + fn) if tp + fn else np.nan
    spec = tn / (tn + fp) if tn + fp else np.nan
    prec = tp / (tp + fp) if tp + fp else np.nan
    mcc = matthews_corrcoef(y, p) if len(set(y)) > 1 and len(set(p)) > 1 else np.nan
    acc = (tp + tn) / len(y)
    return dict(n=len(y), TP=tp, FN=fn, FP=fp, TN=tn, sensitivity=sens, specificity=spec,
                precision=prec, MCC=mcc, accuracy=acc)


def boot(y, x, f):
    y, x = np.asarray(y), np.asarray(x)
    out = []
    n = len(y)
    for _ in range(B):
        i = rng.integers(0, n, n)
        if len(set(y[i])) < 2:
            continue
        out.append(f(y[i], x[i]))
    return np.nanpercentile(out, [2.5, 97.5])


def youden(y, s):
    best_k, best_j = None, -1
    for k in range(0, 34):
        p = s >= k
        se = (p & y).sum() / y.sum(); sp = (~p & ~y).sum() / (~y).sum()
        if se + sp - 1 > best_j:
            best_j, best_k = se + sp - 1, k
    return best_k, best_j


rows, numbers = [], {}
cal = d[(d.set == 'calibration') & d.scored].copy()
cal['lineage'] = cal.apply(lineage, axis=1)
y = cal.label.values.astype(bool)
numbers['n_calibration'] = int(len(cal)); numbers['n_pos'] = int(y.sum()); numbers['n_neg'] = int((~y).sum())
numbers['n_not_scored'] = int(((d.set == 'calibration') & ~d.scored).sum())

for score in ('spo_strict', 'spo_relaxed'):
    s = cal[score].values
    auc = roc_auc_score(y, s)
    lo, hi = boot(y, s, roc_auc_score)
    rows.append(dict(subset='calibration (all Bacillota)', metric=f'AUC {score}', value=auc, ci_lo=lo, ci_hi=hi, n=len(y)))
    numbers[f'auc_{score}'] = [auc, lo, hi]
    k, j = youden(y, s)
    rows.append(dict(subset='calibration (all Bacillota)', metric=f'Youden-optimal threshold {score} (score>=k)', value=k, n=len(y)))
    rows.append(dict(subset='calibration (all Bacillota)', metric=f'Youden J at optimum {score}', value=j, n=len(y)))
    # Youden with Spo0A required
    kbest = max(range(34), key=lambda kk: (lambda p: (p & y).sum() / y.sum() + (~p & ~y).sum() / (~y).sum())(
        (cal.spo0A.values == 1) & (s >= kk)))
    rows.append(dict(subset='calibration (all Bacillota)', metric=f'Youden-optimal k for Spo0A & {score}>=k', value=kbest, n=len(y)))
    numbers[f'youden_{score}'] = [k, j, kbest]

for name, f in RULES.items():
    p = f(cal).values
    c = conf(y, p)
    for mname in ('sensitivity', 'specificity', 'precision', 'MCC'):
        fn = {'sensitivity': lambda yy, pp: (yy & pp).sum() / yy.sum(),
              'specificity': lambda yy, pp: (~yy & ~pp).sum() / (~yy).sum(),
              'precision': lambda yy, pp: (yy & pp).sum() / pp.sum() if pp.sum() else np.nan,
              'MCC': lambda yy, pp: matthews_corrcoef(yy, pp)}[mname]
        lo, hi = boot(y, p, fn)
        rows.append(dict(subset='calibration (all Bacillota)', metric=f'{name}: {mname}', value=c[mname],
                         ci_lo=lo, ci_hi=hi, n=c['n'], TP=c['TP'], FN=c['FN'], FP=c['FP'], TN=c['TN']))
    numbers[name] = c

pre = RULES['Spo0A & >=25 (pre-registered)']
# label-quality sensitivity (pre-specified subsets; no genome removed from the primary metrics)
lab_q = {'tier A (IJSEM description present)': cal.label_source == 'Madin2020_tierA',
         '>=2 curated sources agree': cal.evidence.str.count('=') - cal.evidence.str.count('protraits=') >= 2,
         'complete/chromosome assemblies': cal.assembly_level.isin(['Complete Genome', 'Chromosome'])}
for qn, mask in lab_q.items():
    g = cal[mask.values]
    yy = g.label.values.astype(bool)
    c = conf(yy, pre(g).values)
    auc = roc_auc_score(yy, g.spo_strict)
    lo, hi = boot(yy, g.spo_strict.values, roc_auc_score)
    rows.append(dict(subset=f'label quality: {qn}', metric='AUC spo_strict', value=auc, ci_lo=lo, ci_hi=hi, n=len(g)))
    for k in ('sensitivity', 'specificity', 'MCC'):
        rows.append(dict(subset=f'label quality: {qn}', metric=f'Spo0A & >=25 (pre-registered): {k}', value=c[k], n=len(g),
                         TP=c['TP'], FN=c['FN'], FP=c['FP'], TN=c['TN']))
    numbers[f'quality: {qn}'] = dict(c, auc=[auc, lo, hi])

# per-lineage: pre-registered rule, plus AUC where both classes present
pre = RULES['Spo0A & >=25 (pre-registered)']
lin_tab = []
for ln, g in cal.groupby('lineage'):
    yy = g.label.values.astype(bool)
    c = conf(yy, pre(g).values)
    auc = roc_auc_score(yy, g.spo_strict) if len(set(yy)) == 2 else np.nan
    lin_tab.append(dict(lineage=ln, n=len(g), n_pos=int(yy.sum()), n_neg=int((~yy).sum()), AUC_strict=auc,
                        **{k: c[k] for k in ('TP', 'FN', 'FP', 'TN', 'sensitivity', 'specificity', 'MCC')},
                        median_strict_pos=g[yy].spo_strict.median(), median_strict_neg=g[~yy].spo_strict.median()))
    for k in ('sensitivity', 'specificity', 'MCC'):
        rows.append(dict(subset=f'lineage: {ln}', metric=f'Spo0A & >=25 (pre-registered): {k}', value=c[k], n=len(g),
                         TP=c['TP'], FN=c['FN'], FP=c['FP'], TN=c['TN']))
    rows.append(dict(subset=f'lineage: {ln}', metric='AUC spo_strict', value=auc, n=len(g)))
lin_tab = pd.DataFrame(lin_tab)
numbers['lineage'] = lin_tab.to_dict('records')

# other sets
for st in ('browne_experimental', 'specificity', 'qc_nonendospore'):
    g = d[(d.set == st) & d.scored]
    if not len(g):
        continue
    yy = g.label.values.astype(bool)
    c = conf(yy, pre(g).values)
    rows.append(dict(subset=st, metric='Spo0A & >=25 (pre-registered): confusion', value=c['accuracy'], n=c['n'],
                     TP=c['TP'], FN=c['FN'], FP=c['FP'], TN=c['TN']))
    rows.append(dict(subset=st, metric='genomes with spo_strict > 0', value=int((g.spo_strict > 0).sum()), n=len(g)))
    rows.append(dict(subset=st, metric='max spo_strict', value=g.spo_strict.max(), n=len(g)))
    numbers[st] = dict(c, n_nonzero=int((g.spo_strict > 0).sum()), max=float(g.spo_strict.max()),
                       nonzero=g[g.spo_strict > 0][['organism', 'phylum', 'spo_strict', 'spo0A']].to_dict('records'))

met = pd.DataFrame(rows)
met.to_csv('validation/calibration_metrics.tsv', sep='\t', index=False, float_format='%.4g')
lin_tab.to_csv('validation/lineage_breakdown.tsv', sep='\t', index=False, float_format='%.3g')

# ---------------------------------------------------------------- misclassified genomes
mis = []
for st in ('calibration', 'browne_experimental', 'specificity'):
    g = d[(d.set == st) & d.scored]
    pred = pre(g).values
    for r, pp in zip(g.itertuples(), pred):
        if bool(pp) == bool(r.label):
            continue
        sgs, rgs = genes(r.acc, 'strict'), genes(r.acc, 'relaxed')
        missing = [x for x in scored_genes if x not in sgs]
        mis.append(dict(set=st, organism=r.organism, accession=r.accession, family=r.family,
                        lineage=lineage(lab.loc[lab.accession == r.accession].iloc[0]),
                        label=r.label, error='FN' if r.label else 'FP', spo0A=int(r.spo0A),
                        spo_strict=int(r.spo_strict), spo_relaxed=int(r.spo_relaxed),
                        assembly_level=r.assembly_level, label_source=r.label_source,
                        n_missing=len(missing),
                        missing_strict=','.join(missing) if len(missing) <= 20 else f'{len(missing)} genes (all but: ' + ','.join(x for x in scored_genes if x in sgs) + ')',
                        rescued_by_pass_nc=','.join(x for x in missing if x in rgs),
                        below_cutoff_hits='; '.join(filter(None, (near_miss(r.acc, x) for x in missing if x in rgs or x == 'sigK'))),
                        sigK_fragments=sigk_fragments(r.acc) if 'sigK' in missing else 0,
                        evidence=r.evidence))
mis = pd.DataFrame(mis).sort_values(['set', 'error', 'spo_strict'], ascending=[True, True, False])
mis.to_csv('validation/misclassified.tsv', sep='\t', index=False)
numbers['n_misclassified'] = mis.groupby(['set', 'error']).size().to_dict() if len(mis) else {}
numbers['n_misclassified'] = {f'{k[0]}:{k[1]}': int(v) for k, v in numbers['n_misclassified'].items()}

keep = ['organism', 'strain', 'accession', 'set', 'label', 'label_source', 'family', 'class', 'order',
        'assembly_level', 'type_material', 'n_proteins', 'spo0A', 'spo_strict', 'spo_relaxed', 'sporulator', 'spo_reduced']
d[d.scored][keep].to_csv('validation/per_genome_scores.tsv', sep='\t', index=False)

# gene-level: how often each scored gene is missing in label-positive genomes (strict)
pos = cal[cal.label == 1]
gm = pd.DataFrame({'gene': scored_genes,
                   'missing_in_pos_strict': [sum(gg not in genes(a, 'strict') for a in pos.acc) for gg in scored_genes],
                   'missing_in_pos_relaxed': [sum(gg not in genes(a, 'relaxed') for a in pos.acc) for gg in scored_genes],
                   'present_in_neg_strict': [sum(gg in genes(a, 'strict') for a in cal[cal.label == 0].acc) for gg in scored_genes]})
gm['n_pos'] = len(pos); gm['n_neg'] = int((cal.label == 0).sum())
gm.sort_values('missing_in_pos_strict', ascending=False).to_csv('validation/gene_missingness.tsv', sep='\t', index=False)

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), facecolor='white')
ax = axes[0]
for score, col, ls in (('spo_strict', BLUE, '-'), ('spo_relaxed', ORANGE, '--')):
    fpr, tpr, _ = roc_curve(y, cal[score])
    a, lo, hi = numbers[f'auc_{score}']
    ax.plot(fpr, tpr, color=col, lw=2, ls=ls, label=f'{score}  AUC {a:.3f} [{lo:.3f}–{hi:.3f}]')
c = numbers['Spo0A & >=25 (pre-registered)']
ax.plot(1 - c['specificity'], c['sensitivity'], 'o', color=BLUE, ms=8, mfc='white', mew=2,
        label=f'pre-registered rule (Spo0A & ≥25)')
ax.plot([0, 1], [0, 1], color='#bbbbbb', lw=1)
ax.set_xlabel('1 − specificity'); ax.set_ylabel('sensitivity')
ax.set_title(f'Bacillota, n = {len(cal)} ({int(y.sum())} sporeformers)', fontsize=10)
ax.legend(loc='lower right', fontsize=8, frameon=False)
ax.set_xlim(-0.01, 1); ax.set_ylim(0, 1.01)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
ax = axes[1]
bins = np.arange(-0.5, 34.5, 1)
ax.hist(cal[cal.label == 1].spo_strict, bins=bins, color=BLUE, alpha=0.85, label='sporeformer (label)')
ax.hist(cal[cal.label == 0].spo_strict, bins=bins, color=ORANGE, alpha=0.65, label='non-sporeformer (label)')
ax.axvline(24.5, color='#444444', lw=1, ls=':')
ax.text(24.7, ax.get_ylim()[1] * 0.92, '≥25', fontsize=8, color='#444444')
ax.set_xlabel('spo_strict (scored core genes, 0–33)'); ax.set_ylabel('genomes')
ax.legend(loc='upper center', fontsize=8, frameon=False)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig('validation/fig_sporulation_roc.png', dpi=200, facecolor='white')

json.dump(numbers, open('validation/calibration_numbers.json', 'w'), indent=1, default=float)
print(met.to_string())
print(lin_tab.to_string())
print(mis[['set', 'organism', 'error', 'spo0A', 'spo_strict', 'spo_relaxed', 'missing_strict']].to_string() if len(mis) else 'no errors')
