#!/usr/bin/env python3
"""Clade-level generalisation of the sporulation call (added 18 Nov 2026).
The dictionary and the pre-registered rule (Spo0A + >= 25 of 33 core genes) were fixed BEFORE these labels were
assembled, so no genome here was used to build or tune them. Two further checks address phylogenetic
non-independence of the labelled genomes:
  1. Genus-clustered bootstrap (2,000 resamples of whole genera) for AUC and for the pre-registered rule's
     sensitivity and specificity.
  2. Leave-one-order-out: choose the Youden-optimal spo_strict threshold on all other orders, then evaluate it
     (and the fixed pre-registered rule) on the held-out order.
usage (repo root): python3 validation/06_clade_holdout.py -> validation/clade_holdout.tsv, validation/clade_holdout_summary.tsv"""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

d = pd.read_csv('validation/per_genome_scores.tsv', sep='\t')
d = d[(d.set == 'calibration') & d.label.isin([0, 1])].copy()
d['genus'] = d.organism.str.split().str[0]
d['rule'] = ((d.spo0A > 0) & (d.spo_strict >= 25)).astype(int)
rng = np.random.default_rng(20261118)

def metrics(x, pred):
    tp = ((pred == 1) & (x.label == 1)).sum(); fn = ((pred == 0) & (x.label == 1)).sum()
    tn = ((pred == 0) & (x.label == 0)).sum(); fp = ((pred == 1) & (x.label == 0)).sum()
    return tp / max(tp + fn, 1), tn / max(tn + fp, 1)

# 1. genus-clustered bootstrap
genera = d.genus.unique(); groups = {g: d[d.genus == g] for g in genera}
boot = []
for _ in range(2000):
    s = pd.concat([groups[g] for g in rng.choice(genera, len(genera), replace=True)])
    if s.label.nunique() < 2: continue
    se, sp = metrics(s, s.rule)
    boot.append((roc_auc_score(s.label, s.spo_strict), se, sp))
boot = np.array(boot)
se0, sp0 = metrics(d, d.rule)
summ = [dict(check='genus-clustered bootstrap (2,000)', n_genomes=len(d), n_genera=len(genera),
             AUC=roc_auc_score(d.label, d.spo_strict), AUC_lo=np.percentile(boot[:, 0], 2.5), AUC_hi=np.percentile(boot[:, 0], 97.5),
             sens_rule=se0, sens_lo=np.percentile(boot[:, 1], 2.5), sens_hi=np.percentile(boot[:, 1], 97.5),
             spec_rule=sp0, spec_lo=np.percentile(boot[:, 2], 2.5), spec_hi=np.percentile(boot[:, 2], 97.5))]

# 2. leave-one-order-out threshold selection
rows = []
for o in sorted(d['order'].dropna().unique()):
    tr, te = d[d['order'] != o], d[d['order'] == o]
    best = max(range(0, 34), key=lambda t: sum(metrics(tr, ((tr.spo0A > 0) & (tr.spo_strict >= t)).astype(int))) - 1)
    pred = ((te.spo0A > 0) & (te.spo_strict >= best)).astype(int)
    se, sp = metrics(te, pred); se_r, sp_r = metrics(te, te.rule)
    rows.append(dict(held_out_order=o, n=len(te), n_pos=int(te.label.sum()), threshold_chosen_on_others=best,
                     acc_heldout=(pred == te.label).mean(), sens_heldout=se if te.label.sum() else np.nan,
                     spec_heldout=sp if (te.label == 0).sum() else np.nan,
                     acc_prereg_rule=(te.rule == te.label).mean()))
lo = pd.DataFrame(rows)
lo.to_csv('validation/clade_holdout.tsv', sep='\t', index=False, float_format='%.4g')
w = lo.n / lo.n.sum()
summ.append(dict(check='leave-one-order-out (weighted by n)', n_genomes=int(lo.n.sum()), n_genera=len(lo),
                 AUC=np.nan, sens_rule=np.nan, spec_rule=np.nan,
                 acc_heldout=(lo.acc_heldout * w).sum(), acc_prereg_rule=(lo.acc_prereg_rule * w).sum(),
                 thresholds_chosen=','.join(map(str, sorted(lo.threshold_chosen_on_others.unique())))))
pd.DataFrame(summ).to_csv('validation/clade_holdout_summary.tsv', sep='\t', index=False, float_format='%.4g')
print(pd.DataFrame(summ).round(3).T); print(lo.round(3).to_string())
