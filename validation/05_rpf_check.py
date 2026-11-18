#!/usr/bin/env python3
"""Compare predicted n_rpf / n_sps (scripts/08) with literature-documented counts.

usage: python3 validation/05_rpf_check.py WORKDIR NAME [WORKDIR2 NAME2]
  (optional second run = Puspita 2013 tier, validation/rpf_labels_puspita.tsv)
   -> validation/rpf_check_table.tsv,
       validation/rpf_check_numbers.json, and the per-protein Rpf hits for the Bacillota discordances
"""
import sys, json
import pandas as pd

W, NAME = sys.argv[1], sys.argv[2]
l = pd.read_csv('validation/rpf_labels.tsv', sep='\t')
l['acc'] = [('RS_' if a.startswith('GCF') else 'GB_') + a for a in l.accession]
m = pd.read_csv(f'{W}/{NAME}_modules.tsv', sep='\t')
d = l.merge(m[['accession', 'n_proteins', 'n_rpf', 'n_rpfB', 'n_sps']], left_on='acc', right_on='accession',
            how='left', suffixes=('', '_m')).drop(columns=['accession_m'])
d['rpf_exact'] = d.expected_rpf.isna() | (d.n_rpf == d.expected_rpf)
d['rpf_presence_ok'] = d.expected_rpf.isna() & (d.n_rpf > 0) | ((d.expected_rpf > 0) == (d.n_rpf > 0))
d['sps_diff'] = d.n_sps - d.expected_sps
d.drop(columns=['acc']).to_csv('validation/rpf_check_table.tsv', sep='\t', index=False)

act = d[d.phylum == 'Actinomycetota']; non = d[d.phylum != 'Actinomycetota']
cnt = act[act.expected_rpf.notna()]
num = dict(
    n_actino=len(act), n_actino_count=len(cnt), n_actino_count_pos=int((cnt.expected_rpf > 0).sum()),
    actino_exact=int((cnt.n_rpf == cnt.expected_rpf).sum()),
    actino_presence_ok=int(((cnt.expected_rpf > 0) == (cnt.n_rpf > 0)).sum()),
    n_non=len(non), non_zero_ok=int((non.n_rpf == 0).sum()),
    non_false_pos=non[non.n_rpf > 0][['organism', 'strain', 'n_rpf', 'evidence']].to_dict('records'),
    sps_n=int(d.expected_sps.notna().sum() - (d.phylum != 'Bacillota').sum() * 0),
    sps_bacillota=int((d.phylum == 'Bacillota').sum()),
    sps_exact=int(((d.phylum == 'Bacillota') & (d.sps_diff == 0)).sum()),
    sps_off=d[(d.phylum == 'Bacillota') & (d.sps_diff != 0)][['organism', 'strain', 'expected_sps', 'n_sps']].to_dict('records'),
    sps_outside_bacillota_nonzero=int(((d.phylum != 'Bacillota') & (d.n_sps > 0)).sum()),
)
# protein-level detail for the Bacillota PF06737 hits
h = pd.read_csv(f'{W}/{NAME}_hits.tsv.gz', sep='\t')
fp = non[non.n_rpf > 0].acc.tolist()
det = h[h.accession.isin(fp) & (h.model == 'Transglycosylas')][['accession', 'contig', 'gene_idx', 'start', 'end', 'score', 'dom_score', 'pass']]
det['len_aa'] = (det.end - det.start + 1) // 3
thr = pd.read_csv('dictionary/thresholds.tsv', sep='\t').set_index('model_name')
num['pf06737_threshold'] = float(thr.threshold.get('Transglycosylas'))
num['bacillota_pf06737_hits'] = det.to_dict('records')
if len(sys.argv) > 4:
    from scipy.stats import spearmanr
    W2, N2 = sys.argv[3], sys.argv[4]
    p = pd.read_csv('validation/rpf_labels_puspita.tsv', sep='\t')
    p['acc'] = 'RS_' + p.accession
    m2 = pd.read_csv(f'{W2}/{N2}_modules.tsv', sep='\t')
    p = p.merge(m2[['accession', 'n_rpf', 'n_rpfB']], left_on='acc', right_on='accession', how='left', suffixes=('', '_m')).drop(columns=['accession_m', 'acc'])
    p.to_csv('validation/rpf_check_table_puspita.tsv', sep='\t', index=False)
    q = p[p.n_rpf.notna()]
    diff = q.n_rpf - q.expected_rpf
    num['tier2'] = dict(n=len(p), n_scored=len(q), exact=int((diff == 0).sum()), within1=int((diff.abs() <= 1).sum()),
                        over=int((diff > 0).sum()), under=int((diff < 0).sum()),
                        spearman=float(spearmanr(q.expected_rpf, q.n_rpf)[0]),
                        presence_TP=int(((q.expected_rpf > 0) & (q.n_rpf > 0)).sum()),
                        presence_FN=int(((q.expected_rpf > 0) & (q.n_rpf == 0)).sum()),
                        presence_FP=int(((q.expected_rpf == 0) & (q.n_rpf > 0)).sum()),
                        presence_TN=int(((q.expected_rpf == 0) & (q.n_rpf == 0)).sum()),
                        crosstab=pd.crosstab(q.expected_rpf, q.n_rpf).to_dict(),
                        largest_discordances=q.assign(diff=diff).reindex(diff.abs().sort_values(ascending=False).index)
                        .head(12)[['organism', 'accession', 'expected_rpf', 'n_rpf', 'diff']].to_dict('records'))
json.dump(num, open('validation/rpf_check_numbers.json', 'w'), indent=1, default=float)
print(json.dumps(num, indent=1, default=float))
