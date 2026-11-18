#!/usr/bin/env python3
"""Phase 5c: turn HMM hits into per-genome dormancy module scores, following dictionary/README.md.

usage: 08_aggregate_modules.py HITS.tsv.gz NPROT.tsv META.tsv OUT.tsv
  META must have accession and phylum (for the Sps rule).

Outputs per genome (counts are numbers of distinct proteins unless noted):
  spo_strict (0-33 distinct scored core genes), spo_relaxed, spo0A, sporulator (spo0A & strict>=25),
  spo_reduced (spo0A & 17-24); n_rpf, n_rpfB, n_sps (Bacillota only);
  hpf_long, hpf_short, n_rmf, n_balon, n_mrf, hib_systems (0-5), n_hib;
  n_rsh_long, n_sas, n_sah, n_dksA, n_sr (= long+SAS+SAH);
  hipA, rpoS; ta_toxins, ta_antitoxins, ta_pairs, ta_pairs_plasmid, ta_orphan_tox;
  n_transposase, n_integrase.
"""
import sys, numpy as np, pandas as pd

D = 'dictionary'
hits = pd.read_csv(sys.argv[1], sep='\t')
nprot = pd.read_csv(sys.argv[2], sep='\t')
meta = pd.read_csv(sys.argv[3], sep='\t')
out = sys.argv[4]
dic = pd.read_csv(f'{D}/dormancy_dictionary.tsv', sep='\t')
core = pd.read_csv(f'{D}/sporulation_core.tsv', sep='\t')

m2 = dic.dropna(subset=['model_name']).drop_duplicates('model_name').set_index('model_name')
# a model can occur in two modules (HipA_C: persistence + TA); keep all memberships
mod_of = dic.dropna(subset=['model_name']).groupby('model_name').module.apply(set).to_dict()
gene_of = {(r.module, r.model_name): r.gene for r in dic.dropna(subset=['model_name']).itertuples()}

hits['prot'] = hits.accession + '|' + hits.contig.astype(str) + '|' + hits.gene_idx.astype(str)
hits['len_aa'] = (hits.end - hits.start + 1) // 3
P = hits[hits['pass'] == 1]


def has(df, models):
    return set(df[df.model.isin(models)].prot)


rows = {}
# KOfam K06381 (SpoIID) also detects SpoIID-domain amidases outside Bacillota (Deinococcus,
# cyanobacteria) in the reference check, so spoIID is scored only by NCBIfam spore_II_D.
NONSPECIFIC = {'K06381'}
spor_models = {m for m, s in mod_of.items() if any(x.startswith('sporulation') for x in s)} - NONSPECIFIC
scored = core[core.in_score == 'yes']
scored_genes = scored.gene.tolist()

for acc, g in hits.groupby('accession', sort=False):
    r = {}
    gp = g[g['pass'] == 1]
    # ---- sporulation: competitive assignment over all sporulation hits (any score)
    sg = g[g.model.isin(spor_models)].sort_values('score', ascending=False)
    best = sg.drop_duplicates('prot').set_index('prot').model
    def spor_gene(m):
        for mod in ('sporulation_core', 'sporulation_extended'):
            if (mod, m) in gene_of:
                return gene_of[(mod, m)]
    sg = sg[sg.apply(lambda x: spor_gene(best[x.prot]) == spor_gene(x.model), axis=1)] if len(sg) else sg
    strict = {spor_gene(m) for m in sg[sg['pass'] == 1].model}
    relaxed = strict | {spor_gene(m) for m in sg[sg.pass_nc == 1].model}
    r['spo_strict'] = sum(gn in strict for gn in scored_genes)
    r['spo_relaxed'] = sum(gn in relaxed for gn in scored_genes)
    r['spo0A'] = int('spo0A' in strict)
    r['sporulator'] = int(r['spo0A'] and r['spo_strict'] >= 25)
    r['spo_reduced'] = int(r['spo0A'] and 17 <= r['spo_strict'] <= 24)
    # ---- resuscitation
    rpf = has(gp, ['Transglycosylas'])
    r['n_rpf'] = len(rpf)
    r['n_rpfB'] = len(rpf & has(gp, ['G5']) & has(gp, ['DUF348']))
    r['n_sps'] = len(has(gp, ['3D']) - rpf)
    # ---- hibernation
    s30 = has(gp, ['yfiA', 'Ribosomal_S30AE', 'K05808', 'K05809'])
    cterm = has(gp, ['Ribosom_S30AE_C'])
    r['hpf_long'] = len(s30 & cterm)
    r['hpf_short'] = len(s30 - cterm)
    r['n_rmf'] = len(has(gp, ['RMF', 'K03812']))
    r['n_balon'] = len(has(gp, ['AOC03_06830_fam_RHF', 'Rv2629_fam_RHF', 'MSMEG_1130_fam_RHF', 'Balon_mid']))
    r['n_mrf'] = len(has(gp, ['hiber_recruit']))
    hib = ['hpf_long', 'hpf_short', 'n_rmf', 'n_balon', 'n_mrf']
    r['hib_systems'] = sum(r[k] > 0 for k in hib)
    r['n_hib'] = sum(r[k] for k in hib)
    # ---- stringent response
    synth = has(gp, ['RelA_SpoT']); tgs = has(gp, ['TGS']); hd = has(gp, ['HD_4']); act = has(gp, ['ACT_4'])
    longrsh = has(gp, ['spoT_relA', 'K00951', 'K01139']) | (synth & tgs)
    L = g.drop_duplicates('prot').set_index('prot').len_aa
    sas = (has(gp, ['K07816']) | {p for p in synth - tgs - hd - act if L.get(p, 999) <= 350}) - longrsh
    sah = {p for p in hd - synth if L.get(p, 999) <= 300} - longrsh
    r['n_rsh_long'] = len(longrsh); r['n_sas'] = len(sas); r['n_sah'] = len(sah)
    r['n_dksA'] = len(has(gp, ['dksA']))
    r['n_sr'] = r['n_rsh_long'] + r['n_sas'] + r['n_sah']
    # ---- persistence / stress covariates
    r['hipA'] = len(has(gp, ['HipA_C', 'K07154']) & (has(gp, ['Couple_hipA']) | has(gp, ['K07154'])))
    r['rpoS'] = int(len(has(gp, ['K03087', 'rpoS_proteo'])) > 0)
    # ---- toxin-antitoxin: pairs within 3 genes on the same contig
    tox_models = {m for m, s in mod_of.items() if 'TA_typeII_toxin' in s}
    at_models = {m for m, s in mod_of.items() if 'TA_typeII_antitoxin' in s}
    bestall = g.sort_values('score', ascending=False).drop_duplicates('prot').set_index('prot').model
    excl = {p for p, m in bestall.items() if m in ('spore_V_T', 'K04769')}
    tox = gp[gp.model.isin(tox_models) & ~gp.prot.isin(excl)].drop_duplicates('prot')
    at = gp[gp.model.isin(at_models) & ~gp.prot.isin(excl)].drop_duplicates('prot')
    at = at[~at.prot.isin(tox.prot)]
    pairs, pl = 0, 0
    if len(tox) and len(at):
        for t in tox.itertuples():
            a = at[(at.contig == t.contig) & ((at.gene_idx - t.gene_idx).abs() <= 3)]
            if len(a):
                pairs += 1; pl += int(t.plasmid == 1)
    r['ta_toxins'] = len(tox); r['ta_antitoxins'] = len(at); r['ta_pairs'] = pairs
    r['ta_pairs_plasmid'] = pl; r['ta_orphan_tox'] = len(tox) - pairs
    # ---- mobile elements
    tp_models = [m for m in m2.index if str(m2.at[m, 'gene']).startswith('transposase:')]
    r['n_transposase'] = len(has(gp, tp_models))
    integ = has(gp, ['Phage_integrase', 'Resolvase', 'Recombinase'])
    xer = {p for p, m in bestall.items() if m in ('recomb_XerC', 'recomb_XerD')}
    r['n_integrase'] = len(integ - xer)
    rows[acc] = r

res = pd.DataFrame.from_dict(rows, orient='index')
res.index.name = 'accession'
res = nprot.set_index('accession')[['n_proteins']].join(res, how='left').fillna(0)
ph = meta.set_index('accession').phylum.reindex(res.index).astype(str)
res.loc[~ph.str.startswith('Bacillota'), 'n_sps'] = 0
res = res.astype({c: int for c in res.columns})
res.to_csv(out, sep='\t')
print(res.describe().T[['mean', 'min', 'max']].to_string())
