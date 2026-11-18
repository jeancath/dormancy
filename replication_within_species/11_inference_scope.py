#!/usr/bin/env python3
"""Added 18 Nov 2026: the actual scope of inference of the within-species analysis (HQ sample, species fixed
effects). With species fixed effects, a slope is identified only from variation in P0 among genomes of the SAME
species, so the counts that matter are the species (and genomes) with within-species P0 variation, by genome source.
usage (repo root): python3 replication_within_species/11_inference_scope.py -> replication_within_species/inference_scope.tsv"""
import pandas as pd, numpy as np
d = pd.read_csv('replication_within_species/data/within_species_table_hq.tsv.gz', sep='\t', low_memory=False)
d['is_mag'] = d.is_mag.astype(str).str.lower().isin(['1', 'true'])
d['P0w'] = d.P0 - d.groupby('species').P0.transform('mean')
def varying(g):  # >= 2 genomes at different P0 values
    return g.P0.round(6).nunique() >= 2
R = []
def add(k, v, note=''): R.append(dict(quantity=k, value=v, note=note))
sp = d.groupby('species')
add('genomes', len(d)); add('species', d.species.nunique())
add('MAG genomes', int(d.is_mag.sum())); add('isolate genomes', int((~d.is_mag).sum()))
vs = sp.apply(varying); add('species with within-species P0 variation', int(vs.sum()), 'contribute to any within-species slope')
for lab, m in (('MAG', True), ('isolate', False)):
    s = d[d.is_mag == m].groupby('species')
    v = s.apply(varying); n = s.size()
    add(f'species with >= 2 {lab} genomes', int((n >= 2).sum()))
    add(f'species with {lab} genomes at >= 2 P0 values', int(v.sum()), f'identify the within-species {lab} slope')
    add(f'{lab} genomes in those species', int(n[v[v].index].sum()))
mix = sp.is_mag.agg(lambda x: x.any() and (~x).any()); mixs = mix[mix].index
add('species with both MAG and isolate genomes', len(mixs))
add('genomes in those species', int(d.species.isin(mixs).sum()))
mm = d[d.species.isin(mixs)]
add('  of which MAG genomes', int(mm.is_mag.sum())); add('  of which isolate genomes', int((~mm.is_mag).sum()))
both_vary = sum(1 for s_, g in mm.groupby('species') if varying(g[g.is_mag]) and varying(g[~g.is_mag]))
add('mixed species with P0 variation within BOTH sources', both_vary)
# share of the within-species P0 variance (sum of squared deviations) by species type
tot = (d.P0w ** 2).sum()
kind = d.species.map(lambda s: 'mixed' if s in set(mixs) else None)
kind = kind.fillna(d.groupby('species').is_mag.transform(lambda x: 'MAG-only' if x.all() else 'isolate-only'))
for k in ('MAG-only', 'isolate-only', 'mixed'):
    add(f'share of within-species P0 variation from {k} species', round(float((d.P0w[kind == k] ** 2).sum() / tot), 3))
o = pd.DataFrame(R); o['value'] = o.value.map(lambda v: f'{v:g}'); o.to_csv('replication_within_species/inference_scope.tsv', sep='\t', index=False); print(o.to_string(index=False))

# ---- coverage table by species type (MAG-only, isolate-only, mixed)
T = []
for k in ('MAG-only', 'isolate-only', 'mixed'):
    ss = sorted(set(d.species[kind == k])); dd = d[d.species.isin(ss)]
    var_any = sum(1 for s_, g in dd.groupby('species') if varying(g))
    var_mag = sum(1 for s_, g in dd.groupby('species') if varying(g[g.is_mag]))
    var_iso = sum(1 for s_, g in dd.groupby('species') if varying(g[~g.is_mag]))
    var_both = sum(1 for s_, g in dd.groupby('species') if varying(g[g.is_mag]) and varying(g[~g.is_mag]))
    T.append(dict(species_type=k, species=len(ss), genomes=len(dd), MAG_genomes=int(dd.is_mag.sum()), isolate_genomes=int((~dd.is_mag).sum()),
                  species_with_P0_variation=var_any, species_with_P0_variation_among_MAGs=var_mag,
                  species_with_P0_variation_among_isolates=var_iso, species_with_P0_variation_in_both_sources=var_both,
                  share_of_within_species_P0_variation=round(float((dd.P0w ** 2).sum() / tot), 3)))
T = pd.DataFrame(T); T.loc[len(T)] = ['all'] + [T[c].sum() for c in T.columns[1:]]
T.to_csv('replication_within_species/inference_scope_by_type.tsv', sep='\t', index=False); print(T.to_string(index=False))
