#!/usr/bin/env python3
"""step 1: select genomes of species with >= 2 QC-passing genomes in the same domain.
Domain rules copied exactly from scripts/04_select_dataset.py (marine: marine habitat, not deep,
sst_dist_km <= 50; rescued water_unspecified on an ocean cell and not a land cell; land: soil or
freshwater on a valid tair cell). A species present in both domains is kept in the domain with more
genomes (ties dropped), as in 04. Additional filters per the replication plan: completeness >= 90,
contamination < 5.
usage (repo root): python3 replication_within_species/01_select.py OUTDIR
"""
import sys, os, numpy as np, pandas as pd
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
hab = pd.read_csv('data_small/habitat_classified.tsv.gz', sep='\t', low_memory=False)
env = pd.read_csv('data_small/env_variability.tsv.gz', sep='\t')
d = hab.merge(env.drop(columns=['lat', 'lon']), on='accession', how='left')
marine = (d.habitat == 'marine') & (~d.is_deep) & (d.sst_dist_km <= 50)
rescued = (d.habitat == 'water_unspecified') & (d.sst_dist_km == 0) & (d.tair_dist_km > 0) & (~d.is_deep)
land = d.habitat.isin(['soil', 'freshwater']) & (d.tair_dist_km == 0)
d['domain'] = np.select([marine | rescued, land], ['marine', 'land'], default='')
d = d[(d.domain != '') & (d.checkm2_completeness >= 90) & (d.checkm2_contamination < 5)].copy()
cnt = d.groupby(['species', 'domain']).size().unstack(fill_value=0)
tie = cnt[(cnt.marine > 0) & (cnt.land > 0) & (cnt.marine == cnt.land)].index
keep_dom = pd.Series(np.where(cnt.marine >= cnt.land, 'marine', 'land'), index=cnt.index)
d = d[~d.species.isin(tie)]
d = d[d.domain.values == keep_dom.reindex(d.species).values]
d['T_sd_total'] = np.where(d.domain == 'marine', d.sst_sd_total, d.tair_sd_total)
d = d[np.isfinite(d.T_sd_total)]
n = d.groupby('species').accession.transform('size')
d = d[n >= 2].copy()
d['n_sp'] = d.groupby('species').accession.transform('size')
main = pd.read_csv('data_small/analysis_set.tsv.gz', sep='\t', usecols=['accession', 'pilot'])
d['in_main'] = d.accession.isin(main.accession).astype(int)
d = d.merge(main.rename(columns={'pilot': 'main_pilot'}), on='accession', how='left')
d['is_mag'] = d.is_mag.astype(bool).astype(int)
mixed = d.groupby('species').is_mag.agg(lambda v: v.min() == 0 and v.max() == 1)
d['mixed_species'] = d.species.map(mixed).astype(int)
d.to_csv(f'{out}/selected_genomes.tsv.gz', sep='\t', index=False, compression='gzip')
print('genomes', len(d), 'species', d.species.nunique(), 'already annotated', d.in_main.sum(),
      'to annotate', (d.in_main == 0).sum())
print(d.groupby('domain').agg(genomes=('accession', 'size'), species=('species', 'nunique')))
print('mixed MAG+isolate species:', mixed.sum(), 'genomes in them', d.mixed_species.sum())
print('MAGs', d.is_mag.sum(), 'isolates', (1 - d.is_mag).sum())
print('genomes per species quantiles', d.groupby('species').size().quantile([.5, .9, .99, 1]).to_dict())
print('largest species', d.groupby('species').size().sort_values().tail(8).to_dict())
