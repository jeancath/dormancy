#!/usr/bin/env python3
"""EXTENSION beyond the replication plan (labelled as such everywhere): add MIMAG medium-quality genomes
(CheckM2 completeness >= 50 and contamination < 10, not already high quality) to the within-species
set. Inputs: GTDB r232 bac120 metadata (classified with scripts/01_habitat_classify.py rules, saved
as /data/gtdb/mq_classified.tsv.gz) and site variability computed with scripts/02_env_variability.py
on the NOAA grids (validated: 300 committed genomes reproduce sd_total exactly). Domain rules and the
majority-domain rule as in 01_select.py. Output: the combined HQ+MQ selection.
usage (repo root): python3 replication_within_species/05_select_mq.py /data/rep /data/mq"""
import sys, numpy as np, pandas as pd
W, Q = sys.argv[1:3]
hq = pd.read_csv(f'{W}/selected_genomes.tsv.gz', sep='\t', low_memory=False)
hab = pd.read_csv('data_small/habitat_classified.tsv.gz', sep='\t', low_memory=False)
env = pd.read_csv('data_small/env_variability.tsv.gz', sep='\t')
H = hab.merge(env.drop(columns=['lat', 'lon']), on='accession', how='left'); H['quality'] = 'HQ'
mq = pd.read_csv('/data/gtdb/mq_classified.tsv.gz', sep='\t', low_memory=False)
menv = pd.read_csv(f'{Q}/env_variability_mq.tsv.gz', sep='\t').drop_duplicates('accession')
M = mq.merge(menv.drop(columns=['lat', 'lon']), on='accession', how='inner'); M['quality'] = 'MQ'
M['species'] = M.species.str.replace('s__', '', regex=False)
M['phylum'] = M.gtdb_taxonomy.str.split(';').str[1].str.replace('p__', '', regex=False)
meta = pd.read_csv('/data/gtdb/bac120_metadata.tsv.gz', sep='\t', low_memory=False,
                   usecols=['accession', 'genome_size', 'ncbi_genome_category', 'ncbi_lat_lon', 'gtdb_taxonomy', 'checkm2_completeness'])
M = M.merge(meta[['accession', 'genome_size']], on='accession', how='left')
d = pd.concat([H, M], ignore_index=True)
marine = (d.habitat == 'marine') & (~d.is_deep.astype(bool)) & (d.sst_dist_km <= 50)
rescued = (d.habitat == 'water_unspecified') & (d.sst_dist_km == 0) & (d.tair_dist_km > 0) & (~d.is_deep.astype(bool))
land = d.habitat.isin(['soil', 'freshwater']) & (d.tair_dist_km == 0)
d['domain'] = np.select([marine | rescued, land], ['marine', 'land'], default='')
d = d[d.domain != ''].copy()
cnt = d.groupby(['species', 'domain']).size().unstack(fill_value=0)
tie = cnt[(cnt.marine > 0) & (cnt.land > 0) & (cnt.marine == cnt.land)].index
keep = pd.Series(np.where(cnt.marine >= cnt.land, 'marine', 'land'), index=cnt.index)
d = d[~d.species.isin(tie)]; d = d[d.domain.values == keep.reindex(d.species).values]
d['T_sd_total'] = np.where(d.domain == 'marine', d.sst_sd_total, d.tair_sd_total)
d = d[np.isfinite(d.T_sd_total)]
d = d[d.groupby('species').accession.transform('size') >= 2].copy()
main = pd.read_csv('data_small/analysis_set.tsv.gz', sep='\t', usecols=['accession', 'pilot'])
d['in_main'] = d.accession.isin(main.accession).astype(int)
d = d.merge(main.rename(columns={'pilot': 'main_pilot'}), on='accession', how='left')
d['is_mag'] = d.is_mag.astype(bool).astype(int)
d.to_csv(f'{Q}/selected_genomes_hq_mq.tsv.gz', sep='\t', index=False, compression='gzip')
new = d[~d.accession.isin(hq.accession)]
new[['accession']].assign(ncbi_assembly_name='', ncbi_translation_table=np.where(new.gtdb_taxonomy.str.contains('c__Gracilibacteria'), 25, 11)
    ).to_csv(f'{Q}/to_annotate.tsv', sep='\t', index=False)
mixed = d.groupby('species').is_mag.agg(lambda v: v.min() == 0 and v.max() == 1)
print('HQ+MQ genomes', len(d), 'species', d.species.nunique(), 'of which MQ', (d.quality == 'MQ').sum(),
      '| to annotate (not in HQ run)', len(new), '| mixed species', mixed.sum(), '(HQ-only: 74)')
print('HQ genomes no longer in set (majority-domain flips):', (~hq.accession.isin(d.accession)).sum())
