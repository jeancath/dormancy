#!/usr/bin/env python3
"""Phase 2c/3: assemble the analysis set (one genome per GTDB species per domain) and draw the
held-out pilot.

Domains (surface, free-living environments where a surface climate series is meaningful):
  marine : habitat == marine, not deep-water keywords, nearest valid OISST cell <= 50 km;
           plus 'water_unspecified' genomes whose coordinates fall on an OISST ocean cell and not on
           a GHCN_CAMS land cell (flag rescued=1).
  land   : habitat in {soil, freshwater} and the site falls on a valid GHCN_CAMS land cell.
Excluded by design: host/clinical, engineered, subsurface, hydrothermal, hypersaline, unknown, other.

Per species and domain the genome with the highest quality score
  Q = completeness - 5*contamination + 5*(isolate) + 2*(complete assembly)
is retained. A species present in both domains is kept in the domain where it has more genomes
(ties -> dropped from both; recorded).

Pilot: 400 species (seed 20261117) drawn from the extremes (lowest and highest 20%) of the primary
index (log SD of the 1991-2020 monthly temperature series: SST for marine, air T for land),
200 from each tail, proportional by domain. Pilot species are excluded from all confirmatory tests.
"""
import numpy as np, pandas as pd

hab = pd.read_csv('data_small/habitat_classified.tsv.gz', sep='\t', low_memory=False)
env = pd.read_csv('data_small/env_variability.tsv.gz', sep='\t')
d = hab.merge(env.drop(columns=['lat', 'lon']), on='accession', how='left')
extra = pd.read_csv('/home/user/data/gtdb/bac120_metadata.tsv.gz', sep='\t', low_memory=False,
                    usecols=['accession', 'ncbi_assembly_name', 'ncbi_translation_table', 'trna_count',
                             'lsu_23s_count', 'lsu_5s_count', 'n50_contigs', 'gtdb_type_designation_ncbi_taxa'])
d = d.merge(extra, on='accession', how='left')

marine = (d.habitat == 'marine') & (~d.is_deep) & (d.sst_dist_km <= 50)
rescued = (d.habitat == 'water_unspecified') & (d.sst_dist_km == 0) & (d.tair_dist_km > 0) & (~d.is_deep)
land = d.habitat.isin(['soil', 'freshwater']) & (d.tair_dist_km == 0)
d['domain'] = np.select([marine | rescued, land], ['marine', 'land'], default='')
d['rescued'] = rescued.astype(int)
d = d[d.domain != ''].copy()

iso = (~d.is_mag).astype(int)
comp = d.ncbi_assembly_level.astype(str).str.contains('Complete', case=False).astype(int)
d['Q'] = d.checkm2_completeness - 5 * d.checkm2_contamination + 5 * iso + 2 * comp

# species with genomes in both domains
cnt = d.groupby(['species', 'domain']).size().unstack(fill_value=0)
both = cnt[(cnt.get('marine', 0) > 0) & (cnt.get('land', 0) > 0)]
tie = both[both.marine == both.land].index
keep_dom = pd.Series(np.where(cnt.get('marine', 0) >= cnt.get('land', 0), 'marine', 'land'), index=cnt.index)
d = d[~d.species.isin(tie)]
d = d[d.domain.values == keep_dom.reindex(d.species).values]
d['n_genomes_species_domain'] = d.groupby('species').accession.transform('size')

# species-median environment (sensitivity analysis)
envcols = [c for c in d.columns if c.startswith(('sst_', 'tair_', 'prec_'))]
med = d.groupby('species')[envcols].median().add_suffix('_spmed')

best = d.sort_values('Q', ascending=False).groupby('species', as_index=False).head(1).copy()
best = best.merge(med, left_on='species', right_index=True, how='left')

# primary index
best['T_sd_total'] = np.where(best.domain == 'marine', best.sst_sd_total, best.tair_sd_total)
best = best[np.isfinite(best.T_sd_total)]
best['log_T_sd_total'] = np.log(best.T_sd_total)

rng = np.random.default_rng(20261117)
best['pilot'] = 0
for dom, n_each in (('marine', None), ('land', None)):
    sub = best[best.domain == dom]
    lo, hi = sub.log_T_sd_total.quantile([0.2, 0.8])
    share = len(sub) / len(best)
    k = int(round(200 * share))
    for tail in (sub[sub.log_T_sd_total <= lo], sub[sub.log_T_sd_total >= hi]):
        pick = rng.choice(tail.index.values, size=min(k, len(tail)), replace=False)
        best.loc[pick, 'pilot'] = 1

best.to_csv('data_small/analysis_set.tsv.gz', sep='\t', index=False, compression='gzip')
best[['accession', 'ncbi_assembly_name', 'ncbi_translation_table']].to_csv(
    'data_small/analysis_accessions.tsv', sep='\t', index=False)
best[best.pilot == 1][['accession', 'ncbi_assembly_name', 'ncbi_translation_table']].to_csv(
    'data_small/pilot_accessions.tsv', sep='\t', index=False)
print('dropped tie species:', len(tie))
print(best.groupby(['domain', 'pilot']).size())
print(best.groupby('domain').phylum.value_counts().groupby(level=0).head(8))
print('isolates:', (~best.is_mag).sum(), 'MAG/SAG:', best.is_mag.sum(), 'rescued:', best.rescued.sum())
