#!/usr/bin/env python3
"""step 3: merge selected genomes + module counts + HE counts into the analysis table.
DORM_z uses the MAIN-ANALYSIS z parameters (dormz_params_main.tsv: mean/SD of log1p components over
the 7,867 main genomes, exactly as scripts/10_build_table.py computes them; verified to reproduce the
committed DORM_z to 5e-6). P0 = (log T_sd_total - main mean_domain) / main SD_domain, with the
confirmatory-set (pilot == 0, post-QC) parameters (p0_params_main.tsv), i.e. the same scale as 09_models.R.
QC as pre-registered: n_proteins >= 500 and n_HE >= 20 (NA counts as failure, unlike 10_build_table.py).
usage (repo root): python3 replication_within_species/03_build_table.py MODULES.tsv HE.tsv SELECTED.tsv.gz OUT.tsv.gz"""
import sys, numpy as np, pandas as pd
mods, he, sel, out = sys.argv[1:5]
R = 'replication_within_species/'
m = pd.read_csv(mods, sep='\t'); h = pd.read_csv(he, sep='\t'); s = pd.read_csv(sel, sep='\t', low_memory=False)
d = s.merge(m, on='accession', how='inner').merge(h[['accession', 'n_HE']], on='accession', how='left')
par = pd.read_csv(R + 'dormz_params_main.tsv', sep='\t').set_index('component')
comp = ['spo_strict', 'n_rpf', 'n_hib', 'n_sr']
d['DORM_z'] = np.mean([(np.log1p(d[c]) - par.at[c, 'mean']) / par.at[c, 'sd'] for c in comp], axis=0)
p0 = pd.read_csv(R + 'p0_params_main.tsv', sep='\t').set_index('domain')
d['log_T_sd_total'] = np.log(d.T_sd_total)
d['P0'] = (d.log_T_sd_total - d.domain.map(p0['mean'])) / d.domain.map(p0['std'])
d['T_mean'] = np.where(d.domain == 'marine', d.sst_mean, d.tair_mean)
d['hib_extra'] = (d.n_hib - d.hib_systems).clip(lower=0)
d['lgs'] = np.log10(d.genome_size)
d['is_sag'] = d.ncbi_genome_category.astype(str).str.contains('single cell').astype(int)
d['site'] = d.domain + ' ' + (np.round(d.lat * 2) / 2).astype(str) + ' ' + (np.round(d.lon * 2) / 2).astype(str)
d['qc_pass'] = ((d.n_proteins >= 500) & (d.n_HE >= 20)).astype(int)
print('genomes', len(d), 'QC fail', (d.qc_pass == 0).sum())
d = d[d.qc_pass == 1].copy()
n = d.groupby('species').accession.transform('size'); d = d[n >= 2].copy()
d['n_sp'] = d.groupby('species').accession.transform('size')
d['n_sites_sp'] = d.groupby('species').site.transform('nunique')
d['mixed_species'] = d.groupby('species').is_mag.transform(lambda v: int(v.min() == 0 and v.max() == 1))
keep = ['accession', 'species', 'phylum', 'gtdb_taxonomy', 'domain', 'habitat', 'is_mag', 'is_sag', 'in_main', 'main_pilot',
        'lat', 'lon', 'site', 'checkm2_completeness', 'checkm2_contamination', 'genome_size', 'lgs', 'n_proteins', 'n_HE',
        'T_sd_total', 'log_T_sd_total', 'P0', 'T_mean', 'DORM_z', 'hib_extra', 'n_sp', 'n_sites_sp', 'mixed_species'] + \
       [c for c in m.columns if c not in ('accession',)]
d[keep].to_csv(out, sep='\t', index=False, compression='gzip', float_format='%.6g')
print('analysed genomes', len(d), 'species', d.species.nunique(), 'mixed species', d[d.mixed_species == 1].species.nunique(),
      'species with >1 site', d[d.n_sites_sp > 1].species.nunique())
