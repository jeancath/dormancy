#!/usr/bin/env python3
"""Merge analysis set, module scores, growth-rate and chlorophyll layers into one table.
usage: 10_build_table.py MODULES.tsv GRODON.tsv OUT.tsv.gz
Derived variables:
  log_T_sd_total/seasonal/resid, T_ac1_resid, T_colwell_P/C/M : domain-matched temperature layer
  (SST for marine, air temperature for land)
  DORM_z : primary composite = mean of the four module z-scores of log1p(spo_strict), log1p(n_rpf),
           log1p(n_hib), log1p(n_sr) (z over the full analysis set; equal module weights)
  DORM_sum : spo_strict + n_rpf + n_hib + n_sr (sensitivity)
QC exclusions (pre-specified): n_proteins < 500 or < 20 ribosomal-protein HE genes (gRodon n_HE).
"""
import sys, numpy as np, pandas as pd

mods = pd.read_csv(sys.argv[1], sep='\t')
gro = pd.read_csv(sys.argv[2], sep='\t')
a = pd.read_csv('data_small/analysis_set.tsv.gz', sep='\t', low_memory=False)
chl = pd.read_csv('data_small/chl_variability.tsv.gz', sep='\t')
a = a.copy()
d = a.merge(mods, on='accession', how='inner').merge(gro[['accession', 'd_hours', 'n_HE']], on='accession', how='left')
d = d.merge(chl, on='accession', how='left')
nbm = pd.read_csv('data_small/niche_breadth_map.tsv.gz', sep='\t')[['accession', 'hao_niche_breadth', 'snb_genus']]
d = d.merge(nbm, on='accession', how='left')
d['log_hao_nb'] = np.log(d.hao_niche_breadth)
tx = d.gtdb_taxonomy.astype(str)
d['streamlined'] = tx.str.contains('o__Pelagibacterales|g__Prochlorococcus|o__Actinomarinales|o__SAR86', regex=True).astype(int)
d['is_actino'] = d.phylum.astype(str).str.startswith('Actinomycetota').astype(int)
mar = d.domain == 'marine'
for k in ('sd_total', 'sd_seasonal', 'sd_resid'):
    d[f'log_T_{k}'] = np.log(np.where(mar, d[f'sst_{k}'], d[f'tair_{k}']))
d['T_ac1_resid'] = np.where(mar, d.sst_ac1_resid, d.tair_ac1_resid)
d['T_frac_seasonal'] = np.where(mar, d.sst_frac_seasonal, d.tair_frac_seasonal)  # exploratory (Marshall & Burgess 2015 parameterisation)
for k in ('P', 'C', 'M'):
    d[f'T_colwell_{k}'] = np.where(mar, d[f'sst_colwell_{k}'], d[f'tair_colwell_{k}'])
d['T_mean'] = np.where(mar, d.sst_mean, d.tair_mean)
d['log_T_sd_total_spmed'] = np.log(np.where(mar, d.sst_sd_total_spmed, d.tair_sd_total_spmed))
d['log_prec_cv_total'] = np.log(d.prec_cv_total)
d['log_prec_cv_resid'] = np.log(d.prec_cv_resid)
d['log_chl_sd_total'] = np.log(d.chl_sd_total)
d['log_chl_sd_resid'] = np.log(d.chl_sd_resid)
d['log_chl_sd_seasonal'] = np.log(d.chl_sd_seasonal)
z = lambda v: (v - v.mean()) / v.std()
comp = ['spo_strict', 'n_rpf', 'n_hib', 'n_sr']
d['DORM_z'] = np.mean([z(np.log1p(d[c])) for c in comp], axis=0)
d['DORM_sum'] = d[comp].sum(axis=1)
d['log_d_hours'] = np.log(pd.to_numeric(d.d_hours, errors='coerce'))
d['log_ssu'] = np.log(d.ssu_count.where(d.ssu_count > 0))
d['tp_per_kcds'] = 1000 * d.n_transposase / d.n_proteins
d['log_mobile'] = np.log1p(d.n_transposase + d.n_integrase)
d['is_bacillota'] = d.phylum.str.startswith('Bacillota').astype(int)
d['is_cpr'] = (d.phylum == 'Patescibacteriota').astype(int)
d['complete_asm'] = d.ncbi_assembly_level.astype(str).str.contains('Complete').astype(int)
d['rpf_any'] = (d.n_rpf > 0).astype(int)
d['balon_any'] = (d.n_balon > 0).astype(int)
d['rmf_any'] = (d.n_rmf > 0).astype(int)
d['hpf_long_any'] = (d.hpf_long > 0).astype(int)
d['sas_any'] = (d.n_sas > 0).astype(int)
# exploratory additions (see results/deviations.md)
d['sporulator27'] = ((d.spo0A > 0) & (d.spo_strict >= 27)).astype(int)
d['abs_lat'] = d.lat.abs()
d['is_soil'] = (d.habitat == 'soil').astype(int)
d['is_sed'] = d.is_sediment.astype(str).isin(['True', '1']).astype(int)
# follow-up exploratory columns (results/deviations.md item 10)
_dec = d.ncbi_lat_lon.astype(str).str.findall(r'\d+\.(\d+)').apply(lambda v: min(len(x) for x in v) if len(v) else 0)
d['coord_decimals'] = _dec
d['coord_precise'] = (_dec >= 3).astype(int)
d['hib_extra'] = (d.n_hib - d.hib_systems).clip(lower=0)
d['rsh_extra'] = (d.n_rsh_long - 1).clip(lower=0)
d['is_isolate'] = 1 - d.is_mag.astype(str).isin(['True', '1']).astype(int)
d['log_n50'] = np.log10(d.n50_contigs)
d['log_contigs'] = np.log10(d.contig_count)
for c in ('Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota', 'Cyanobacteriota'):
    d['is_' + c[:6].lower()] = d.phylum.str.startswith(c).astype(int)
# Pre-registered QC: >= 500 proteins and >= 20 ribosomal-protein HE genes. Genomes whose gRodon run failed
# ('fewer than 10 HE genes', n_HE = NA) are excluded too (audit M7). The committed data_small/final table and all
# reported models were built before this fix and retain 6 such genomes; the corrected rerun (audit/qc_fix_rerun.tsv,
# n = 7,321) changes no decision.
qc = (d.n_proteins >= 500) & (pd.to_numeric(d.n_HE, errors='coerce') >= 20)
print('QC excluded:', (~qc).sum())
d = d[qc]
keep = [c for c in d.columns if not c.startswith(('ssu_gg', 'ssu_silva'))]
d = d.copy()
d[keep].to_csv(sys.argv[3], sep='\t', index=False, compression='gzip', float_format='%.6g')
print(d.groupby(['domain', 'pilot']).size())
