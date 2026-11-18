#!/usr/bin/env python3
"""Marine resource-variability layer: ESA OC-CCI v6 monthly chlorophyll-a (1998-2022, ~0.25 deg).
Metrics on log10(chl) at each marine site (nearest cell with >= 60% valid months within 4 cells):
  chl_mean, chl_sd_total, chl_sd_seasonal, chl_sd_resid, chl_ac1_resid, chl_frac_valid, chl_colwell_P/C/M
Gaps (cloud, polar night) are left missing; climatology uses available months (>= 10 of 12 months
must have a climatological value).
usage: 02c_chl_variability.py analysis_set.tsv.gz OUT.tsv.gz
"""
import sys, glob, numpy as np, pandas as pd, xarray as xr
sys.path.insert(0, 'scripts')
from importlib import import_module
ev = import_module('02_env_variability')

files = sorted(glob.glob('/home/user/data/env/chl/chl_*.nc'))
ds = xr.concat([xr.open_dataset(f)['chlor_a'].load() for f in files], dim='time')
ds = ds.sortby('latitude')
x = np.log10(ds.values.astype('float32'))
lat = ds.latitude.values; lon = ds.longitude.values
nt = x.shape[0]
valid_frac = np.isfinite(x).mean(axis=0)
valid = valid_frac >= 0.6

a = pd.read_csv(sys.argv[1], sep='\t', low_memory=False)
a = a[a.domain == 'marine'][['accession', 'lat', 'lon']]
rows = []
bins = np.arange(-2.0, 1.6, 0.25)
for r in a.itertuples():
    res = ev.nearest_valid(valid, lat, lon, r.lat, r.lon, 4)
    rec = {'accession': r.accession}
    if res is None:
        rows.append(rec); continue
    i, j, d = res
    s = x[:, i, j]
    clim = np.array([np.nanmean(s[m::12]) if np.isfinite(s[m::12]).any() else np.nan for m in range(12)])
    rec.update(chl_dist_km=d, chl_frac_valid=float(np.isfinite(s).mean()), chl_mean=float(np.nanmean(s)),
               chl_sd_total=float(np.nanstd(s)))
    if np.isfinite(clim).sum() >= 10:
        an = s - np.tile(clim, nt // 12)
        rec['chl_sd_seasonal'] = float(np.nanstd(clim))
        rec['chl_sd_resid'] = float(np.nanstd(an))
        ok = np.isfinite(an[:-1]) & np.isfinite(an[1:])
        if ok.sum() > 24:
            rec['chl_ac1_resid'] = float(np.corrcoef(an[:-1][ok], an[1:][ok])[0, 1])
    if np.isfinite(s).all():
        P, C, M = ev.colwell(s, bins)
        rec.update(chl_colwell_P=P, chl_colwell_C=C, chl_colwell_M=M)
    rows.append(rec)
out = pd.DataFrame(rows)
out.to_csv(sys.argv[2], sep='\t', index=False, compression='gzip', float_format='%.5g')
print(out.describe().T.to_string())
