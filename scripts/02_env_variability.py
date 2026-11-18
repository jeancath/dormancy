#!/usr/bin/env python3
"""Phase 2b: organism-independent physical variability indices at isolation sites.

Grids (NOAA PSL), WMO normal period 1991-01 .. 2020-12 (360 months):
  * SST   : NOAA OISST v2.1 high-res monthly, 0.25 deg   (marine genomes)
  * T_air : GHCN_CAMS monthly 2 m air temperature, 0.5 deg (land / freshwater genomes)
  * Precip: NOAA PREC/L monthly precipitation, 0.5 deg     (land / freshwater genomes)

Per grid cell and variable:
  mean          long-term mean
  sd_total      SD of all 360 monthly values (seasonal + non-seasonal variability)
  sd_seasonal   SD of the 12-month climatology (predictable, seasonal amplitude)
  range_seas    max - min of the climatology
  sd_resid      SD of de-seasonalised, linearly de-trended anomalies (unpredictable part)
  ac1_resid     lag-1 autocorrelation of anomalies (reddening / persistence)
  frac_seasonal sd_seasonal^2 / sd_total^2 (share of variance that is seasonal)
  sd_interann   SD of the 30 annual means
Point-level only: Colwell (1974) predictability P, constancy C, contingency M on monthly values.
For precipitation, variability is expressed as coefficients of variation (divide by mean).

Points: each genome's lat/lon. If the nearest cell is masked (e.g. coastal isolate on an OISST
land cell) the nearest valid cell within `max_cells` is used and the distance recorded.
"""
import sys, numpy as np, pandas as pd, xarray as xr

ENV = '/home/user/data/env'
T0, T1 = '1991-01-01', '2020-12-31'


def load(fname, var):
    d = xr.open_dataset(f'{ENV}/{fname}')[var].sel(time=slice(T0, T1))
    assert d.sizes['time'] == 360, d.sizes
    if float(d.lat[0]) > float(d.lat[-1]):
        d = d.isel(lat=slice(None, None, -1))
    return d.load().astype('float32')


def grid_metrics(da):
    x = da.values  # time, lat, lon
    nt = x.shape[0]
    clim = np.stack([np.nanmean(x[m::12], axis=0) for m in range(12)])
    anom = x - np.tile(clim, (nt // 12, 1, 1))
    t = np.arange(nt, dtype='float32') - (nt - 1) / 2
    slope = np.nansum(anom * t[:, None, None], axis=0) / np.sum(t ** 2)
    anom = anom - slope[None] * t[:, None, None]
    out = {}
    out['mean'] = np.nanmean(x, axis=0)
    out['sd_total'] = np.nanstd(x, axis=0)
    out['sd_seasonal'] = np.nanstd(clim, axis=0)
    out['range_seas'] = np.nanmax(clim, axis=0) - np.nanmin(clim, axis=0)
    out['sd_resid'] = np.nanstd(anom, axis=0)
    a0, a1 = anom[:-1], anom[1:]
    num = np.nanmean(a0 * a1, axis=0)
    out['ac1_resid'] = num / np.nanvar(anom, axis=0)
    out['frac_seasonal'] = out['sd_seasonal'] ** 2 / out['sd_total'] ** 2
    ann = x.reshape(nt // 12, 12, *x.shape[1:]).mean(axis=1)
    out['sd_interann'] = np.nanstd(ann, axis=0)
    valid = np.isfinite(x).sum(axis=0) >= 300
    for k in out:
        out[k] = np.where(valid, out[k], np.nan).astype('float32')
    return out, valid


def colwell(series, bins):
    """Colwell 1974 predictability for a monthly series (len multiple of 12)."""
    s = np.asarray(series)
    if np.isnan(s).any() or s.size % 12:
        return np.nan, np.nan, np.nan
    st = np.digitize(s, bins)
    nst = len(bins) + 1
    N = np.zeros((nst, 12))
    for i, v in enumerate(st):
        N[v, i % 12] += 1
    Z = N.sum()
    X = N.sum(axis=0); Y = N.sum(axis=1)
    def H(p):
        p = p[p > 0] / Z
        return -(p * np.log(p)).sum()
    HX, HY, HXY = H(X), H(Y), H(N.ravel())
    s_ = np.log(nst)
    P = 1 - (HXY - HX) / s_
    C = 1 - HY / s_
    return P, C, P - C


def nearest_valid(valid, lat_axis, lon_axis, lat, lon, max_cells):
    lon360 = lon % 360 if lon_axis.min() >= 0 else lon
    i = int(np.abs(lat_axis - lat).argmin()); j = int(np.abs(lon_axis - lon360).argmin())
    if valid[i, j]:
        return i, j, 0.0
    best = None
    nlat, nlon = valid.shape
    for r in range(1, max_cells + 1):
        for di in range(-r, r + 1):
            for dj in range(-r, r + 1):
                if max(abs(di), abs(dj)) != r:
                    continue
                ii = i + di; jj = (j + dj) % nlon
                if 0 <= ii < nlat and valid[ii, jj]:
                    dlat = np.radians(lat_axis[ii] - lat)
                    dlon = np.radians(((lon_axis[jj] - lon360 + 180) % 360) - 180)
                    a = np.sin(dlat / 2) ** 2 + np.cos(np.radians(lat)) * np.cos(np.radians(lat_axis[ii])) * np.sin(dlon / 2) ** 2
                    d = 6371 * 2 * np.arcsin(np.sqrt(a))
                    if best is None or d < best[2]:
                        best = (ii, jj, d)
        if best is not None:
            return best
    return None


def extract(name, da, pts, max_cells, colwell_bins, cv=False):
    m, valid = grid_metrics(da)
    lat_axis = da.lat.values; lon_axis = da.lon.values
    rows = []
    x = da.values
    for idx, (lat, lon) in pts[['lat', 'lon']].iterrows():
        r = nearest_valid(valid, lat_axis, lon_axis, lat, lon, max_cells)
        rec = {'accession': pts.at[idx, 'accession']}
        if r is None:
            rows.append(rec); continue
        i, j, d = r
        rec[f'{name}_dist_km'] = d
        for k, v in m.items():
            rec[f'{name}_{k}'] = float(v[i, j])
        if cv:
            mu = rec[f'{name}_mean']
            for k in ('sd_total', 'sd_seasonal', 'sd_resid', 'sd_interann'):
                rec[f'{name}_cv_{k[3:]}'] = rec[f'{name}_{k}'] / mu if mu > 0 else np.nan
        P, C, M = colwell(x[:, i, j], colwell_bins(x[:, i, j]))
        rec[f'{name}_colwell_P'] = P; rec[f'{name}_colwell_C'] = C; rec[f'{name}_colwell_M'] = M
        rows.append(rec)
    return pd.DataFrame(rows)


def main():
    pts = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else 'data_small/habitat_classified.tsv.gz', sep='\t',
                      usecols=['accession', 'lat', 'lon'])
    out = sys.argv[2] if len(sys.argv) > 2 else 'data_small/env_variability.tsv.gz'
    # unique coordinates to save time
    u = pts[['lat', 'lon']].drop_duplicates().reset_index(drop=True)
    u['accession'] = u.index.astype(str)
    print('unique sites', len(u), flush=True)
    temp_bins = lambda s: np.arange(-40, 45, 2.0)  # 2 degC classes, fixed across sites
    prec_bins = lambda s: np.array([0.1, 0.5, 1, 2, 4, 8, 16, 32])  # mm/day, log-like classes
    res = []
    sst = load('sst.mon.mean.nc', 'sst')
    res.append(extract('sst', sst, u, 4, temp_bins)); del sst
    print('sst done', flush=True)
    air = load('air.mon.mean.nc', 'air')
    if float(air.mean()) > 100:  # Kelvin -> Celsius
        air = air - 273.15
    res.append(extract('tair', air, u, 2, temp_bins)); del air
    print('tair done', flush=True)
    pr = load('precip.mon.mean.0.5x0.5.nc', 'precip')
    res.append(extract('prec', pr, u, 2, prec_bins, cv=True)); del pr
    print('precip done', flush=True)
    site = u.copy()
    for r in res:
        site = site.merge(r, on='accession', how='left')
    site = site.drop(columns='accession')
    full = pts.merge(site, on=['lat', 'lon'], how='left')
    full.to_csv(out, sep='\t', index=False, compression='gzip', float_format='%.5g')
    print(full.describe().T.to_string())


if __name__ == '__main__':
    main()
