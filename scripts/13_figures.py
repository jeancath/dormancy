#!/usr/bin/env python3
"""Figures for the confirmatory analysis.
usage: 13_figures.py TABLE.tsv.gz MODELS.tsv OUTDIR [EXPL_MODELS.tsv SITE_CAPPED_BETAS.tsv]
Palette: reference categorical slots 1-2 (blue = marine, orange = land), validated all-pairs.
"""
import sys, os, numpy as np, pandas as pd, xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

tab, mod, outd = sys.argv[1:4]
os.makedirs(outd, exist_ok=True)
d = pd.read_csv(tab, sep='\t', low_memory=False)
d = d[d.pilot == 0]
m = pd.read_csv(mod, sep='\t')
C = {'marine': '#2a78d6', 'land': '#eb6834'}
CM = {'MAG': '#7b4fc9', 'iso': '#1f9e74'}  # genome type; distinct from the domain hues
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight'})

# ---- Fig 1: isolation sites over a land mask, plus predictor distributions
mask = xr.open_dataset('/home/user/data/env/lsmask.oisst.nc').lsmask.isel(time=0)
fig = plt.figure(figsize=(7.2, 4.6))
ax = fig.add_axes([0.0, 0.32, 1.0, 0.68])
lm = np.roll(mask.values, 720, axis=1)
ax.imshow(np.where(lm == 0, 1, np.nan), extent=[-180, 180, -90, 90], origin='lower', cmap='Greys', vmin=0, vmax=4,
          interpolation='nearest')
for dom in ('land', 'marine'):
    s = d[d.domain == dom]
    ax.scatter(s.lon, s.lat, s=5, c=C[dom], lw=0, alpha=0.75, label=f'{dom} (n = {len(s):,})')
ax.set_xlim(-180, 180); ax.set_ylim(-80, 85); ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
ax.legend(loc='lower left', frameon=False, markerscale=3)
ax.set_title('Isolation sites of the confirmatory species (one genome per GTDB species)', color=INK, loc='left')
for i, (dom, var, lab) in enumerate((('marine', 'sst_sd_total', 'SST monthly SD, °C'), ('land', 'tair_sd_total', 'Air temperature monthly SD, °C'))):
    a2 = fig.add_axes([0.07 + 0.5 * i, 0.06, 0.4, 0.2])
    a2.hist(d[d.domain == dom][var], bins=40, color=C[dom], edgecolor='white', linewidth=0.5)
    a2.set_xlabel(lab); a2.set_ylabel('species')
fig.savefig(f'{outd}/fig1_sites.png'); plt.close(fig)

# ---- Fig 2: forest plot of pre-registered slopes
order = [('H1_DORM_P0', 'H1 DORM, direct effect (primary)'), ('H1_DORM_P0_total', 'H1 DORM, total effect'),
         ('H1a_spo_strict', 'H1a sporulation genes (Bacillota)'), ('H1b_rpf', 'H1b Rpf'),
         ('H1b_rpf_actino', 'H1b Rpf (Actinomycetota)'), ('H1c_hib', 'H1c hibernation'), ('H1d_sr', 'H1d stringent response'),
         ('H6_seasonal', 'H6 seasonal SD'), ('H6_resid', 'H6 residual SD'), ('H6_ac1', 'H6 anomaly lag-1 AC'),
         ('H6_colwellP', 'H6 Colwell predictability'), ('H6_colwellM', 'H6 Colwell contingency'),
         ('H6_colwellC', 'H6 Colwell constancy'),
         ('S_isolates', 'sens: isolates only'), ('S_comp95', 'sens: completeness ≥95%'), ('S_spmed', 'sens: species-median env'),
         ('S_marine', 'sens: marine only'), ('S_land', 'sens: land only'), ('S_noBacCPR', 'sens: no Bacillota/CPR'),
         ('S_growth', 'sens: + gRodon doubling time'), ('S_mobile', 'sens: + mobile elements'), ('S_Tmean', 'sens: + mean temperature'),
         ('X_chl_total', 'secondary: chlorophyll SD (marine)'), ('X_prec_cv', 'secondary: precipitation CV (land)'),
         ('X_hao_nb', 'secondary: habitat niche breadth'), ('X_snb', 'secondary: social niche breadth (genus)')]
f = m.set_index('id')
rows = [(lab, f.at[i, 'estimate'], f.at[i, 'se'], f.at[i, 'n']) for i, lab in order if i in f.index and pd.notna(f.at[i, 'estimate'])]
fig, ax = plt.subplots(figsize=(6.4, 0.24 * len(rows) + 0.8))
y = np.arange(len(rows))[::-1]
ax.axvspan(-0.015, 0.015, color=GRID, zorder=0, label='±SESOI (0.015)')
ax.axvline(0, color=MUTED, lw=0.8)
for yi, (lab, e, se, n) in zip(y, rows):
    col = INK if lab.startswith('H1 ') else '#2a78d6'
    ax.plot([e - 1.96 * se, e + 1.96 * se], [yi, yi], color=col, lw=2, solid_capstyle='round')
    ax.plot(e, yi, 'o', color=col, ms=5, mec='white', mew=1)
ax.set_yticks(y); ax.set_yticklabels([f'{r[0]}  (n={int(r[3]):,})' for r in rows])
ax.set_xlabel('slope (response units per 1 SD of predictor), 95% CI')
ax.legend(frameon=False, loc='upper right')
ax.set_title('Phylogenetic regression slopes (confirmatory set)', loc='left', color=INK)
fig.savefig(f'{outd}/fig2_forest.png'); plt.close(fig)

# ---- Fig 3: binned DORM_z and module prevalence against P0, by domain
z = lambda v: (v - v.mean()) / v.std()
d['P0'] = d.groupby('domain').log_T_sd_total.transform(z)
fig, axs = plt.subplots(1, 4, figsize=(7.4, 2.3), sharex=True)
for ax, (col, lab, fn) in zip(axs, (('DORM_z', 'DORM_z (mean)', 'mean'), ('sporulator', 'sporulators (%)', 'pct'),
                                     ('rpf_any', 'Rpf present (%)', 'pct'), ('hib_systems', 'hibernation systems (mean)', 'mean'))):
    for dom in ('land', 'marine'):
        s = d[d.domain == dom].copy()
        s['bin'] = pd.qcut(s.P0, 10, labels=False, duplicates='drop')
        g = s.groupby('bin').agg(x=('P0', 'mean'), y=(col, 'mean'), se=(col, 'sem'))
        k = 100 if fn == 'pct' else 1
        ax.fill_between(g.x, k * (g.y - 1.96 * g.se), k * (g.y + 1.96 * g.se), color=C[dom], alpha=0.15, lw=0)
        ax.plot(g.x, k * g.y, '-o', color=C[dom], lw=2, ms=3, label=dom)
    ax.set_title(lab, fontsize=8, color=INK, loc='left'); ax.set_xlabel('P0 (z, within domain)')
axs[0].legend(frameon=False, fontsize=7)
fig.suptitle('Raw (not phylogenetically corrected) deciles of temperature variability', x=0.01, ha='left', fontsize=9, color=INK)
fig.tight_layout(); fig.savefig(f'{outd}/fig3_deciles.png'); plt.close(fig)

# ---- Fig 4: TA loci vs genome size by domain, isolates and MAGs separately (H3)
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
for ax, (mag, lab) in zip(axs, ((False, 'isolate genomes'), (True, 'MAGs'))):
    for dom in ('land', 'marine'):
        s = d[(d.domain == dom) & (d.is_mag.astype(str).isin(['True', '1']) == mag)]
        ax.scatter(s.genome_size / 1e6, s.ta_pairs, s=4, c=C[dom], alpha=0.3, lw=0)
        s2 = s.assign(b=pd.cut(s.genome_size / 1e6, np.arange(0.5, 12.5, 0.5))).groupby('b', observed=True).agg(
            x=('genome_size', 'median'), y=('ta_pairs', 'mean'), n=('ta_pairs', 'size'))
        s2 = s2[s2.n >= 15]
        ax.plot(s2.x / 1e6, s2.y, '-', color=C[dom], lw=2, label=f'{dom} (n = {len(s):,})')
    ax.set_yscale('symlog', linthresh=1); ax.set_xlabel('genome size (Mb)')
    ax.set_title(lab, loc='left', color=INK, fontsize=9); ax.legend(frameon=False, fontsize=7, loc='upper left')
axs[0].set_ylabel('paired type II TA loci')
fig.suptitle('TA loci scale with genome size; marine deficit is largest in MAGs', x=0.01, ha='left', fontsize=9, color=INK)
fig.tight_layout(); fig.savefig(f'{outd}/fig4_ta_genome_size.png'); plt.close(fig)

# ---- Fig 5 (exploratory): MAG vs isolate slopes and site-capped resampling
if len(sys.argv) > 5:
    ex = pd.read_csv(sys.argv[4], sep='\t').set_index('id')
    cb = pd.read_csv(sys.argv[5], sep='\t')
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.8), gridspec_kw={'width_ratios': [1.1, 1]})
    ax = axs[0]
    ids = [('H1_DORM_P0', 'all species (H1)', m.set_index('id')), ('E_MAGs_only', 'MAGs', ex), ('S_isolates', 'isolates', m.set_index('id')),
           ('E_MAGs_marine', 'MAGs, marine', ex), ('E_MAGs_land', 'MAGs, land', ex), ('E_iso_marine', 'isolates, marine', ex),
           ('E_iso_land', 'isolates, land', ex)]
    ids = [t for t in ids if t[0] in t[2].index]
    yy = np.arange(len(ids))[::-1]
    ax.axvspan(-0.015, 0.015, color=GRID, zorder=0); ax.axvline(0, color=MUTED, lw=0.8)
    for yi, (i, lab, src) in zip(yy, ids):
        e, se = src.at[i, 'estimate'], src.at[i, 'se']
        col = INK if i == 'H1_DORM_P0' else (CM['MAG'] if 'MAG' in lab else CM['iso'])
        ax.plot([e - 1.96 * se, e + 1.96 * se], [yi, yi], color=col, lw=2, solid_capstyle='round')
        ax.plot(e, yi, 'o', color=col, ms=5, mec='white', mew=1)
    ax.set_yticks(yy); ax.set_yticklabels([t[1] for t in ids]); ax.set_xlabel('DORM_z slope per SD of P0 (phylolm 95% CI)')
    ax.set_title('a  The pooled effect comes from MAGs', loc='left', color=INK, fontsize=9)
    ax = axs[1]
    for nm, col, lab in (('H1_direct', INK, 'all species'), ('MAGs', CM['MAG'], 'MAGs'), ('isolates', CM['iso'], 'isolates')):
        v = cb[cb.set == nm].beta
        if len(v):
            ax.hist(v, bins=30, histtype='step', color=col, lw=2, label=lab)
    ax.axvline(0, color=MUTED, lw=0.8); ax.axvspan(-0.015, 0.015, color=GRID, zorder=0)
    ax.set_xlabel('DORM_z slope, one random species per 0.5° site'); ax.set_ylabel('resamples')
    ax.legend(frameon=False, fontsize=7); ax.set_title('b  Site-capped resampling', loc='left', color=INK, fontsize=9)
    fig.tight_layout(); fig.savefig(f'{outd}/fig5_mag_isolate_sites.png'); plt.close(fig)
print('figures written to', outd)
