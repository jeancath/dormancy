#!/usr/bin/env python3
"""Figure script: (a) within-species P0 slopes in MAG and isolate genomes (species fixed effects,
site-clustered 95 % CI) next to the between-species main-analysis slopes; (b) mixed species:
MAG − isolate difference in DORM_z (species means) against the P0 difference.
Light background, two hues: purple #7b4fc9 = MAG, green #1f9e74 = isolate.
usage (repo root): python3 replication_within_species/06_figure.py RESULTS.tsv MIXED_PAIRS.tsv OUT.png"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
res, pairs, out = sys.argv[1:4]
r = pd.read_csv(res, sep='\t'); p = pd.read_csv(pairs, sep='\t')
MAG, ISO, INK, MUTED = '#7b4fc9', '#1f9e74', '#0b0b0b', '#52514e'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.dpi': 300, 'savefig.bbox': 'tight'})
main = {('DORM_z', 'MAG'): (0.01807, 0.00339), ('DORM_z', 'iso'): (-0.00792, 0.00365),
        ('hib_extra', 'MAG'): (0.01320, 0.00313), ('hib_extra', 'iso'): (-0.01119, 0.00353),
        ('log_hpf_short', 'MAG'): (0.01342, 0.00348), ('log_hpf_short', 'iso'): (-0.00572, 0.00278),
        ('log_n_balon', 'MAG'): (0.00804, 0.00256), ('log_n_balon', 'iso'): (-0.00166, 0.00338),
        ('log_n_sr', 'MAG'): (0.00884, 0.00319), ('log_n_sr', 'iso'): (-0.00213, 0.00360)}
rows = [('DORM_z', 'DORM_z'), ('hib_extra', 'extra hibernation copies'), ('log_hpf_short', 'short HPF / RaiA'),
        ('log_n_balon', 'Balon'), ('log_n_sr', 'stringent response')]
fig, axs = plt.subplots(1, 2, figsize=(7.8, 3.3), gridspec_kw={'width_ratios': [1.35, 1]})
ax = axs[0]; ax.axvline(0, color=MUTED, lw=0.8); y0 = np.arange(len(rows))[::-1]
sub = r[r.model == 'b_speciesFE_clsite']
for yi, (k, lab) in zip(y0, rows):
    for g, col, off, term in (('MAG', MAG, 0.17, 'P0 (MAGs)'), ('iso', ISO, -0.17, 'P0 (isolates)')):
        x = sub[(sub.response == k) & (sub.term == term)]
        if len(x):
            e, s = x.estimate.iloc[0], x.se.iloc[0]
            ax.plot([e - 1.96 * s, e + 1.96 * s], [yi + off] * 2, color=col, lw=2.2, solid_capstyle='round')
            ax.plot(e, yi + off, 'o', color=col, ms=5, mec='white', mew=0.8,
                    label=(('MAG genomes' if g == 'MAG' else 'isolate genomes') + ', within species') if yi == y0[0] else None)
        if (k, g) in main:
            e, s = main[(k, g)]
            ax.plot(e, yi + off, 'D', color='white', mec=col, mew=1.2, ms=4.5,
                    label=('between species (main analysis)' if (g == 'MAG' and yi == y0[0]) else None))
ax.set_yticks(y0); ax.set_yticklabels([l for _, l in rows])
ax.set_xlabel('slope per SD of P0 (95 % CI, site-clustered)')
h, l = ax.get_legend_handles_labels()
ax.set_title('a  Within-species slopes by genome type', loc='left', fontsize=9, color=INK)
ax = axs[1]
ax.axhline(0, color=MUTED, lw=0.8); ax.axvline(0, color=MUTED, lw=0.8)
sz = 8 + 4 * np.sqrt(p.nM + p.nI)
ax.scatter(p.dP0, p.dDORM, s=sz, color=MAG, alpha=0.55, edgecolor='white', lw=0.5)
if len(p) > 3:
    b = np.polyfit(p.dP0, p.dDORM, 1); xx = np.linspace(p.dP0.min(), p.dP0.max(), 10)
    ax.plot(xx, np.polyval(b, xx), color=MAG, lw=1.5, label='unadjusted OLS fit')
    ax.legend(frameon=False, fontsize=7, loc='upper left')
ax.set_xlabel('ΔP0 (MAG sites − isolate sites, SD)')
ax.set_ylabel('Δ DORM_z (MAG − isolate)')
ax.set_title(f'b  Species with both genome types (n = {len(p)})', loc='left', fontsize=9, color=INK)
fig.tight_layout(rect=(0, 0.08, 1, 1))
fig.legend(h, l, frameon=False, fontsize=7, loc='lower center', ncol=3, bbox_to_anchor=(0.5, 0))
fig.savefig(out); print('written', out)
