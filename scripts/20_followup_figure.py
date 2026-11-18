#!/usr/bin/env python3
"""Fig 6 (exploratory follow-ups): which hibernation / stringent-response genes carry the MAG-isolate reversal,
and Balon presence vs mean temperature in marine isolates.
usage: 20_followup_figure.py TABLE.tsv.gz FOLLOWUP_MODELS.tsv OUTDIR"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

tab, mod, outd = sys.argv[1:4]
m = pd.read_csv(mod, sep='\t').set_index('id')
d = pd.read_csv(tab, sep='\t', low_memory=False)
d = d[d.pilot == 0]
CM = {'MAG': '#7b4fc9', 'iso': '#1f9e74'}
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight'})
fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.0), gridspec_kw={'width_ratios': [1.25, 1]})
ax = axs[0]
rows = [('hib_systems', 'hibernation systems present'), ('hib_extra', 'extra hibernation copies'),
        ('hpfshort', 'short HPF / RaiA'), ('hpflong', 'long HPF'), ('rmf', 'RMF'), ('balon', 'Balon')]
y = np.arange(len(rows))[::-1]
ax.axvline(0, color=MUTED, lw=0.8)
for yi, (k, lab) in zip(y, rows):
    for g, off in (('MAG', 0.15), ('iso', -0.15)):
        i = f'F_{g}_{k}'
        if i in m.index:
            e, se = m.at[i, 'estimate'], m.at[i, 'se']
            ax.plot([e - 1.96 * se, e + 1.96 * se], [yi + off] * 2, color=CM[g], lw=2, solid_capstyle='round')
            ax.plot(e, yi + off, 'o', color=CM[g], ms=4.5, mec='white', mew=0.8,
                    label=('MAGs' if g == 'MAG' else 'isolates') if yi == y[0] else None)
ax.set_yticks(y); ax.set_yticklabels([r[1] for r in rows])
ax.set_xlabel('slope per SD of P0 (log1p count, 95% CI)')
ax.legend(frameon=False, fontsize=7, loc='lower right')
ax.set_title('a  The reversal sits in hibernation gene copies', loc='left', color=INK, fontsize=9)
ax = axs[1]
s = d[(d.domain == 'marine') & (d.is_mag.astype(str).isin(['False', '0']))].copy()
s['bin'] = pd.qcut(s.sst_mean, 8, duplicates='drop')
g = s.groupby('bin', observed=True).agg(x=('sst_mean', 'mean'), y=('balon_any', 'mean'), n=('balon_any', 'size'))
se = np.sqrt(g.y * (1 - g.y) / g.n)
ax.fill_between(g.x, 100 * (g.y - 1.96 * se), 100 * (g.y + 1.96 * se), color='#2a78d6', alpha=0.15, lw=0)
ax.plot(g.x, 100 * g.y, '-o', color='#2a78d6', lw=2, ms=3.5)
ax.set_xlabel('mean SST at isolation site (°C)'); ax.set_ylabel('genomes with Balon (%)')
ax.set_title(f'b  Balon in marine isolates (n = {len(s):,})', loc='left', color=INK, fontsize=9)
fig.tight_layout(); fig.savefig(f'{outd}/fig6_followups.png'); plt.close(fig)
print('written', f'{outd}/fig6_followups.png')
