#!/usr/bin/env python3
"""WITHDRAWN (18 Nov 2026). The premise below is wrong: the real-data SE of 0.0022 came from one unusually
favourable random grouping (rank 4 of 200, scripts/30_grouping_sensitivity.R). The typical real-data SE is 0.0034-0.0035,
equal to the simulated slope SD, so no rescaling is needed and the 'real-data noise level' power figures in
results/tost_calibration_noise_scale.tsv do not apply. The invariance argument itself (linearity) remains correct.

Added 18 Nov 2026: relate the calibration of procedure P to the noise level of the real data.
The simulated data were noisier than the real data: the simulation SD of the slope was 0.0035 (mean jackknife SE
0.0037), whereas the jackknife SE of procedure P on the real data is 0.0022. The data-generating model added a site
effect on top of a phylogenetic residual variance that had been estimated with site variation already included,
so it is deliberately noisier than the real data.
Because the tree is fixed, the GLS slope is linear in the response and each jackknife replicate is linear in it too,
so for y = beta * P0 + k * e:  estimate - beta = k * d(e)  and  SE_jackknife = k * s(e).
At a boundary (beta = -0.015) P declares equivalence iff d/s > t_{0.95,19}, which does not involve k: the Type I
rate at the boundaries is exactly invariant to the overall noise scale. Power at interior slopes does depend on k.
This script pools the 14,000 simulated error draws (d, s), which do not depend on the true slope, and evaluates
P(declare equivalence) at the simulated noise level and at the real-data noise level (k chosen so that the mean
jackknife SE equals the real-data SE).
usage (repo root): python3 scripts/29_calibration_noise_scale.py -> results/tost_calibration_noise_scale.tsv"""
import numpy as np, pandas as pd
from scipy.stats import t as T, beta as B
s = pd.read_csv('results/tost_calibration_sims.tsv', sep='\t', header=None, skiprows=1)
s.columns = ['true_beta', 'seed', 'estimate', 'se', 'lo90', 'hi90', 'p_tost', 'equiv']
d, se = (s.estimate - s.true_beta).to_numpy(), s.se.to_numpy()
real_se = pd.read_csv('results/tost_calibration_real.tsv', sep='\t').set_index('set').at['H1_direct', 'se']
q, M = T.ppf(0.95, 19), 0.015
rows = []
for lab, k in (('simulated noise level', 1.0), ('real-data noise level', real_se / se.mean())):
    for b in (-0.015, -0.010, -0.005, 0.0, 0.005, 0.010, 0.015):
        est, sek = b + k * d, k * se
        x = int(((est - q * sek > -M) & (est + q * sek < M)).sum()); n = len(d)
        rows.append(dict(noise_level=lab, scale_k=k, mean_jackknife_se=k * se.mean(), sd_slope=k * d.std(ddof=1),
                         true_beta=b, nsim=n, P_declare_equivalence=x / n,
                         mc_lo95=B.ppf(0.025, x, n - x + 1) if x else 0.0, mc_hi95=B.ppf(0.975, x + 1, n - x) if x < n else 1.0))
o = pd.DataFrame(rows); o.to_csv('results/tost_calibration_noise_scale.tsv', sep='\t', index=False, float_format='%.5g')
print(o.round(4).to_string())
