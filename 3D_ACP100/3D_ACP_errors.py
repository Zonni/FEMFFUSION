# -*- coding: utf-8 -*-
"""
Postprocessing – ACP-100 SMR neutron noise (sinusoidal Sigma_a perturbation)
Three ROM strategies: Dynamic, Static time-sampled, Parametric (XS-sampled)
"""

import sys
sys.path.append('../postprocess')
from utils import (get_td_power, get_td_time, parse_time_file,
                   power_at_time, spatial_L2_error_timeavg, power_errors)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

plt.close('all')

# =============================================================================
# Plot settings
# =============================================================================
params = {
    'backend': 'pdf',
    'font.size': 14,
    'axes.labelsize': 12,
    'legend.fontsize': 10,
    'xtick.labelsize': 13,
    'ytick.labelsize': 13,
    'text.usetex': False,
    'lines.linewidth': 1.4,
    'lines.markersize': 3,
    'lines.markeredgewidth': 1,
    'legend.numpoints': 1,
    'axes.formatter.useoffset': False,
    'figure.autolayout': True,
}
rcParams.update(params)

# =============================================================================
# File names and labels
# Adapt ROM filenames once you have them
# =============================================================================

labels = [
    'FOM',
    'Dynamic',
     'Static time-sampled',
     'Static xs-sampled',
]

out_files = [
    'ACP100_FE2_td_1cell_FOM.out',
    'ACP100_FE2_td_1cell_rom_dynamic.out',
     'ACP100_FE2_td_1cell_rom_static_time_sampled.out',
    'ACP100_FE2_td_1cell_rom_static_xs_sampled.out',
]

# =============================================================================
# Read output files
# =============================================================================
time       = []
powers     = []
local_times = []
local_pows  = []

for i, fname in enumerate(out_files):
    print(f'Reading {fname} ...')
    powers.append(np.array(get_td_power(fname)))
    time.append(get_td_time(fname))
    _, t_loc, lp = parse_time_file(fname + '.time')
    local_times.append(np.asarray(t_loc, dtype=float))
    local_pows.append(lp)

# =============================================================================
# Figure 1 – Total power vs time (FOM + ROMs)
# =============================================================================
colors = ['black', 'tab:blue', 'tab:orange', 'tab:green']
ls     = ['-',     '--',       '-.',          ':']

fig, ax = plt.subplots(figsize=(7, 4))
for i in range(len(out_files)):
    ax.plot(time[i], powers[i],
            color=colors[i], ls=ls[i], label=labels[i])

ax.set_xlabel('t (s)')
ax.set_ylabel('Relative power  P(t)')
# ax.set_title('ACP-100 – Sinusoidal $\\Sigma_a$ perturbation (f = 1 Hz, $\\varepsilon$ = 5%)')
ax.grid(True)
ax.legend()
fig.savefig('ACP100_power.pdf', format='pdf')
print('Saved ACP100_power.pdf')

# =============================================================================
# Compute errors  (index 0 = FOM reference)
# =============================================================================
t_fom  = time[0]
P_fom  = powers[0]

results = {}
for i in range(1, len(out_files)):
    t, e, e_mean, e_inf, e_rms = power_errors(t_fom, P_fom, time[i], powers[i])

    t_sp, eL2_t, e_bar_L2, e_max_L2 = spatial_L2_error_timeavg(
        local_times[0], local_pows[0],
        local_times[i], local_pows[i],
        V=None)

    results[labels[i]] = {
        't': t, 'e': e,
        'e_mean': e_mean, 'e_inf': e_inf,
        'e_rms':  e_rms,
        'e_bar_L2': e_bar_L2, 'e_max_L2': e_max_L2,
    }

    print(f"{labels[i]:<24s}  "
          f"e_mean={e_mean:.3e}   e_inf={e_inf:.3e}   "
          f"e_bar_L2={e_bar_L2:.3e}   e_max_L2={e_max_L2:.3e}")

# =============================================================================
# Figure 2 – Relative power error e(t) for the three strategies
# =============================================================================
fig, ax = plt.subplots(figsize=(7, 4))
for (label, r), c, l in zip(results.items(), colors[1:], ls[1:]):
    ax.semilogy(r['t'], r['e'], color=c, ls=l, label=label)

ax.set_xlabel('t (s)')
ax.set_ylabel('Relative power error  $e(t)$')
ax.set_title('ACP-100 – Power error vs time')
ax.grid(True, which='both')
ax.legend()
fig.savefig('ACP100_power_error.pdf', format='pdf')
print('Saved ACP100_power_error.pdf')

# =============================================================================
# Figure 3 – Summary bar chart (e_mean and e_bar_L2)
# =============================================================================
rom_labels = list(results.keys())
e_mean_v   = np.array([r['e_mean']   for r in results.values()])
e_inf_v    = np.array([r['e_inf']    for r in results.values()])
e_bar_L2_v = np.array([r['e_bar_L2'] for r in results.values()])

x   = np.arange(len(rom_labels))
w   = 0.25

fig, ax = plt.subplots(figsize=(7, 4))
ax.bar(x - w,  e_mean_v,   w, label=r'$\bar{e}$ (power)',   color='tab:blue',  alpha=0.85)
ax.bar(x,      e_inf_v,    w, label=r'$e_\infty$ (power)',  color='tab:orange', alpha=0.85)
ax.bar(x + w,  e_bar_L2_v, w, label=r'$\bar{e}_{L^2}$',    color='tab:green',  alpha=0.85)

ax.set_yscale('log')
ax.set_xticks(x)
ax.set_xticklabels(rom_labels)
ax.set_ylabel('Relative error')
ax.set_title('ACP-100 – Strategy comparison')
ax.grid(True, which='both', axis='y')
ax.legend()
fig.savefig('ACP100_strategy_comparison.pdf', format='pdf')
print('Saved ACP100_strategy_comparison.pdf')

# =============================================================================
# Figure 4 (OPTIONAL) – Snapshot convergence for one strategy
# Uncomment and fill in if you run multiple Ns for the same strategy
# =============================================================================

labels_conv = ['Dynamic-Ns=3', 'Dynamic-Ns=5', 'Dynamic-Ns=8', 'Dynamic-Ns=10', 'Dynamic-Ns=15' , 'Dynamic-Ns=20']
files_conv  = [f'ACP100_FE2_td_1cell_rom_dynamic_{n}.out'
                for n in [3, 5, 8, 10, 15, 20]]

Ns_v, e_mean_c, e_inf_c = [], [], []
for fname, lab in zip(files_conv, labels_conv):
    pw = np.array(get_td_power(fname))
    ti = get_td_time(fname)
    _, e, e_m, e_i, _ = power_errors(t_fom, P_fom, ti, pw)
    Ns_v.append(int(lab.split('=')[1]))
    e_mean_c.append(e_m)
    e_inf_c.append(e_i)

fig, ax = plt.subplots(figsize=(6, 4))
ax.semilogy(Ns_v, e_inf_c,  'o-',  label=r'$e_\infty$')
ax.semilogy(Ns_v, e_mean_c, '^-.', label=r'$\bar{e}$')
ax.set_xlabel(r'Number of snapshots $N_s$')
ax.set_ylabel('Relative power error')
ax.set_title('ACP-100 – Dynamic strategy convergence')
ax.set_xticks(Ns_v)
ax.grid(True, which='both')
ax.legend()
fig.savefig('ACP100_snapshot_convergence.pdf', format='pdf')

plt.show()
