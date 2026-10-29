# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal.
"""
import sys
sys.path.append('../postprocess')
from utils import get_td_power, get_td_time, parse_time_file, power_at_time, spatial_L2_error_timeavg, spatial_L2_error, power_errors
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.interpolate import CubicSpline
#import scipy.io as sio
plt.close('all')
#%% ===========================================================================

params = {'backend': 'pdf',
#          'font.family': 'serif',
          'font.size': 14,
          'axes.labelsize': 12,
          'legend.fontsize': 10,
          'xtick.labelsize': 13,
          'ytick.labelsize': 13,
          'text.usetex': False,
          'lines.linewidth': 1.2,
          'lines.markersize': 3,
          'lines.markeredgewidth': 1,
          'legend.numpoints': 1, 
          'axes.formatter.useoffset': False,
          'figure.autolayout': True,
          }

rcParams.update(params)


# labels = ['FOM','Dynamic-Ns=5','Dynamic-Ns=8','Dynamic-Ns=10','Dynamic-Ns=15','Dynamic-Ns=20','Dynamic-Ns=40']
# out_files =['3D_Langenbuch_ds_fe3_01.out','3D_Langenbuch_fe3_rom_dynamic_5.out','3D_Langenbuch_fe3_rom_dynamic_8.out',
#             '3D_Langenbuch_fe3_rom_dynamic_10.out','3D_Langenbuch_fe3_rom_dynamic_15.out','3D_Langenbuch_fe3_rom_dynamic_20.out',
#             '3D_Langenbuch_fe3_rom_dynamic_40.out'
#     ]

# POD
# labels = ['FOM','Dynamic-Ns=10','Static Time-Sampled-Ns=10','Static Bank-Swept-Ns=10','Static IHS-Sampled-Ns=10']
# out_files =['3D_Langenbuch_ds_fe3_01.out','3D_Langenbuch_fe3_rom_dynamic.out','3D_Langenbuch_fe3_rom_static_time_sampled.out',
#             '3D_Langenbuch_fe3_rom_static_bank_swept.out','3D_Langenbuch_fe3_rom_static_bank_ihs.out']


# RPOD
labels = ['FOM', 'POD fixed (Static-Bank-Swept)', 'RPOD (Static-Bank-Swept)']
out_files = [
    '3D_Langenbuch_ds_fe3_01.out',
    '3D_Langenbuch_fe3_rom_static_bank_swept.out',
    '3D_Langenbuch_fe3_rrom_static_bank.out'    
]

# labels = ['FOM','POD-3','POD-5','POD-10']
# out_files =['3D_Langenbuch_ds_bars_fe3.out',
#     '3D_Langenbuch_rom3_bars_fe3.out',
#         '3D_Langenbuch_rom5_bars_fe3.out',
#         '3D_Langenbuch_rom10_bars_fe3.out'
#     ]

# labels = ['FOM','ROM-10-Static']
# out_files =['3D_Langenbuch_ds_fe3_01.out','3D_Langenbuch_rom_gg.out'
#     ]




# labels = ['FOM','ROM-10-Static','ROM-10-Dyn']
# out_files =['3D_Langenbuch_ds_bars_fe1.out','3D_Langenbuch_rom_bars.out','3D_Langenbuch_rom10_dyn_fe1.out'
#     ]

# ---------------------------------------------------------------------------
# Read files .out and .out.time
# ---------------------------------------------------------------------------

    
time = []
powers = []
local_times = []
local_pows = []

for i in range(len(out_files)):
    print(i)
    power = get_td_power(out_files[i])
    powers.append(np.array(power))
    time.append(get_td_time(out_files[i]))

    _, t_loc, local_pow = parse_time_file(out_files[i]+ '.time')
    local_times.append(np.asarray(t_loc, dtype=float)) 
    local_pows.append(local_pow)


# ---------------------------------------------------------------------------
# Check sizes
# ---------------------------------------------------------------------------
# print("=== Comprobación de dimensiones ===")
# for i in range(len(out_files)):
#     lp = np.asarray(local_pows[i])
#     lt = np.asarray(local_times[i])
#     print(f"\n{out_files[i]}")
#     print(f"  local_pow shape : {lp.shape}")
#     print(f"  local_time len  : {lt.shape}")
#     print(f"  time(.out) len  : {len(time[i])}")
#     # ¿el nº de filas de local_pow coincide con su eje temporal?
#     print(f"  filas==tiempo?  : {lp.shape[0] == lt.shape[0]}")

# # ¿FOM y ROM tienen el mismo nº de celdas?
# print("\n=== Celdas ===")
# n_cells_fom = np.asarray(local_pows[0]).shape[1]
# for i in range(len(out_files)):
#     nc = np.asarray(local_pows[i]).shape[1]
#     print(f"{out_files[i]}: n_cells = {nc}  ", 
#           "OK" if nc == n_cells_fom else "DISTINTO al FOM")
    


#%% ===========================================================================



# ---------------------------------------------------------------------------
# Plot power
# ---------------------------------------------------------------------------

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
error_power=[];
# ax.plot(time, power, label='Δt = 0.1 s')
for i in range(len(out_files)):
    ax.plot(time[i], powers[i],  label=labels[i])

    
# ax.axvspan(0.0, 26.7, alpha=0.12, color='tab:blue',  label='C1 withdrawal')
# ax.axvspan(7.5, 47.5, alpha=0.12, color='tab:red',   label='C2 insertion')
ax.grid(True)
ax.legend()
ax.set_xlabel('t (s)')
ax.set_ylabel('Relative Power P(t)')



fig.savefig('langenbuch_fom_rom.pdf', format='pdf')


# ---------------------------------------------------------------------------
# Bucle único de errores (índice 0 = FOM, referencia)
# ---------------------------------------------------------------------------
t_fom  = time[0]
P_fom  = powers[0]
lp_fom = np.asarray(local_pows[0])          # (n_time_fom, n_cells)

   
results = {}
for i in range(1, len(out_files)):
    # --- métricas temporales sobre la potencia total ---
    t, e, e_mean, e_inf, e_rms = power_errors(t_fom, P_fom, time[i], powers[i])

    # --- instante de máximo error y métrica espacial en ese instante ---
    t_star = t[np.argmax(e)]
    p_fom_star = power_at_time(local_times[0], np.asarray(local_pows[0]), t_star, method='linear')
    p_rom_star = power_at_time(local_times[i], np.asarray(local_pows[i]), t_star, method='linear')
    # eL2 = spatial_L2_error(p_rom_star, p_fom_star, V=None)
    t_sp, eL2_t, e_bar_L2, e_max_L2 = spatial_L2_error_timeavg(
    local_times[0], local_pows[0],
    local_times[i], local_pows[i], V=None)

    results[labels[i]] = {'t': t, 'e': e,
                          'e_mean': e_mean, 'e_inf': e_inf,
                          'e_rms': e_rms, 'e_bar_L2': e_bar_L2, 'e_max_L2': e_max_L2 }

    print(f"{labels[i]:<16s}  "
          f"e_mean = {e_mean:.3e}   e_inf = {e_inf:.3e}   e_rms = {e_rms:.3e}   "
             f"e_bar_L2={e_bar_L2:.3e}    e_max_L2={e_max_L2:.3e}" )



    
# ---------------------------------------------------------------------------
#  PLOT e(t)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots()
for label, r in results.items():
    ax.semilogy(r['t'], r['e'], label=label)
ax.grid(True, which='both')
ax.legend()
ax.set_xlabel('t (s)')
ax.set_ylabel('Relative power error  e(t)')
fig.savefig('langenbuch_power_error.pdf', format='pdf')


# ---------------------------------------------------------------------------
#  PLOT convergencia de snapshots (Static IHS-Sampled, epsilon = 1e-5)
# ---------------------------------------------------------------------------
# Numero de snapshots extraido de cada label 'Static IHS-Sampled-Ns=XX'
N_s     = np.array([int(lab.split('=')[1]) for lab in results.keys()])
e_inf_v  = np.array([r['e_inf']  for r in results.values()])
e_rms_v  = np.array([r['e_rms']  for r in results.values()])
e_mean_v = np.array([r['e_mean'] for r in results.values()])

# ordenar por N_s por si el dict no viene en orden
order = np.argsort(N_s)
N_s, e_inf_v, e_rms_v, e_mean_v = N_s[order], e_inf_v[order], e_rms_v[order], e_mean_v[order]

fig, ax = plt.subplots()
ax.semilogy(N_s, e_inf_v,  'o-',  label=r'$e_{\infty}$')
ax.semilogy(N_s, e_rms_v,  's--', label=r'$e_{\mathrm{rms}}$')
ax.semilogy(N_s, e_mean_v, '^-.', label=r'$\bar{e}$')
ax.set_ylim(1e-3, 1e-2) 
ax.set_xticks(N_s)
ax.set_xticklabels([str(n) for n in N_s])
ax.grid(True, which='both')
ax.legend()
ax.set_xlabel(r'Number of snapshots $N_s$')
ax.set_ylabel('Relative power error')
fig.savefig('langenbuch_snapshot_convergence.pdf', format='pdf')


