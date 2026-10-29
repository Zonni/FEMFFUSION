# -*- coding: utf-8 -*-
"""
Comparacion FOM / POD base fija / RPOD para el reactor Roseton (hexagonal 3D).

Transitorio:
  Banco 1: 125 -> 275 cm en 0.15 s   (extraccion rapida, tipo eyeccion)
  Banco 2: 225 -> ...  desde 1.0 s   (insercion lenta)
  Banco 3: fijo en 100 cm
  t_end = 7.0 s

Reinicios naturales del RPOD: t = 0.15 s y t = 1.0 s

Metricas:
  e(t)      error relativo instantaneo de la potencia total
  e_mean    (1/T) int e(t) dt
  e_inf     max_t e(t)
  e_rms     sqrt((1/T) int e(t)^2 dt)
  e_bar_L2  (1/T) int e_L2(t) dt   (error espacial L2 por celdas, normalizado)
  e_max_L2  max_t e_L2(t)
"""

import sys
sys.path.append('../postprocess')
from utils import get_td_power, get_td_time, parse_time_file

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

plt.close('all')

# ---------------------------------------------------------------------------
# Parametros de figura
# ---------------------------------------------------------------------------
params = {'backend': 'pdf',
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

# ---------------------------------------------------------------------------
# CASOS  (indice 0 = FOM, referencia)
# ---------------------------------------------------------------------------
# labels = ['FOM',
#           'POD fixed (Static-Bank-Swept)',
#           'RPOD (Static-Bank-Swept)']

# out_files = ['roseton_FE2_60p.out',
#               'roseton_FE2_60p_rom.out',
#               'roseton_FE2_60p_rrom.out']

labels = ['FOM',
          'POD fixed (Static-Bank-Swept)',
          'RPOD (Static-Bank-Swept)']

out_files = ['roseton_FE3_60p.out',
              'roseton_FE3_60p_rom_bs_ns20.out',
              'roseton_FE3_60p_rrom.out']


# labels = ['FOM',
#           'POD fixed (Static-Bank-Swept)',
#           'RPOD (Static-Bank-Swept)']

# out_files = ['roseton_FE2_60p_c2.out',
#              'roseton_FE2_60p_rom_c2.out',
#              'roseton_FE2_60p_rrom_c2.out']

# Instantes de reinicio del RPOD (para anotarlos en las figuras)
RESTART_TIMES = [1.0]

# Escala logaritmica en el eje de tiempos (util aqui: el pico esta en 0.15 s
# y el transitorio dura 7 s). Pon False para escala lineal.
LOG_TIME = False

# ---------------------------------------------------------------------------
# Lectura de ficheros
# ---------------------------------------------------------------------------
powers = []
time = []
local_times = []
local_pows = []

for f in out_files:
    powers.append(np.asarray(get_td_power(f), dtype=float))
    time.append(np.asarray(get_td_time(f), dtype=float))
    _, t_loc, lp = parse_time_file(f + '.time')
    local_times.append(np.asarray(t_loc, dtype=float))
    local_pows.append(np.asarray(lp, dtype=float))

# ---------------------------------------------------------------------------
# Comprobacion de dimensiones
# ---------------------------------------------------------------------------
print("=== Comprobacion de dimensiones ===")
for i, f in enumerate(out_files):
    lp = local_pows[i]
    lt = local_times[i]
    print(f"\n{f}")
    print(f"  local_pow shape : {lp.shape}")
    print(f"  local_time len  : {lt.shape}")
    print(f"  time(.out) len  : {len(time[i])}")
    print(f"  filas==tiempo?  : {lp.shape[0] == lt.shape[0]}")
    print(f"  n_cells         : {lp.shape[1] if lp.ndim == 2 else 'N/A'}")

n_cells_fom = local_pows[0].shape[1]
print("\n=== Celdas ===")
for i, f in enumerate(out_files):
    nc = local_pows[i].shape[1]
    print(f"{f}: n_cells = {nc}  ", "OK" if nc == n_cells_fom else "DISTINTO al FOM")

# ---------------------------------------------------------------------------
# Funciones de error
# ---------------------------------------------------------------------------
def power_errors(t_fom, P_fom, t_rom, P_rom):
    """Metricas temporales del error relativo de la potencia total."""
    t_fom = np.asarray(t_fom, dtype=float)
    P_fom = np.asarray(P_fom, dtype=float)
    P_rom_i = np.interp(t_fom, np.asarray(t_rom, dtype=float),
                        np.asarray(P_rom, dtype=float))

    t_min = max(t_fom[0], t_rom[0])
    t_max = min(t_fom[-1], t_rom[-1])
    mask = (t_fom >= t_min) & (t_fom <= t_max)
    t = t_fom[mask]
    e = np.abs(P_fom[mask] - P_rom_i[mask]) / np.abs(P_fom[mask])

    T = t[-1] - t[0]
    e_mean = np.trapz(e, t) / T
    e_inf = np.max(e)
    e_rms = np.sqrt(np.trapz(e**2, t) / T)
    return t, e, e_mean, e_inf, e_rms


def spatial_L2_error_timeavg(t_fom, lp_fom, t_rom, lp_rom, V=None):
    """
    Error L2 espacial a lo largo del transitorio, sobre distribuciones
    normalizadas por su potencia total (mide la FORMA espacial).
    """
    t_fom = np.asarray(t_fom, dtype=float)
    lp_fom = np.asarray(lp_fom, dtype=float)
    t_rom = np.asarray(t_rom, dtype=float)
    lp_rom = np.asarray(lp_rom, dtype=float)

    if lp_fom.ndim != 2 or lp_rom.ndim != 2:
        raise ValueError(f"lp_fom {lp_fom.shape} y lp_rom {lp_rom.shape} deben ser 2D")
    if lp_fom.shape[1] != lp_rom.shape[1]:
        raise ValueError(f"n_cells distinto: FOM={lp_fom.shape[1]}, ROM={lp_rom.shape[1]}")

    n_cells = lp_fom.shape[1]
    V = np.ones(n_cells) if V is None else np.asarray(V, dtype=float)

    t_min = max(t_fom[0], t_rom[0])
    t_max = min(t_fom[-1], t_rom[-1])
    mask = (t_fom >= t_min) & (t_fom <= t_max)
    t_common = t_fom[mask]

    lp_rom_i = np.zeros((t_common.size, n_cells))
    for c in range(n_cells):
        lp_rom_i[:, c] = np.interp(t_common, t_rom, lp_rom[:, c])
    lp_fom_c = lp_fom[mask]

    # normalizar por la potencia total de cada instante (compara solo la forma)
    fom_sum = lp_fom_c.sum(axis=1, keepdims=True)
    rom_sum = lp_rom_i.sum(axis=1, keepdims=True)
    fom_sum = np.where(fom_sum > 1e-30, fom_sum, 1.0)
    rom_sum = np.where(rom_sum > 1e-30, rom_sum, 1.0)
    lp_fom_n = lp_fom_c / fom_sum
    lp_rom_n = lp_rom_i / rom_sum

    num = np.sqrt(np.sum((lp_rom_n - lp_fom_n)**2 * V, axis=1))
    den = np.sqrt(np.sum(lp_fom_n**2 * V, axis=1))

    good = den > 1e-30
    if not np.any(good):
        raise ValueError("Denominador nulo en todos los instantes.")
    eL2_t = np.full(t_common.shape, np.nan)
    eL2_t[good] = num[good] / den[good]

    t_g = t_common[good]
    e_g = eL2_t[good]
    T = t_g[-1] - t_g[0]
    e_bar_L2 = np.trapz(e_g, t_g) / T
    e_max_L2 = np.max(e_g)
    return t_common, eL2_t, e_bar_L2, e_max_L2


def mark_restarts(ax):
    """Lineas verticales en los instantes de reinicio del RPOD."""
    for k, tr in enumerate(RESTART_TIMES):
        ax.axvline(tr, color='gray', ls=':', lw=1.0,
                   label='RPOD restarts' if k == 0 else None)


# ---------------------------------------------------------------------------
# Bucle de errores
# ---------------------------------------------------------------------------
t_fom = time[0]
P_fom = powers[0]

results = {}
print("\n=== Errores ===")
for i in range(1, len(out_files)):
    t, e, e_mean, e_inf, e_rms = power_errors(t_fom, P_fom, time[i], powers[i])
    t_sp, eL2_t, e_bar_L2, e_max_L2 = spatial_L2_error_timeavg(
        local_times[0], local_pows[0],
        local_times[i], local_pows[i], V=None)

    results[labels[i]] = {'t': t, 'e': e,
                          'e_mean': e_mean, 'e_inf': e_inf, 'e_rms': e_rms,
                          't_sp': t_sp, 'eL2_t': eL2_t,
                          'e_bar_L2': e_bar_L2, 'e_max_L2': e_max_L2}

    print(f"{labels[i]:<32s} e_mean={e_mean:.3e}  e_inf={e_inf:.3e}  "
          f"e_rms={e_rms:.3e}  e_bar_L2={e_bar_L2:.3e}  e_max_L2={e_max_L2:.3e}")

# ---------------------------------------------------------------------------
# FIGURA 1: potencia total
# ---------------------------------------------------------------------------
fig, ax = plt.subplots()
for i in range(len(out_files)):
    ax.plot(time[i], powers[i], label=labels[i])
mark_restarts(ax)
if LOG_TIME:
    ax.set_xscale('log')
ax.grid(True, which='both', alpha=0.3)
ax.legend()
ax.set_xlabel('t (s)')
ax.set_ylabel('Relative Power P(t)')
fig.savefig('roseton_fom_rom_rrom.pdf', format='pdf')

# ---------------------------------------------------------------------------
# FIGURA 2: error relativo de potencia
# ---------------------------------------------------------------------------
fig, ax = plt.subplots()
for label, r in results.items():
    ax.semilogy(r['t'], r['e'], label=label)
mark_restarts(ax)
if LOG_TIME:
    ax.set_xscale('log')
ax.grid(True, which='both', alpha=0.3)
ax.legend()
ax.set_xlabel('t (s)')
ax.set_ylabel('Relative power error  e(t)')
fig.savefig('roseton_power_error.pdf', format='pdf')

# ---------------------------------------------------------------------------
# FIGURA 3: error espacial L2
# ---------------------------------------------------------------------------
fig, ax = plt.subplots()
for label, r in results.items():
    ax.semilogy(r['t_sp'], r['eL2_t'], label=label)
mark_restarts(ax)
if LOG_TIME:
    ax.set_xscale('log')
ax.grid(True, which='both', alpha=0.3)
ax.legend()
ax.set_xlabel('t (s)')
ax.set_ylabel(r'Spatial $L^2$ power error  $e_{L^2}(t)$')
fig.savefig('roseton_spatial_error.pdf', format='pdf')

print("\nFiguras guardadas: roseton_fom_rom_rrom.pdf, roseton_power_error.pdf, "
      "roseton_spatial_error.pdf")
