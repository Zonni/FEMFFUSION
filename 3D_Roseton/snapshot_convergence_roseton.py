# -*- coding: utf-8 -*-
"""
Roseton - Convergencia del error con el numero de snapshots.

Genera una figura con el error medio relativo de la potencia total, e_bar,
frente al numero de snapshots N_s, con una curva por estrategia de muestreo.
Sirve a la vez para justificar el N_s adoptado y para comparar las estrategias.

La tolerancia del criterio de energia esta fijada en eps = 1e-5 en todas las
corridas, de modo que el numero de modos lo determina el propio criterio.

Uso:
  1) Ajusta FOM_FILE y el diccionario CASES con los nombres de tus ficheros.
  2) Ejecuta. Imprime una tabla con e_bar, e_inf y e_rms de cada caso y
     guarda la figura roseton_snapshot_convergence.pdf
"""

import sys
sys.path.append('../postprocess')
from utils import get_td_power, get_td_time

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
          'lines.markersize': 5,
          'legend.numpoints': 1,
          'axes.formatter.useoffset': False,
          'figure.autolayout': True,
          }
rcParams.update(params)

# ---------------------------------------------------------------------------
# FICHEROS
#   FOM_FILE : solucion de referencia
#   CASES    : por cada estrategia, un diccionario {N_s: fichero .out}
#              Comenta o borra las entradas que no tengas corridas.
# ---------------------------------------------------------------------------
FOM_FILE = 'roseton_FE2_60p.out'

CASES = {
    'Dynamic': {
        5:  'roseton_FE3_60p_rom_dyn_ns6.out',
        10: 'roseton_FE3_60p_rom_dyn_ns10.out',
        15: 'roseton_FE3_60p_rom_dyn_ns16.out',
        20: 'roseton_FE3_60p_rom_dyn_ns20.out',
        40: 'roseton_FE3_60p_rom_dyn_ns40.out',
    },
    'Static time-sampled': {
        5:  'roseton_FE3_60p_rom_ts_ns6.out',
        10: 'roseton_FE3_60p_rom_ts_ns10.out',
        15: 'roseton_FE3_60p_rom_ts_ns16.out',
        20: 'roseton_FE3_60p_rom_ts_ns20.out',
        40: 'roseton_FE3_60p_rom_ts_ns40.out',
    },
    'Static bank-swept': {
        5:  'roseton_FE3_60p_rom_bs_ns6.out',
        10: 'roseton_FE3_60p_rom_bs_ns10.out',
        15: 'roseton_FE3_60p_rom_bs_ns16.out',
        20: 'roseton_FE3_60p_rom_bs_ns20.out',
        40: 'roseton_FE3_60p_rom_bs_ns40.out',
    },
    'Static IHS-sampled': {
        5:  'roseton_FE3_60p_rom_ihs_ns6.out',
        10: 'roseton_FE3_60p_rom_ihs_ns10.out',
        15: 'roseton_FE3_60p_rom_ihs_ns16.out',
        20: 'roseton_FE3_60p_rom_ihs_ns20.out',
        40: 'roseton_FE3_60p_rom_ihs_ns40.out',
    },
}

# Valor de N_s que se adopta en el resto del estudio (se marca en la figura).
# Pon None si no quieres la linea vertical.
NS_ADOPTED = 20

# ---------------------------------------------------------------------------
# Metricas de error
# ---------------------------------------------------------------------------
def power_errors(t_fom, P_fom, t_rom, P_rom):
    """Error relativo de la potencia total. Devuelve e_bar, e_inf, e_rms."""
    t_fom = np.asarray(t_fom, dtype=float)
    P_fom = np.asarray(P_fom, dtype=float)
    t_rom = np.asarray(t_rom, dtype=float)
    P_rom = np.asarray(P_rom, dtype=float)

    P_rom_i = np.interp(t_fom, t_rom, P_rom)

    t_min = max(t_fom[0], t_rom[0])
    t_max = min(t_fom[-1], t_rom[-1])
    mask = (t_fom >= t_min) & (t_fom <= t_max)
    t = t_fom[mask]
    e = np.abs(P_fom[mask] - P_rom_i[mask]) / np.abs(P_fom[mask])

    T = t[-1] - t[0]
    e_bar = np.trapz(e, t) / T
    e_inf = np.max(e)
    e_rms = np.sqrt(np.trapz(e**2, t) / T)
    return e_bar, e_inf, e_rms


# ---------------------------------------------------------------------------
# Lectura de la referencia
# ---------------------------------------------------------------------------
t_fom = np.asarray(get_td_time(FOM_FILE), dtype=float)
P_fom = np.asarray(get_td_power(FOM_FILE), dtype=float)
print(f"FOM: {FOM_FILE}  ({t_fom.size} instantes)")

# ---------------------------------------------------------------------------
# Bucle sobre estrategias y numero de snapshots
# ---------------------------------------------------------------------------
results = {}

print(f"\n{'Strategy':<22s} {'N_s':>5s} {'e_bar':>11s} {'e_inf':>11s} {'e_rms':>11s}")
for strategy, files in CASES.items():
    ns_list, ebar_list, einf_list, erms_list = [], [], [], []
    for ns in sorted(files):
        fname = files[ns]
        try:
            t_rom = np.asarray(get_td_time(fname), dtype=float)
            P_rom = np.asarray(get_td_power(fname), dtype=float)
        except Exception as err:
            print(f"{strategy:<22s} {ns:>5d}   (no disponible: {err})")
            continue
        if t_rom.size < 2:
            print(f"{strategy:<22s} {ns:>5d}   (fichero vacio o incompleto)")
            continue

        e_bar, e_inf, e_rms = power_errors(t_fom, P_fom, t_rom, P_rom)
        ns_list.append(ns)
        ebar_list.append(e_bar)
        einf_list.append(e_inf)
        erms_list.append(e_rms)
        print(f"{strategy:<22s} {ns:>5d} {e_bar:11.3e} {e_inf:11.3e} {e_rms:11.3e}")

    results[strategy] = {'ns': np.array(ns_list),
                         'e_bar': np.array(ebar_list),
                         'e_inf': np.array(einf_list),
                         'e_rms': np.array(erms_list)}

# ---------------------------------------------------------------------------
# FIGURA: e_bar frente a N_s, una curva por estrategia
# ---------------------------------------------------------------------------
markers = ['o', 's', '^', 'd']

fig, ax = plt.subplots()
for (strategy, r), mk in zip(results.items(), markers):
    if r['ns'].size == 0:
        continue
    ax.semilogy(r['ns'], r['e_bar'], marker=mk, label=strategy)

if NS_ADOPTED is not None:
    ax.axvline(NS_ADOPTED, color='gray', ls=':', lw=1.0)

ax.grid(True, which='both', alpha=0.3)
ax.legend()
ax.set_xlabel('Number of snapshots $N_s$')
ax.set_ylabel(r'Mean relative power error  $\bar{e}$')
fig.savefig('roseton_FE3_60p_snapshot_convergence.pdf', format='pdf')

print("\nFigura guardada en roseton_snapshot_convergence.pdf")
