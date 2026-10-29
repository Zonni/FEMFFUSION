# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
plt.close('all')

params = {'backend': 'pdf', 'font.size': 14, 'axes.labelsize': 12,
          'legend.fontsize': 9, 'xtick.labelsize': 13, 'ytick.labelsize': 13,
          'text.usetex': False, 'lines.linewidth': 1.2, 'lines.markersize': 4,
          'legend.numpoints': 1, 'axes.formatter.useoffset': False,
          'figure.autolayout': True}
rcParams.update(params)

EPS = 1e-5

# Rango de ejes de la figura
M_MAX   = 5       # cortar el eje x en este nº de modos (None = todos)
Y_FLOOR = 1e-12     # suelo del eje y (un par de ordenes por debajo de EPS)


SINGULAR_VALUES_RAW = {
    'Dynamic': {
        'g1': "3.208e+04 1.467e+02 2.828e-01 1.335e-01 5.044e-03 ",
        'g2': "8.823e+03 4.342e+01 1.488e-01 1.110e-01 8.068e-03 ",
    },
    'Static time-sampled': {
        'g1': "2.156e+00 1.089e-02 1.898e-05 1.151e-16 6.314e-17",
        'g2': "5.929e-01 3.202e-03 1.147e-05 3.154e-17 1.998e-17",
    },
    'Static XS-sampled': {
        'g1': "2.156e+00 1.216e-02 1.629e-05 1.287e-07 2.840e-09",
        'g2': "5.929e-01 3.575e-03 9.836e-06 9.316e-08 1.029e-09",
    },

}

# SINGULAR_VALUES_RAW = {
#     'Dynamic': {
#         'g1': "4.531e+04 2.113e+02 2.323e-01 4.520e-02 6.067e-03 4.126e-05 1.890e-05 7.274e-06 2.715e-08 2.249e-08 ",
#         'g2': "1.246e+04 5.417e+01 9.340e-02 1.289e-02 6.017e-03 2.969e-05 1.201e-05 1.910e-06 1.233e-08 5.993e-09 ",
#     },
#     'Static time-sampled': {
#         'g1': "3.049e+00 1.577e-02 2.142e-05 1.352e-16 1.344e-16 7.976e-17 7.271e-17 6.603e-17 6.199e-17 6.077e-17",
#         'g2': "8.386e-01 4.111e-03 6.036e-06 3.660e-17 2.773e-17 2.575e-17 2.288e-17 2.274e-17 2.004e-17 1.677e-17",
#     },

# }

SINGULAR_VALUES = {}
for name, groups in SINGULAR_VALUES_RAW.items():
    SINGULAR_VALUES[name] = {gk: np.array([float(x) for x in txt.split()])
                             for gk, txt in groups.items()}

def captured_energy(sigmas):
    s = np.sort(np.asarray(sigmas, dtype=float))[::-1]
    I = np.cumsum(s**2) / np.sum(s**2)
    return np.arange(1, s.size + 1), I

def modes_for_tolerance(sigmas, eps):
    M, I = captured_energy(sigmas)
    mask = I > (1.0 - eps)
    if not np.any(mask):
        return None
    return int(M[np.argmax(mask)])

print(f"=== Dimension de la base a eps = {EPS:.0e} ===")
print(f"{'Strategy':<22s} {'n_s,g1':>7s} {'n_s,g2':>7s} {'M_1':>5s} {'M_2':>5s}")
for name, groups in SINGULAR_VALUES.items():
    s1, s2 = groups['g1'], groups['g2']
    M1 = modes_for_tolerance(s1, EPS); M2 = modes_for_tolerance(s2, EPS)
    M1s = 'NO' if M1 is None else str(M1); M2s = 'NO' if M2 is None else str(M2)
    print(f"{name:<22s} {s1.size:>7d} {s2.size:>7d} {M1s:>5s} {M2s:>5s}")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for ax, gk, gtitle in zip(axes, ['g1', 'g2'], ['Group 1 (fast)', 'Group 2 (thermal)']):
    for name, groups in SINGULAR_VALUES.items():
        M, I = captured_energy(groups[gk])
        neg = 1.0 - I
        # recortar al rango visible para no arrastrar valores ~1e-200
        if M_MAX is not None:
            sel = M <= M_MAX
            M, neg = M[sel], neg[sel]
        # clip por debajo del suelo para que no se salga del eje
        neg = np.clip(neg, Y_FLOOR, None)
        ax.semilogy(M, neg, marker='o', ms=3, label=name)
    ax.axhline(EPS, color='k', ls='--', lw=1.0, label=fr'$\varepsilon = {EPS:.0e}$')
    ax.set_xlabel('Number of modes $M$')
    ax.set_ylabel(r'Neglected energy $1 - I_g(M)$')
    ax.set_title(gtitle)
    ax.set_ylim(Y_FLOOR, 1.5)
    if M_MAX is not None:
        ax.set_xlim(0.5, M_MAX + 0.5)
    ax.grid(True, which='both', alpha=0.3)
axes[0].legend()
fig.savefig('captured_energy_5.pdf', format='pdf')
print("\nFigura guardada.")