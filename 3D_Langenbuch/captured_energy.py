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
M_MAX   = 10       # cortar el eje x en este nº de modos (None = todos)
Y_FLOOR = 1e-12     # suelo del eje y (un par de ordenes por debajo de EPS)

# SINGULAR_VALUES_RAW = {
#     'Dynamic': {
#         'g1': "1.768e+05 1.905e+03 9.080e+02 7.813e+02 1.998e+02 1.516e+02 7.150e+01 5.255e+01 1.198e+01 6.098e+00 5.778e+00 1.332e+00 3.099e-01 1.852e-01 1.541e-01 9.446e-02 8.186e-02 7.246e-02 6.193e-02 2.285e-02 1.747e-02 7.287e-03 5.356e-03 4.292e-03 2.198e-03 1.325e-03 9.180e-04 2.595e-04 2.314e-04 1.561e-04 1.133e-04 7.281e-05 2.429e-05 1.189e-05 8.828e-06 5.099e-06 3.725e-07 3.582e-08 8.780e-09 7.655e-10",
#         'g2': "3.516e+04 4.392e+02 2.065e+02 1.924e+02 5.862e+01 4.404e+01 2.464e+01 2.093e+01 3.548e+00 2.221e+00 1.963e+00 3.969e-01 1.157e-01 9.752e-02 7.911e-02 6.511e-02 5.439e-02 4.394e-02 4.216e-02 1.833e-02 8.870e-03 3.655e-03 2.621e-03 1.628e-03 9.349e-04 5.535e-04 3.769e-04 1.452e-04 1.158e-04 6.347e-05 4.702e-05 3.491e-05 1.496e-05 6.890e-06 6.161e-06 4.229e-06 2.202e-07 1.733e-08 3.122e-09 2.164e-10",
#     },
#     'Static time-sampled': {
#         'g1': "6.202e+00 7.727e-02 5.402e-02 2.577e-02 7.488e-03 5.402e-03 2.789e-03 1.625e-03 3.130e-04 1.840e-04 1.499e-04 1.494e-05 6.201e-06 4.782e-06 3.860e-06 2.733e-06 1.883e-06 1.710e-06 5.129e-07 5.301e-08 2.100e-08 1.607e-08 1.343e-08 1.062e-08 6.545e-09 3.667e-09 3.527e-09 2.344e-09 1.352e-09 6.692e-10 7.634e-11 1.873e-11 2.452e-15 6.146e-16 6.146e-16 6.146e-16 6.146e-16 6.146e-16 6.146e-16 6.146e-16",
#         'g2': "1.233e+00 1.870e-02 1.177e-02 6.216e-03 2.182e-03 1.533e-03 9.593e-04 6.581e-04 8.880e-05 6.912e-05 5.109e-05 5.890e-06 3.524e-06 2.941e-06 2.338e-06 2.120e-06 1.570e-06 1.331e-06 4.262e-07 3.708e-08 1.207e-08 1.038e-08 7.575e-09 3.474e-09 3.131e-09 2.202e-09 1.274e-09 6.110e-10 4.104e-10 1.552e-10 2.488e-11 1.743e-11 4.863e-16 1.220e-16 1.220e-16 1.220e-16 1.220e-16 1.220e-16 1.220e-16 1.220e-16",
#     },
#     'Static bank-swept': {
#         'g1': "6.202e+00 1.001e-01 7.234e-02 2.823e-02 1.942e-02 1.294e-02 8.015e-03 5.929e-03 4.441e-03 3.852e-03 3.184e-03 2.735e-03 1.985e-03 1.486e-03 1.110e-03 8.635e-04 6.607e-04 6.412e-04 4.806e-04 9.646e-06 7.111e-06 5.636e-06 3.213e-06 3.040e-06 2.575e-06 2.325e-06 2.186e-06 2.144e-06 1.885e-06 1.681e-06 1.609e-06 1.548e-06 1.363e-06 1.223e-06 1.105e-06 8.938e-07 8.278e-07 4.277e-08 3.418e-08 2.630e-15 ",
#         'g2': "1.233e+00 2.351e-02 1.541e-02 7.479e-03 4.864e-03 3.344e-03 2.179e-03 1.646e-03 1.203e-03 1.152e-03 8.834e-04 8.193e-04 6.312e-04 5.444e-04 4.141e-04 3.627e-04 2.907e-04 2.697e-04 2.121e-04 6.895e-05 4.979e-05 2.417e-06 2.264e-06 2.190e-06 1.768e-06 1.694e-06 1.683e-06 1.588e-06 1.342e-06 1.296e-06 1.158e-06 1.074e-06 9.594e-07 9.104e-07 8.874e-07 6.882e-07 6.540e-07 3.887e-08 2.775e-08 5.277e-16",
#     },
#     'Static IHS-sampled': {
#         'g1': "6.202e+00 1.082e-01 8.064e-02 3.918e-02 2.326e-02 1.835e-02 1.030e-02 8.584e-03 5.241e-03 4.322e-03 3.853e-03 3.124e-03 2.560e-03 2.216e-03 1.517e-03 1.159e-03 9.229e-04 7.136e-04 6.257e-04 3.383e-04 2.439e-04 6.728e-05 5.011e-05 2.195e-05 1.903e-05 1.203e-05 1.001e-05 8.906e-06 7.434e-06 6.233e-06 5.501e-06 4.673e-06 4.244e-06 3.870e-06 3.435e-06 2.957e-06 2.882e-06 1.943e-06 1.713e-06 1.367e-06",
#         'g2': "1.231e+00 2.265e-02 1.918e-02 1.035e-02 5.850e-03 4.665e-03 2.826e-03 2.300e-03 1.566e-03 1.332e-03 1.050e-03 9.123e-04 8.037e-04 7.407e-04 5.778e-04 4.892e-04 3.878e-04 2.937e-04 2.610e-04 8.681e-05 7.067e-05 4.364e-05 3.548e-05 1.803e-05 1.385e-05 6.529e-06 5.950e-06 4.251e-06 4.087e-06 3.274e-06 3.162e-06 2.979e-06 2.251e-06 2.230e-06 2.033e-06 1.881e-06 1.609e-06 1.354e-06 9.657e-07 9.185e-07",
#     },
# }

SINGULAR_VALUES_RAW = {
    'Dynamic': {
        'g1': "1.008e+05 1.036e+03 5.308e+02 4.235e+02 1.357e+02 1.005e+02 5.419e+01 3.672e+01 7.796e-01 3.144e-02",
        'g2': "2.005e+04 2.372e+02 1.219e+02 1.042e+02 3.913e+01 2.984e+01 2.007e+01 1.351e+01 2.167e-01 8.940e-03",
    },
    'Static time-sampled': {
        'g1': "3.101e+00 3.975e-02 2.744e-02 1.426e-02 4.323e-03 3.371e-03 1.747e-03 1.102e-03 6.722e-06 6.362e-17",
        'g2': "6.163e-01 9.644e-03 5.943e-03 3.445e-03 1.247e-03 9.530e-04 6.070e-04 4.511e-04 2.272e-06 1.248e-17",
    },
    'Static bank-swept': {
        'g1': "3.101e+00 4.961e-02 3.542e-02 1.399e-02 9.088e-03 6.977e-03 3.989e-03 2.956e-03 2.175e-03 7.672e-17 ",
        'g2': "6.163e-01 1.171e-02 7.532e-03 3.802e-03 2.299e-03 1.767e-03 1.102e-03 9.344e-04 6.920e-04 1.570e-17",
    },
    'Static IHS-sampled': {
        'g1': "3.101e+00 5.788e-02 2.756e-02 2.343e-02 1.128e-02 8.535e-03 5.871e-03 3.082e-03 1.613e-03 7.074e-04",
        'g2': "6.154e-01 1.215e-02 6.766e-03 6.146e-03 2.798e-03 2.157e-03 1.584e-03 9.361e-04 5.455e-04 2.720e-04",
    },
}

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
fig.savefig('captured_energy_10.pdf', format='pdf')
print("\nFigura guardada.")