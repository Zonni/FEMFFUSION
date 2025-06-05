# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
#from scipy.interpolate import CubicSpline
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
#%% ===========================================================================

singular_mono= [1.775436e+02,
6.891334e+01,
1.890058e+01,
8.952434e+00,
5.009933e+00,
3.981258e+00,
3.404655e+00,
2.180436e+00,
1.742657e+00,
1.138686e+00,
1.016963e+00,
8.139106e-01,
6.786340e-01,
5.702274e-01,
5.217114e-01,
3.997852e-01,
3.780542e-01,
2.929930e-01,
2.774712e-01,
2.558286e-01,
1.869820e-01,
1.769939e-01,
1.666856e-01,
1.346003e-01,
1.328715e-01,
1.233331e-01,
1.038643e-01,
8.285659e-02,
7.661183e-02,
6.829239e-02,
5.237106e-02,
4.728663e-02,
3.830436e-02,
3.568602e-02,
3.208144e-02,
2.894635e-02,
2.104945e-02,
1.922755e-02,
1.596846e-02,
9.892759e-03]


singular_g1= [1.739250e+02,
6.736036e+01,
1.849202e+01,
8.745760e+00,
4.310716e+00,
3.326564e+00,
1.677485e+00,
1.047450e+00,
9.824588e-01,
7.006491e-01,
5.207938e-01,
4.621116e-01,
3.356224e-01,
2.688243e-01,
2.084048e-01,
1.730168e-01,
1.541189e-01,
1.299479e-01,
1.233722e-01,
1.138998e-01,
9.886145e-02,
7.826918e-02,
6.123798e-02,
5.634990e-02,
5.187163e-02,
4.685903e-02,
4.038163e-02,
3.866690e-02,
2.935237e-02,
2.350924e-02,
1.808187e-02,
1.618317e-02,
1.350525e-02,
1.162401e-02,
9.185913e-03,
7.859572e-03,
6.448614e-03,
4.914308e-03,
3.926436e-03,
2.968753e-03]

singular_g2=[3.593438e+01,1.471815e+01,3.956834e+00,1.963859e+00,
1.483564e+00,9.730401e-01,
6.922082e-01,
4.499916e-01,
3.482754e-01,
3.139915e-01,
2.055090e-01,
1.721768e-01,
1.627093e-01,
1.339779e-01,
1.167284e-01,
9.902288e-02,
8.798713e-02,
7.264899e-02,
6.676163e-02,
6.394949e-02,
4.757877e-02,
4.048654e-02,
3.456507e-02,
2.923634e-02,
2.623970e-02,
2.126641e-02,
2.046316e-02,
1.589332e-02,
1.459203e-02,
1.248557e-02,
1.205236e-02,
9.202394e-03,
8.530666e-03,
8.047719e-03,
7.618637e-03,
6.595239e-03,
4.957903e-03,
3.598188e-03,
2.835058e-03,
2.708163e-03]

indices = list(range(1,len(singular_g2)+1))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.semilogy(indices, singular_mono, 'o-', label='Monolithic')
ax.semilogy(indices, singular_g1, 'o-', label='Group 1')
ax.semilogy(indices, singular_g2, 'o-', label='Group 2')
ax.grid(True)
ax.legend()
ax.set_xlabel('Index')
ax.set_ylabel('Singular Value')
fig.savefig('snapshots.pdf', format='pdf')

