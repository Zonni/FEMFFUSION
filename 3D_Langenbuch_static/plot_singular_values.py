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

singular_mono= [
       1.687935e+02, 8.975244e+01, 2.364731e+01, 1.343284e+01,
       5.939941e+00, 4.412835e+00, 3.608049e+00, 2.302933e+00, 2.195260e+00, 1.462544e+00,
       1.401199e+00, 1.232452e+00, 1.076079e+00, 9.260312e-01, 6.604482e-01,
       5.950273e-01, 5.127481e-01, 4.687336e-01, 4.147507e-01,
       3.228161e-01, 2.881719e-01, 2.415692e-01, 2.269684e-01,
       1.868520e-01, 1.822841e-01, 1.510321e-01, 1.464070e-01,
       1.337281e-01, 1.302193e-01, 1.087933e-01, 9.317195e-02,
       8.868275e-02, 8.548673e-02, 7.473541e-02, 5.839699e-02,
       5.551290e-02, 4.753830e-02, 4.423433e-02, 4.189071e-02,
       3.936871e-02, 3.412119e-02, 2.957434e-02, 2.489210e-02,
       2.351208e-02, 2.169397e-02, 2.009456e-02, 1.498584e-02,
       1.263817e-02, 1.107247e-02, 9.563602e-03]


singular_g1 = [
 1.652e+02, 8.746e+01, 2.303e+01, 1.306e+01, 5.476e+00, 3.655e+00, 2.197e+00, 1.419e+00, 1.334e+00,
 1.124e+00, 7.495e-01, 6.946e-01, 4.738e-01, 4.474e-01, 4.046e-01, 2.818e-01, 2.727e-01, 2.037e-01,
 1.698e-01, 1.416e-01, 1.250e-01, 1.120e-01, 1.007e-01, 8.822e-02, 8.153e-02, 7.158e-02, 5.896e-02,
 5.215e-02, 4.340e-02, 3.885e-02, 3.660e-02, 2.813e-02, 2.685e-02, 2.570e-02, 2.362e-02, 1.877e-02,
 1.704e-02, 1.505e-02, 1.423e-02, 1.194e-02, 1.084e-02, 9.782e-03, 9.351e-03, 8.129e-03, 6.656e-03,
 5.549e-03, 4.666e-03, 3.780e-03, 2.834e-03, 2.206e-03]

singular_g2= [
3.499e+01, 2.035e+01, 5.029e+00, 3.085e+00, 1.890e+00, 1.208e+00, 6.095e-01, 5.992e-01, 4.745e-01,
4.597e-01, 3.000e-01, 2.473e-01, 2.215e-01, 1.849e-01, 1.555e-01, 1.414e-01, 1.221e-01, 1.057e-01,
8.936e-02, 8.268e-02, 7.480e-02, 5.832e-02, 5.512e-02, 4.411e-02, 3.770e-02, 3.372e-02, 3.072e-02,
2.960e-02, 2.708e-02, 2.539e-02, 2.153e-02, 1.876e-02, 1.713e-02, 1.600e-02, 1.403e-02, 1.151e-02,
9.947e-03, 9.810e-03, 8.711e-03, 8.292e-03, 6.944e-03, 6.630e-03, 6.238e-03, 5.494e-03, 4.609e-03,
4.148e-03, 3.073e-03, 2.835e-03, 2.227e-03, 1.937e-03]



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
fig.savefig('singular_values.pdf', format='pdf')



