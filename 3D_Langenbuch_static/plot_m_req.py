# -*- coding: utf-8 -*-
import sys
sys.path.append('../postprocess')
from utils import parse_file_same_line
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
# LUPOD POINTS 


m_req = np.array(list(range(2, 101, 2)))

max_keff = []
mean_keff = []
mean_phi = []
max_phi = []

for m in m_req:
    out_file = 'mreq/3D_Langenbuch_mreq' + str(m) + '.out'

    
    mean_keff.append(parse_file_same_line(out_file, begin='Mean Delta Keff (pcm):')[0])
    max_keff.append(parse_file_same_line(out_file, begin='Max  Delta Keff (pcm):')[0])
    mean_phi.append(parse_file_same_line(out_file, begin='Mean RMS Phi (%):')[0])
    max_phi.append(parse_file_same_line(out_file, begin='Max  RMS Phi (%):')[0])

## Plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.semilogy(m_req, mean_keff, 'o-', label='Mean Error')
ax.semilogy(m_req, max_keff, 'o-', label='Max Error')
ax.grid(True)
ax.legend()
ax.set_xlabel('Number of modes retained')
ax.set_ylabel('$\Delta K$eff Error (pcm)')
fig.savefig('keff_error_mreq.pdf', format='pdf')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.semilogy(m_req, mean_phi, 'o-', label='Mean RMS Error')
ax.semilogy(m_req, max_phi, 'o-', label='Max RMS Error')
ax.grid(True)
ax.legend()
ax.set_xlabel('Number of modes retained')
ax.set_ylabel('$\phi$ Error (%)')
fig.savefig('phi_error_mreq.pdf', format='pdf')


singular_g1 = [
 1.652e+02, 8.746e+01, 2.303e+01, 1.306e+01, 5.476e+00, 3.655e+00, 2.197e+00, 1.419e+00, 1.334e+00,
 1.124e+00, 7.495e-01, 6.946e-01, 4.738e-01, 4.474e-01, 4.046e-01, 2.818e-01, 2.727e-01, 2.037e-01,
 1.698e-01, 1.416e-01, 1.250e-01, 1.120e-01, 1.007e-01, 8.822e-02, 8.153e-02, 7.158e-02, 5.896e-02,
 5.215e-02, 4.340e-02, 3.885e-02, 3.660e-02, 2.813e-02, 2.685e-02, 2.570e-02, 2.362e-02, 1.877e-02,
 1.704e-02, 1.505e-02, 1.423e-02, 1.194e-02, 1.084e-02, 9.782e-03, 9.351e-03, 8.129e-03, 6.656e-03,
 5.549e-03, 4.666e-03, 3.780e-03, 2.834e-03, 2.206e-03]
singular_g1 = np.array(singular_g1)

singular_g2= [
3.499e+01, 2.035e+01, 5.029e+00, 3.085e+00, 1.890e+00, 1.208e+00, 6.095e-01, 5.992e-01, 4.745e-01,
4.597e-01, 3.000e-01, 2.473e-01, 2.215e-01, 1.849e-01, 1.555e-01, 1.414e-01, 1.221e-01, 1.057e-01,
8.936e-02, 8.268e-02, 7.480e-02, 5.832e-02, 5.512e-02, 4.411e-02, 3.770e-02, 3.372e-02, 3.072e-02,
2.960e-02, 2.708e-02, 2.539e-02, 2.153e-02, 1.876e-02, 1.713e-02, 1.600e-02, 1.403e-02, 1.151e-02,
9.947e-03, 9.810e-03, 8.711e-03, 8.292e-03, 6.944e-03, 6.630e-03, 6.238e-03, 5.494e-03, 4.609e-03,
4.148e-03, 3.073e-03, 2.835e-03, 2.227e-03, 1.937e-03]
singular_g2 = np.array(singular_g2)


tot_g1 = np.linalg.norm(singular_g1)
tot_g2 = np.linalg.norm(singular_g2)

leftout_energy_g1 = []
leftout_energy_g2 = []
for m in m_req:
    leftout_energy_g1.append(1.0 - (np.linalg.norm(singular_g1[0:m])/tot_g1))
    leftout_energy_g2.append(1.0 - (np.linalg.norm(singular_g2[0:m])/tot_g2))
    
    
## Plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.loglog(leftout_energy_g1, mean_keff, 'o-',  label='Group 1 - Mean Error')
ax.loglog(leftout_energy_g1, max_keff, 'o-',   label='Group 1 - Max Error')
ax.loglog(leftout_energy_g2, mean_keff, 'x--', label='Group 2 - Mean Error')
ax.loglog(leftout_energy_g2, max_keff, 'x--',  label='Group 2 - Max Error')
ax.grid(True)
ax.legend()
ax.set_xlabel('Energy Fraction Left Out per Group')
ax.set_ylabel('$\Delta K$eff Error (pcm)')
fig.savefig('keff_error_energy.pdf', format='pdf')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.loglog(leftout_energy_g1, mean_phi, 'o-',  label='Group 1 - Mean RMS Error')
ax.loglog(leftout_energy_g1, max_phi,  'o-',  label='Group 1 - Max RMS Error')
ax.loglog(leftout_energy_g2, mean_phi, 'x--', label='Group 2 - Mean RMS Error')
ax.loglog(leftout_energy_g2, max_phi,  'x--', label='Group 2 - Max RMS Error')
ax.grid(True)
ax.legend()
ax.set_xlabel('Energy Fraction Left Out per Group')
ax.set_ylabel('$\phi$ Error ')
fig.savefig('phi_error_energy.pdf', format='pdf')

