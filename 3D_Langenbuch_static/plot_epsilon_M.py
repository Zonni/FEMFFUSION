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

# Vector de epsilon en notación científica
eps=(
 
)

# Archivo de parámetros
param_file="3D_Langenbuch_static/3D_Langenbuch_POD_group_wise.prm"

# Loop a través del vector de epsilons


eps = np.array([
  3e-1, 1e-1,
  3e-2, 1e-2,
  3e-3, 1e-3,
  3e-4, 1e-4,
  3e-5, 1e-5,
  3e-6, 1e-6,
  3e-7, 1e-7,
  3e-8, 1e-8,
  3e-9, 1e-9,
  3e-10, 1e-10,
  3e-11, 1e-11,
  3e-12, 1e-12,
  3e-13, 1e-13,
  3e-14, 1e-14,
  3e-15, 1e-15,
  3e-16, 1e-16,
])

max_keff = []
mean_keff = []
mean_phi = []
max_phi = []
rom_dim = []
# Iterate over the eps values
for e in eps:
    # Convert float to scientific notation string, then replace '-' with 'm'
    e_label = f"{e:.0e}".replace('-', 'm')
    out_file = f"epsilonM/3D_Langenbuch_eps_{e_label}.out"
    print(out_file)  # Replace with your actual processing logic

    mean_keff.append(parse_file_same_line(out_file, begin='Mean Delta Keff (pcm):')[0])
    max_keff.append(parse_file_same_line(out_file, begin='Max  Delta Keff (pcm):')[0])
    mean_phi.append(parse_file_same_line(out_file, begin='Mean RMS Phi (%):')[0])
    max_phi.append(parse_file_same_line(out_file, begin='Max  RMS Phi (%):')[0])
    rom_dim.append(parse_file_same_line(out_file, begin='ROM_DIM:')[0])


## Plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.loglog(eps, mean_keff, 'o-', label='Mean Error')
ax.loglog(eps, max_keff, 'o-', label='Max Error')
ax.grid(True)
ax.legend()
ax.set_xlabel('$\epsilon_M$')
ax.set_ylabel('$\Delta K$eff Error (pcm)')
fig.savefig('keff_error_epsilonM.pdf', format='pdf')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.loglog(eps, mean_phi, 'o-', label='Mean RMS Error')
ax.loglog(eps, max_phi, 'o-', label='Max RMS Error')
ax.grid(True)
ax.legend()
ax.set_xlabel('$\epsilon_M$')
ax.set_ylabel('$\phi$ Error (%)')
fig.savefig('phi_error_epsilonM.pdf', format='pdf')

