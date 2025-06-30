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

snap_values = [2, 3, 5, 10, 15, 20, 25, 35, 40, 45, 50, 60, 70, 80, 90, 100,
               125, 150, 175, 200, 300, 400, 500]

max_keff = []
mean_keff = []
mean_phi = []
max_phi = []

for snap in snap_values:
    out_file = 'POD_n_snapshots/3D_Langenbuch_POD' + str(snap) + '_group_wise.out'
    n_snapshots = parse_file_same_line(out_file, begin='N_Snapshots:')[0]
    assert(n_snapshots == snap)
    
    mean_keff.append(parse_file_same_line(out_file, begin='Mean Delta Keff (pcm):')[0])
    max_keff.append(parse_file_same_line(out_file, begin='Max  Delta Keff (pcm):')[0])
    mean_phi.append(parse_file_same_line(out_file, begin='Mean RMS Phi (%):')[0])
    max_phi.append(parse_file_same_line(out_file, begin='Max  RMS Phi (%):')[0])

## Plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.semilogy(snap_values, mean_phi, 'o-', label='Mean Error')
ax.semilogy(snap_values, max_keff, 'o-', label='Max Error')
ax.grid(True)
ax.legend()
ax.set_xlabel('Number of Snapshots')
ax.set_ylabel('$\Delta K$eff Error (pcm)')
fig.savefig('keff_error.pdf', format='pdf')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.semilogy(snap_values, mean_keff, 'o-', label='Mean RMS Error')
ax.semilogy(snap_values, max_phi, 'o-', label='Max RMS Error')
ax.grid(True)
ax.legend()
ax.set_xlabel('Number of Snapshots')
ax.set_ylabel('$\phi$ Error (%)')
fig.savefig('phi_error.pdf', format='pdf')

#%% ===========================================================================
# LUPOD POINTS 

LUPOD_points=[]
for i in range(1, 51):
    LUPOD_points.append(int(10773 * 0.02 * i))  # Varying points per iteration
LUPOD_points = np.array(LUPOD_points)

points_percent = LUPOD_points/10773 * 100
max_keff = []
mean_keff = []
mean_phi = []
max_phi = []

for points in LUPOD_points:
    out_file = 'LUPOD_points/3D_Langenbuch_LUPODext' + str(points) + '_group_wise.out'
    n_points = parse_file_same_line(out_file, begin='N_LUPOD_Points:')[0]
    assert(points == n_points)
    
    mean_keff.append(parse_file_same_line(out_file, begin='Mean Delta Keff (pcm):')[0])
    max_keff.append(parse_file_same_line(out_file, begin='Max  Delta Keff (pcm):')[0])
    mean_phi.append(parse_file_same_line(out_file, begin='Mean RMS Phi (%):')[0])
    max_phi.append(parse_file_same_line(out_file, begin='Max  RMS Phi (%):')[0])

## Plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.semilogy(points_percent, mean_keff, 'o-', label='Mean Error')
ax.semilogy(points_percent, max_keff, 'o-', label='Max Error')
ax.grid(True)
ax.legend()
ax.set_xlabel('% of Collocation Points')
ax.set_ylabel('$\Delta K$eff Error (pcm)')
fig.savefig('keff_error_LUPODPoints.pdf', format='pdf')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.semilogy(points_percent, mean_phi, 'o-', label='Mean RMS Error')
ax.semilogy(points_percent, max_phi, 'o-', label='Max RMS Error')
ax.grid(True)
ax.legend()
ax.set_xlabel('% of Collocation Points')
ax.set_ylabel('$\phi$ Error (%)')
ax.set_ylim([1e-1,1e1])
fig.savefig('phi_error_LUPODPoints.pdf', format='pdf')

#%% ===========================================================================
# RANDOM POINTS 
random_points=[]
for i in range(1, 51):
    random_points.append(int(10773 * 0.02 * i))  # Varying points per iteration
random_points = np.array(random_points)

points_percent = random_points/10773 * 100
max_keff_rand = []
mean_keff_rand = []
mean_phi_rand = []
max_phi_rand = []

for points in LUPOD_points:
    print(points)
    out_file = 'random_npoints/3D_Langenbuch_random' + str(points) + '_group_wise.out'
    n_points = parse_file_same_line(out_file, begin='N_LUPOD_Points:')[0]
    assert(points == n_points)
    
    mean_keff_rand.append(parse_file_same_line(out_file, begin='Mean Delta Keff (pcm):')[0])
    max_keff_rand.append(parse_file_same_line(out_file, begin='Max  Delta Keff (pcm):')[0])
    mean_phi_rand.append(parse_file_same_line(out_file, begin='Mean RMS Phi (%):')[0])
    max_phi_rand.append(parse_file_same_line(out_file, begin='Max  RMS Phi (%):')[0])

## Plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.semilogy(points_percent, mean_keff_rand, 'o-', label='Mean Error Random')
ax.semilogy(points_percent, max_keff_rand, 'o-', label='Max Error Random')
ax.semilogy(points_percent, mean_keff, 'o-', label='Mean Error LUPODext')
ax.semilogy(points_percent, max_keff, 'o-', label='Max Error LUPODext')
ax.grid(True)
ax.legend()
ax.set_xlabel('% of Collocation Points')
ax.set_ylabel('$\Delta K$eff Error (pcm)')
fig.savefig('keff_error_randomPoints.pdf', format='pdf')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.semilogy(points_percent, mean_phi_rand, 'o-', label='Mean RMS Error Random')
ax.semilogy(points_percent, max_phi_rand, 'o-', label='Max RMS Error Random')
ax.semilogy(points_percent, mean_phi, 'o-', label='Mean RMS Error LUPODext')
ax.semilogy(points_percent, max_phi, 'o-', label='Max RMS Error LUPODext')
ax.grid(True)
ax.legend()
ax.set_xlabel('% of Collocation Points')
ax.set_ylabel('$\phi$ Error (%)')
ax.set_ylim([1e-1,1e1])
fig.savefig('phi_error_randomPoints.pdf', format='pdf')

