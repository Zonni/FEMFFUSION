# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:33:28 2019

@author: zonni
"""
from utils import parse_vtk_file, parse_vtk_grid
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib import rcParams
import scipy.io
import numpy as np 

plt.close('all')

params = {'backend': 'pdf',
#          'font.family': 'serif',
          'font.size': 14,
          'axes.labelsize': 12,
          'legend.fontsize': 12,
          'xtick.labelsize': 13,
          'ytick.labelsize': 13,
          'text.usetex': True,
          'lines.linewidth': 2,
          'lines.markersize': 5,
          'lines.markeredgewidth': 1,
          'legend.numpoints': 1, 
          'axes.formatter.useoffset': False,
          'figure.autolayout': True,
          }

rcParams.update(params)

file = '../2D_test_hex/2D_test.out.vtk'
file_fft = '../2D_test_hex/time_2D_test_hex_f1.mat'
problem = '2D_test_hex'
#cmap="viridis" 
#cmap="inferno"
cmap="RdBu_r"

# Get From VTK
[x, y, z] = parse_vtk_grid(file)
static_g1 = parse_vtk_file(file, "Static_Flux_g1")
static_g2 = parse_vtk_file(file, "Static_Flux_g2")
noise_g1 = parse_vtk_file(file, "Noise_g1_Magnitude")
noise_g2 = parse_vtk_file(file, "Noise_g2_Magnitude")
phase_g1 = parse_vtk_file(file, "Noise_g1_Phase")
phase_g2 = parse_vtk_file(file, "Noise_g2_Phase")
n_levels = 10




mat = scipy.io.loadmat(file_fft)

x_py = mat['x_py'][0]
y_py = mat['y_py'][0]
py_fs = np.absolute(mat['py_fs'][0])
py_th = np.absolute(mat['py_th'][0])
py_pha_fs = np.angle(mat['py_fs'][0], deg=True)
py_pha_th = np.angle(mat['py_th'][0], deg=True)
diff_mag1 = abs(py_fs - noise_g1) /py_fs *100;
diff_mag2 = abs(py_th - noise_g2) /py_th *100;
diff_pha1 = abs(py_pha_fs - phase_g1) /py_pha_fs *100;
diff_pha2 = abs(py_pha_th - phase_g2) /py_pha_th *100;

## Print noise_g1
#fig1 = plt.figure()
#ax1 = fig1.add_subplot(1, 1, 1)
##ax1.tricontour(x, y, static_g1, levels=n_levels, linewidths=0.5, colors='k')
#cntr2 = ax1.tricontourf(x, y, static_g1, levels=n_levels, cmap=cmap)
#fig1.colorbar(cntr2, ax=ax1)
#ax1.set_ylabel("y (cm)")
#ax1.set_xlabel("x (cm)")
#ax1.set_title("Static Fast Flux")
#fig1.savefig(problem + "_static_g1.pdf", format='pdf')

# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.tricontour(x_py, y_py, diff_mag1, levels=n_levels, linewidths=0.5, colors='k')
cntr2 = ax1.tricontourf(x_py, y_py, diff_mag1, levels=n_levels, cmap=cmap)
fig1.colorbar(cntr2, ax=ax1)
ax1.set_ylabel("y (cm)")
ax1.set_xlabel("x (cm)")
ax1.set_title("Magnitude g1 Difference (\%)")


# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.tricontour(x_py, y_py, diff_mag2, levels=n_levels, linewidths=0.5, colors='k')
cntr2 = ax1.tricontourf(x_py, y_py, diff_mag2, levels=n_levels, cmap=cmap)
fig1.colorbar(cntr2, ax=ax1)
ax1.set_ylabel("y (cm)")
ax1.set_xlabel("x (cm)")
ax1.set_title("Magnitude g2 Difference (\%)")


# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.tricontour(x_py, y_py, diff_pha1, levels=n_levels, linewidths=0.5, colors='k')
cntr2 = ax1.tricontourf(x_py, y_py, diff_pha1, levels=n_levels, cmap=cmap)
fig1.colorbar(cntr2, ax=ax1)
ax1.set_ylabel("y (cm)")
ax1.set_xlabel("x (cm)")
ax1.set_title("Phase g1 Difference (\%)")


# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.tricontour(x_py, y_py, diff_pha2, levels=n_levels, linewidths=0.5, colors='k')
cntr2 = ax1.tricontourf(x_py, y_py, diff_pha2, levels=n_levels, cmap=cmap)
fig1.colorbar(cntr2, ax=ax1)
ax1.set_ylabel("y (cm)")
ax1.set_xlabel("x (cm)")
ax1.set_title("Phase G2 Difference (\%)")



