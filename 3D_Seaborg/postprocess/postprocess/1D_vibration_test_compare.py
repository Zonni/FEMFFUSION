#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 13:03:14 2020

@author: zonni
"""

from utils import parse_vtk_file, parse_vtk_grid
import matplotlib.pyplot as plt
#import matplotlib.colors
from matplotlib import rcParams
from utils import remove_repeated_point, remove_repeated_data_point

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

problem = '1D_vib_test_SPN'
folder   = '../test/1D_vib_test_SPN/'

file_1 = folder + '1D_vib_diff.out.vtk'
file_2 = folder + '1D_vib_sp3.out.vtk'
file_3 = folder + '1D_vib_fsp3.out.vtk'
files = [file_1, file_2, file_3]
labels = ['Diffusion', 'SP3', 'FSP3']
style = ['*-', '.--', '+']
n_files = len(files)

assert(len(files)==len(labels) and len(labels)==n_files)
 
cmap="viridis" 

#cmap="inferno"
#cmap="RdBu_r"

#%% ---------------------------------------------------------------------------
# Frequency domain boundary perturbation data
x_line = []
noise_g1_line = []
noise_g2_line = []
phase_g1_line = []
phase_g2_line = []
static_g1_line= []
static_g2_line= []

for i in range(n_files):
    print(files[i] + '...')
    # Get From VTK
    [x, y, z] = parse_vtk_grid(files[i])
    stati_g1 = parse_vtk_file(files[i], "Static_Flux_g1")
    stati_g2 = parse_vtk_file(files[i], "Static_Flux_g2")
    noise_g1 = parse_vtk_file(files[i], "Noise_g1_Magnitude")
    noise_g2 = parse_vtk_file(files[i], "Noise_g2_Magnitude")
    phase_g1 = parse_vtk_file(files[i], "Noise_g1_Phase")
    phase_g2 = parse_vtk_file(files[i], "Noise_g2_Phase")
    
    # Remove repeated data
    stati_g1 = remove_repeated_data_point(x,y,z, stati_g1)
    stati_g2 = remove_repeated_data_point(x,y,z, stati_g2)
    noise_g1 = remove_repeated_data_point(x,y,z, noise_g1)
    noise_g2 = remove_repeated_data_point(x,y,z, noise_g2)
    phase_g1 = remove_repeated_data_point(x,y,z, phase_g1)
    phase_g2 = remove_repeated_data_point(x,y,z, phase_g2)
    x,y,z = remove_repeated_point(x,y,z)
    
    # Relative noise
    noise_g1 = noise_g1 / stati_g1 * 100
    noise_g2 = noise_g2 / stati_g2 * 100
    
    # Midline
    mid_y = max(y) / 2
    x_line.append([])
    noise_g1_line.append([])
    noise_g2_line.append([])
    phase_g1_line.append([])
    phase_g2_line.append([])
    static_g1_line.append([])
    static_g2_line.append([])
    
    for p in range(len(y)):
        if (y[p] == mid_y):
            x_line[i].append(x[p])
            noise_g1_line[i].append(noise_g1[p])
            noise_g2_line[i].append(noise_g2[p])
            phase_g1_line[i].append(phase_g1[p])
            phase_g2_line[i].append(phase_g2[p])
            static_g1_line[i].append(stati_g1[p])
            static_g2_line[i].append(stati_g2[p])



#%% ---------------------------------------------------------------------------

# Static g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_line_ftd, noise_g1_line_ftd, '-', label='Time-Domain')
for i in range(n_files):
    ax1.plot(x_line[i], static_g1_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Static Flux g1")
fig1.savefig(folder + problem + "_static_g1.pdf", format='pdf')


# Static g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_line_ftd, noise_g1_line_ftd, '-', label='Time-Domain')
for i in range(n_files):
    ax1.plot(x_line[i], static_g2_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Static Flux g2")
fig1.savefig(folder + problem + "_static_g2.pdf", format='pdf')

# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line[i], noise_g1_line[i], style[i], label=labels[i])
#ax1.plot(x_line_ftd, noise_g1_line_ftd, '-', label='Time-Domain')
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Relative Noise Magnitude (\%)")
fig1.savefig(folder + problem + "_noise_line_amp_g1.pdf", format='pdf')

# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_line_ftd, noise_g2_line_ftd, '-', label='Time-Domain')
for i in range(n_files):
    ax1.plot(x_line[i], noise_g2_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Relative Noise Magnitude (\%)")
fig1.savefig(folder + problem + "_noise_line_pha_g1.pdf", format='pdf')


# Print phase_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_line_ftd, noise_g2_line_ftd, '-', label='Time-Domain')
for i in range(n_files):
    ax1.plot(x_line[i], phase_g1_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Phase (deg)")
fig1.savefig(folder + problem + "_noise_line_pha_g1.pdf", format='pdf')


# Print phase_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_line_ftd, noise_g2_line_ftd, '-', label='Time-Domain')
for i in range(n_files):
    ax1.plot(x_line[i], phase_g2_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Phase (deg)")
fig1.savefig(folder + problem + "_noise_line_amp_g2.pdf", format='pdf')


