#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 13:03:14 2020

@author: zonni
"""
import sys
sys.path.insert(1, '../postprocess/')
from utils import parse_vtk_file, parse_vtk_grid
import matplotlib.pyplot as plt
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

problem = '2D_VVER1000'

line_y_plot = 81.753

file_1 =  '2D_VVER1000_SP1.out.vtk'
file_2 =  '2D_VVER1000_SP3.out.vtk'
file_3 =  '2D_VVER1000_SP5.out.vtk'
file_4 =  '2D_VVER1000_SP7.out.vtk'

#
##files = [file_1, file_2, file_3]
##labels = ['Diffusion', 'SP3', 'FSP3']
##style = ['*-', '.--', '+']
#
files = [file_1, file_2, file_3, file_4]
labels = ['SP1', 'SP3', 'SP5', 'SP7']
style = ['*-', '.', '^', '+']
n_files = len(files)
assert(len(files)==len(labels) and len(labels)==n_files)

#%% ---------------------------------------------------------------------------
# Frequency domain boundary perturbation  data
# Get From VTK
x_line = []
noise_g1_line = []
noise_g2_line = []
phase_g1_line = []
phase_g2_line = []
stati_g1_line = []
stati_g2_line = []

for i in range(n_files):
    print(files[i], '...')

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
    x_line_bpr = []
    stati_g1_line_bpr = []
    stati_g2_line_bpr = []
    noise_g1_line_bpr = []
    noise_g2_line_bpr = []
    phase_g1_line_bpr = []
    phase_g2_line_bpr = []
    
    for p in range(len(y)):
        if (abs(y[p] - line_y_plot) < 1e-4):
            x_line_bpr.append(x[p])
            stati_g1_line_bpr.append(stati_g1[p])
            stati_g2_line_bpr.append(stati_g2[p])
            noise_g1_line_bpr.append(noise_g1[p])
            noise_g2_line_bpr.append(noise_g2[p])
            phase_g1_line_bpr.append(phase_g1[p])
            phase_g2_line_bpr.append(phase_g2[p])
            

    # SORT        
    (x_line_bpr2,
     stati_g1_line_bpr2,
     stati_g2_line_bpr2,
     noise_g1_line_bpr2,
     noise_g2_line_bpr2,
     phase_g1_line_bpr2,
     phase_g2_line_bpr2) = zip(*sorted(zip(x_line_bpr, 
                                   stati_g1_line_bpr,
                                   stati_g2_line_bpr,
                                   noise_g1_line_bpr,
                                   noise_g2_line_bpr,
                                   phase_g1_line_bpr,
                                   phase_g2_line_bpr)))
    x_line.append(x_line_bpr2)
    stati_g1_line.append(stati_g1_line_bpr2)
    stati_g2_line.append(stati_g2_line_bpr2)
    noise_g1_line.append(noise_g1_line_bpr2)
    noise_g2_line.append(noise_g2_line_bpr2)
    phase_g1_line.append(phase_g1_line_bpr2)
    phase_g2_line.append(phase_g2_line_bpr2)

    #
#%% ---------------------------------------------------------------------------

## Print stati_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line[i], stati_g1_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Static Flux (AU)")
fig1.savefig(problem + "_stati_line_g1.pdf", format='pdf')

## Print stati_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line[i], stati_g2_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Static Flux (AU)")
fig1.savefig(problem + "_stati_line_g2.pdf", format='pdf')

# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line[i], noise_g1_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Relative Noise Magnitude (\%)")
fig1.savefig(problem + "_noise_line_amp_g1.pdf", format='pdf')

# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line[i], noise_g2_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Relative Noise Magnitude (\%)")
fig1.savefig(problem + "_noise_line_amp_g2.pdf", format='pdf')

# Print phase_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line[i], phase_g1_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Phase (deg)")
fig1.savefig(problem +"_noise_line_pha_g1.pdf", format='pdf')

# Print phase_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line[i], phase_g2_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Phase (deg)")
ax1.legend(loc='best')
fig1.savefig(problem + "_noise_line_pha_g2.pdf", format='pdf')

