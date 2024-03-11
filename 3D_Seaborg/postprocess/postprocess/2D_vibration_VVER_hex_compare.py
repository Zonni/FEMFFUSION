#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 13:03:14 2020

@author: zonni
"""

from utils import parse_vtk_file, parse_vtk_grid
import matplotlib.pyplot as plt
import matplotlib.colors
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
folder   = '../2D_VVER1000_vibration/'
file_bpr = folder + '2D_VVER1000.out.vtk'

file_ftd= folder + 'time_domain/2D_VVER1000_FE2.out'  
static_file = folder + 'time_domain/2D_VVER1000_FE2.out.vtk'  

looking_freq = 1.0
problem = '2D_test_hex'

cmap="viridis" 
#cmap="inferno"
#cmap="RdBu_r"

#%% ---------------------------------------------------------------------------
# Frequency domain boundary perturbation  data
# Get From VTK
[x, y, z] = parse_vtk_grid(file_bpr)
stati_g1 = parse_vtk_file(file_bpr, "Static_Flux_g1")
stati_g2 = parse_vtk_file(file_bpr, "Static_Flux_g2")
noise_g1 = parse_vtk_file(file_bpr, "Noise_g1_Magnitude")
noise_g2 = parse_vtk_file(file_bpr, "Noise_g2_Magnitude")
phase_g1 = parse_vtk_file(file_bpr, "Noise_g1_Phase")
phase_g2 = parse_vtk_file(file_bpr, "Noise_g2_Phase")

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
x_line_bpr = []
noise_g1_line_bpr = []
noise_g2_line_bpr = []
phase_g1_line_bpr = []
phase_g2_line_bpr = []

for p in range(len(y)):
    if (abs(y[p] - mid_y) < 1e-4):
        x_line_bpr.append(x[p])
        noise_g1_line_bpr.append(noise_g1[p])
        noise_g2_line_bpr.append(noise_g2[p])
        phase_g1_line_bpr.append(phase_g1[p])
        phase_g2_line_bpr.append(phase_g2[p])


#%% ---------------------------------------------------------------------------
## FEMFFUSION TIME DOMAIN
#from utils import parse_file
#import numpy as np
#
#
#steps_fem = range(0, 300)
#files =[]
#for i, step in enumerate(steps_fem):
#    files.append(file_ftd+ str(step) +  '.vtk')
#    
#time_fem = parse_file(file_ftd, begin='Time vector', n_max_lines=1)
#if (len(time_fem) == (len(steps_fem) + 1)):
#    time_fem = time_fem[:-1]
#n_nodes_fem = len(parse_vtk_grid(files[0])[0])
#n_steps_fem = len(steps_fem)
#delta_t_fem = time_fem[1]
#
#noise_fs = np.zeros([n_nodes_fem, n_steps_fem])
#noise_th = np.zeros([n_nodes_fem, n_steps_fem])
#
#for t, f in enumerate(files):
#    ffs = parse_vtk_file(f, 'noise_g1')
#    fth = parse_vtk_file(f, 'noise_g2')
#    for node in range(n_nodes_fem):   
#        noise_fs[node][t] = ffs[node] 
#        noise_th[node][t] = fth[node] 
#
#for t, f in enumerate(files):
#    if (t==0.0):
#        continue
#    for node in range(n_nodes_fem):   
#        noise_fs[node][t] -=  noise_fs[node][1]
#        noise_th[node][t] -=  noise_th[node][1]
#
## Compute Fast Fourier Transform
#freq_fem = np.fft.rfftfreq(n_steps_fem, d=delta_t_fem)
#noise_fft_fs = np.fft.rfft(noise_fs) * 2.0 / n_steps_fem
#noise_fft_th = np.fft.rfft(noise_th) * 2.0 / n_steps_fem
#
## We cut at looking_freq Hz
#cut_freq = int (looking_freq * n_steps_fem * delta_t_fem)
#assert(abs(freq_fem[cut_freq] - looking_freq) < 1e-4)
#
## Get noise amplitude and phase
#noise_g1 = np.zeros([n_nodes_fem])
#noise_g2 = np.zeros([n_nodes_fem])
#phase_g1 = np.zeros([n_nodes_fem])
#phase_g2 = np.zeros([n_nodes_fem])
#for node in range(n_nodes_fem):
#    noise_g1[node] = abs(noise_fft_fs[node][cut_freq])
#    noise_g2[node] = abs(noise_fft_th[node][cut_freq])
#    phase_g1[node] = np.angle(noise_fft_fs[node][cut_freq], deg=True)
#    phase_g2[node] = np.angle(noise_fft_th[node][cut_freq], deg=True)
#        
## Get Grid 
#[x, y, z] = parse_vtk_grid(static_file)
## Static Fluxes
#stati_g1 = parse_vtk_file(static_file, 'phi_g1_eig_1')
#stati_g2 = parse_vtk_file(static_file, 'phi_g2_eig_1')
#
## Remove repeated data
#stati_g1 = remove_repeated_data_point(x,y,z, stati_g1)
#stati_g2 = remove_repeated_data_point(x,y,z, stati_g2)
#noise_g1 = remove_repeated_data_point(x,y,z, noise_g1)
#noise_g2 = remove_repeated_data_point(x,y,z, noise_g2)
#phase_g1 = remove_repeated_data_point(x,y,z, phase_g1)
#phase_g2 = remove_repeated_data_point(x,y,z, phase_g2)
#x, y, z  = remove_repeated_point(x,y,z)
#
## Relative noise
#noise_g1 = noise_g1 / stati_g1 * 100
#noise_g2 = noise_g2 / stati_g2 * 100

#%% ---------------------------------------------------------------------------
# Midline
mid_y = max(y) / 2
x_line_ftd = []
noise_g1_line_ftd = []
noise_g2_line_ftd = []
phase_g1_line_ftd = []
phase_g2_line_ftd = []

for p in range(len(y)):
    if (abs(y[p] - mid_y) < 1e-4):
        x_line_ftd.append(x[p])
        noise_g1_line_ftd.append(noise_g1[p])
        noise_g2_line_ftd.append(noise_g2[p])
        phase_g1_line_ftd.append(phase_g1[p])
        phase_g2_line_ftd.append(phase_g2[p])
  
(x_line_bpr2,
 noise_g1_line_bpr2,
 noise_g2_line_bpr2,
 phase_g1_line_bpr2,
 phase_g2_line_bpr2) = zip(*sorted(zip(x_line_bpr, 
                               noise_g1_line_bpr,
                               noise_g2_line_bpr,
                               phase_g1_line_bpr,
                               phase_g2_line_bpr)))
#(x_line_ftd2,
# noise_g1_line_ftd2,
# noise_g2_line_ftd2,
# phase_g1_line_ftd2,
# phase_g2_line_ftd2) = zip(*sorted(zip(x_line_ftd, 
#                               noise_g1_line_ftd,
#                               noise_g2_line_ftd,
#                               phase_g1_line_ftd,
#                               phase_g2_line_ftd)))  
##
#%% ---------------------------------------------------------------------------

## Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_line_ftd2, noise_g1_line_ftd2, '-', label='Time-Domain')
ax1.plot(x_line_bpr2, noise_g1_line_bpr2, 'X', label='Edge wise FD')
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Relative Noise Magnitude (\%)")
fig1.savefig(problem + "_noise_line_amp_g1.pdf", format='pdf')

# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_line_ftd2, noise_g2_line_ftd2, '-', label='Time-Domain')
ax1.plot(x_line_bpr2, noise_g2_line_bpr2, 'X', label='Edge wise FD')
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Relative Noise Magnitude (\%)")
fig1.savefig(problem + "_noise_line_amp_g2.pdf", format='pdf')

# Print phase_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_line_ftd2, phase_g1_line_ftd2, '-', label='Time-Domain')
ax1.plot(x_line_bpr2, phase_g1_line_bpr2, 'X', label='Edge wise FD')
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Phase (deg)")
fig1.savefig(problem + "_noise_line_pha_g1.pdf", format='pdf')

# Print phase_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_line_ftd2, phase_g2_line_ftd2, '-', label='Time-Domain')
ax1.plot(x_line_bpr2, phase_g2_line_bpr2, 'X', label='Edge wise FD')
ax1.grid(True)
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Phase (deg)")
ax1.legend(loc='best')
fig1.savefig(problem + "_noise_line_pha_g2.pdf", format='pdf')

