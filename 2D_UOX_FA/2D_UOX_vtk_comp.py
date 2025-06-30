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
#import matplotlib.colors
from matplotlib import rcParams
from utils import remove_repeated_point, remove_repeated_data_point
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

problem = '2D_UOX_FA_vtk'
# folder   = 'exercise_2/'
folder   = 'exercise_3/'
file_1 = folder + '2D_UOX_FA_diff.out.vtk'
file_2 = folder + '2D_UOX_FA_sp3.out.vtk'
file_3 = folder + '2D_UOX_FA_sp5.out.vtk'

#files = [file_1, file_2, file_3]
#labels = ['Diffusion', 'SP3', 'FSP3']
#style = ['*-', '.--', '+']

files = [file_1, file_2, file_3]
labels = ['Diffusion', 'SP3', 'SP5']
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
    # Get From VTK
    print('Postprocessing file ', files[i] )
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
    
    # Diagonal
#    mid_y = max(y) / 2
    x_line.append([])
    noise_g1_line.append([])
    noise_g2_line.append([])
    phase_g1_line.append([])
    phase_g2_line.append([])
    static_g1_line.append([])
    static_g2_line.append([])
    
    for p in range(len(y)):
        if (y[p] == x[p]):
            x_line[i].append(x[p] * np.sqrt(2)) # Diagonal lenght
            noise_g1_line[i].append(noise_g1[p])
            noise_g2_line[i].append(noise_g2[p])
            phase_g1_line[i].append(phase_g1[p])
            phase_g2_line[i].append(phase_g2[p])
            static_g1_line[i].append(stati_g1[p])
            static_g2_line[i].append(stati_g2[p])
    




      
#%% ---------------------------------------------------------------------------
## FEMFFUSION TIME DOMAIN
#from utils import parse_file
#import numpy as np
#teponer r a la cadena evita que Python interprete las s
#file_out = folder + 'timedomain/2D_test_ref.out'  311s
#static_file = folder + 'timedomain/2D_test_ref.out.vtk'  
#steps_fem = range(0, 300)
#files =[]
#for i, step in enumerate(steps_fem):
#    files.append(file_out+ str(step) +  '.vtk')
#    
#time_fem = parse_file(file_out, begin='Time vector', n_max_lines=1)
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
#
## Midline
#mid_y = max(y) / 2
#x_line_ftd = []
#noise_g1_line_ftd = []
#noise_g2_line_ftd = []
#phase_g1_line_ftd = []
#phase_g2_line_ftd = []
#teponer r a la cadena evita que Python interprete las s
#for p in range(len(y)):
#    if (y[p] == mid_y):
#        x_line_ftd.append(x[p])
#        noise_g1_line_ftd.append(noise_g1[p])
#        noise_g2_line_ftd.append(noise_g2[p])
#        phase_g1_line_ftd.append(phase_g1[p])
#        phase_g2_line_ftd.append(phase_g2[p])
#        

#%% ---------------------------------------------------------------------------

# Static g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_line_ftd, noise_g1_line_ftd, '-', label='Time-Domain')
for i in range(n_files):
    ax1.plot(x_line[i], static_g1_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
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
ax1.set_xlabel("Diagonal length (cm)")
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
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel(r"Relative Noise Magnitude (\%)")
fig1.savefig(folder + problem + "_noise_line_amp_g1.pdf", format='pdf')

# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_line_ftd, noise_g2_line_ftd, '-', label='Time-Domain')
for i in range(n_files):
    ax1.plot(x_line[i], noise_g2_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel(r"Relative Noise Magnitude (\%)")
fig1.savefig(folder + problem + "_noise_line_pha_g1.pdf", format='pdf')


# Print phase_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_line_ftd, noise_g2_line_ftd, '-', label='Time-Domain')
for i in range(n_files):
    ax1.plot(x_line[i], phase_g1_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
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
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Phase (deg)")
fig1.savefig(folder + problem + "_noise_line_amp_g2.pdf", format='pdf')


