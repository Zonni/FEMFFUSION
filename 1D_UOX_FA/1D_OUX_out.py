#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 13:03:14 2020

@author: zonni
"""

import matplotlib.pyplot as plt
#import matplotlib.colors
from matplotlib import rcParams
from utils import parse_file, parse_file_complex
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
          'lines.markersize': 2.5,
          'lines.markeredgewidth': 1,
          'legend.numpoints': 1, 
          'axes.formatter.useoffset': False,
          'figure.autolayout': True,
          }

rcParams.update(params)

colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c',
          u'#d62728', u'#9467bd', u'#8c564b']

#%% FILES

# Frequency-Domain
#problem = '1D_UOX_FA_out'
problem = '1D_UOX_FA_sint'

if (problem == '1D_UOX_FA_out'):
    folder   = '../1D_UOX_FA/'
    file_1_fd = folder + '1D_UOX_FA_ex2_diff.out'
    file_2_fd= folder + '1D_UOX_FA_ex2_sp3.out'
    file_3_fd = folder + '1D_UOX_FA_ex2_fsp3.out'
    
    files_fd = [file_1_fd, file_2_fd, file_3_fd]
    labels_fd = ['DSP1', 'DSP3', 'FSP3']
    style_fd = ['-', '-', '-']
    
    # Time-Domain
    folder   = '../1D_UOX_FA/time_domain/'
    
    file_sta_1 = folder + '1D_UOX_FA_ex2_diffout'  
    file_nos_1 = folder + '1D_UOX_FA_ex2_diff.outnos'  
    file_sta_2 = folder + '1D_UOX_FA_ex2_sp3out'  
    file_nos_2 = folder + '1D_UOX_FA_ex2_sp3.outnos'  
    file_sta_3 = folder + '1D_UOX_FA_ex2_fsp3out'  
    file_nos_3 = folder + '1D_UOX_FA_ex2_fsp3.outnos'  
    
    files_sta_td = [file_sta_1, file_sta_2, file_sta_3]
    files_nos_td = [file_nos_1, file_nos_2, file_nos_3]
    labels_td = ['DSP1-TD', 'DSP3-TD', 'FSP3-TD']
    style_td = ['*', 'v', 's']
    

elif(problem == "1D_UOX_FA_sint"):

    folder   = '../1D_UOX_FA/'
    file_1_fd = folder + '1D_UOX_FA_ex2_diff.out'
    file_2_fd= folder + '1D_UOX_FA_ex2_sp3.out'
    file_3_fd = folder + '1D_UOX_FA_ex2_fsp3.out'
    
    files_fd = [file_1_fd, file_2_fd, file_3_fd]
    labels_fd = ['DSP1', 'DSP3', 'FSP3']
    style_fd = ['-', '-', '-']
    
    # Time-Domain
    folder   = '../1D_UOX_FA/time_domain/'
    
    file_sta_1 = folder + '1D_UOX_FA_ex2_diff_sintout'  
    file_nos_1 = folder + '1D_UOX_FA_ex2_diff_sint.outnos'  
    file_sta_2 = folder + '1D_UOX_FA_ex2_sp3_sintout'  
    file_nos_2 = folder + '1D_UOX_FA_ex2_sp3_sint.outnos'  
    file_sta_3 = folder + '1D_UOX_FA_ex2_fsp3out'  
    file_nos_3 = folder + '1D_UOX_FA_ex2_fsp3.outnos'  
    
    files_sta_td = [file_sta_1, file_sta_2, file_sta_3]
    files_nos_td = [file_nos_1, file_nos_2, file_nos_3]
    labels_td = ['DSP1-TD', 'DSP3-TD', 'FSP3-TD']
    style_td = ['*', 'v', 's']
    



n_files_fd = len(files_fd)
n_files_td = len(files_nos_td)

assert(len(files_sta_td) == len(files_nos_td))
assert(len(files_sta_td) == len(labels_td))
assert(len(labels_fd) == n_files_fd)
ny= 1


#%% ---------------------------------------------------------------------------
# Frequency domain boundary perturbation data
x_line = [0.08,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
          0.08]
x_line = np.cumsum(np.array(x_line))



noise_g1_line = []
noise_g2_line = []
phase_g1_line = []
phase_g2_line = []
static_g1_line= []
static_g2_line= []

for i in range(len(files_fd)):
    print('FD', files_fd[i], '...')
    # Get From .OUT
    sta_g1 = parse_file(files_fd[i], 'Group 1 flux', n_max_lines=ny)
    sta_g2 = parse_file(files_fd[i], 'Group 2 flux', n_max_lines=ny)
    sta_g1 = np.array(sta_g1)
    sta_g2 = np.array(sta_g2)
    
    nos_g1 = parse_file_complex(files_fd[i], 'Flux Noise Group 1', n_max_lines=ny)
    nos_g2 = parse_file_complex(files_fd[i], 'Flux Noise Group 2', n_max_lines=ny)
    
    # Relative noise
    relnos_g1 = nos_g1 / sta_g1 * 100 /100
    relnos_g2 = nos_g2 / sta_g2 * 100 /100 
    
    static_g1_line.append(sta_g1)
    static_g2_line.append(sta_g2)
    
    noise_g1_line.append(np.abs(relnos_g1))
    noise_g2_line.append(np.abs(relnos_g2))
    phase_g1_line.append(np.angle(relnos_g1, deg=True))
    phase_g2_line.append(np.angle(relnos_g2, deg=True))

#%% ---------------------------------------------------------------------------

noise_g1_line_td = []
noise_g2_line_td = []
phase_g1_line_td = []
phase_g2_line_td = []
statc_g1_line_td = []
statc_g2_line_td = []

for i in range(n_files_td):
    print('Time Domain', files_nos_td[i], '...')
    looking_freq = 1.0
    n_cells = 138
    
    sta_g1 = parse_file(files_sta_td[i], 'Group 1 flux', n_max_lines=ny)
    sta_g2 = parse_file(files_sta_td[i], 'Group 2 flux', n_max_lines=ny)
    sta_g1 = np.array(sta_g1)
    sta_g2 = np.array(sta_g2)
        
    time_fem = parse_file(files_sta_td[i], begin='Time vector', n_max_lines=1)
    time_fem = time_fem[:-1]
    n_steps = len(time_fem)
    steps = range(0, n_steps)
    
    noise_g1_fem = np.zeros([n_steps, n_cells])
    noise_g2_fem = np.zeros([n_steps, n_cells])
    for st in steps:
        noise_g1_fem[st] = parse_file(files_nos_td[i],
                                      'Noise of group 1 time step ' + str(st) ,
                                      n_max_lines=ny)
        noise_g2_fem[st] = parse_file(files_nos_td[i],
                                      'Noise of group 2 time step ' + str(st),
                                       n_max_lines=ny)
    
    
    
    # Transpose and normalize
    noise_g1_fem = np.transpose(noise_g1_fem)
    noise_g2_fem = np.transpose(noise_g2_fem) 
    
    freq   = np.fft.rfftfreq(n_steps, d=time_fem[1])
    fft_g1 = np.fft.rfft(noise_g1_fem) * 2.0/ n_steps 
    fft_g2 = np.fft.rfft(noise_g2_fem) * 2.0/ n_steps
    
    
    # We cut at looking_freq Hz
    cut_freq = int (looking_freq * n_steps * time_fem[1])
    assert(freq[cut_freq] == looking_freq)
    noise_g1 = np.zeros([n_cells], dtype='cfloat')
    noise_g2 = np.zeros([n_cells], dtype='cfloat')
    for node in range(n_cells):
        noise_g1[node] = fft_g1[node][cut_freq] 
        noise_g2[node] = fft_g2[node][cut_freq]
    
    amp_g1_fem = np.abs(noise_g1) / sta_g1 * 100
    amp_g2_fem = np.abs(noise_g2) / sta_g2 * 100
    # FEMFUSSION PERTURBATION is a SINE not a cosine
    # This way we convert it
    pha_g1_fem = np.angle(noise_g1, deg=True)  + 90
    pha_g2_fem = np.angle(noise_g2, deg=True)  + 90
    
    statc_g1_line_td.append(sta_g1)
    statc_g2_line_td.append(sta_g2)
    noise_g1_line_td.append(amp_g1_fem)
    noise_g2_line_td.append(amp_g2_fem)
    phase_g1_line_td.append(pha_g1_fem)
    phase_g2_line_td.append(pha_g2_fem)

#%% ---------------------------------------------------------------------------

# Static g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files_fd):
    ax1.plot(x_line, static_g1_line[i], style_fd[i], label=labels_fd[i], color=colors[i])
    ax1.plot(x_line, statc_g1_line_td[i], style_td[i], label=labels_td[i], color=colors[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Static Flux g1")
fig1.savefig(folder + problem + "_static_g1.pdf", format='pdf')


# Static g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files_fd):
    ax1.plot(x_line, static_g2_line[i], style_fd[i], label=labels_fd[i], color=colors[i])
    ax1.plot(x_line, statc_g2_line_td[i], style_td[i], label=labels_td[i], color=colors[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Static Flux g2")
fig1.savefig(folder + problem + "_static_g2.pdf", format='pdf')


# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files_fd):
    ax1.plot(x_line, noise_g1_line[i], style_fd[i], label=labels_fd[i], color=colors[i])
    ax1.plot(x_line, noise_g1_line_td[i], style_td[i], label=labels_td[i], color=colors[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Relative Noise Magnitude (\%)")
fig1.savefig(folder + problem + "_noise_line_amp_g1.pdf", format='pdf')


# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files_fd):
    ax1.plot(x_line, noise_g2_line[i], style_fd[i], label=labels_fd[i], color=colors[i])
    ax1.plot(x_line, noise_g2_line_td[i], style_td[i], label=labels_td[i], color=colors[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Relative Noise Magnitude (\%)")
fig1.savefig(folder + problem + "_noise_line_amp_g2.pdf", format='pdf')


# Print phase_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files_fd):
    ax1.plot(x_line, phase_g1_line[i], style_fd[i], label=labels_fd[i], color=colors[i])
    ax1.plot(x_line, phase_g1_line_td[i], style_td[i], label=labels_td[i], color=colors[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Phase (deg)")
fig1.savefig(folder + problem + "_noise_line_pha_g1.pdf", format='pdf')


# Print phase_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files_fd):
    ax1.plot(x_line, phase_g2_line[i], style_fd[i], label=labels_fd[i], color=colors[i])
    ax1.plot(x_line, phase_g2_line_td[i], style_td[i], label=labels_td[i], color=colors[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Phase (deg)")
fig1.savefig(folder + problem + "_noise_line_pha_g2.pdf", format='pdf')

