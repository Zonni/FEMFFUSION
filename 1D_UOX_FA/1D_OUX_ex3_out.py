#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""



"""
import sys
sys.path.insert(1, '../postprocess/')
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
problem = '1D_UOX_FA_ex3'

file_1_fd =  '1D_UOX_FA_ex3_sp1_fd.out'
file_2_fd=  '1D_UOX_FA_ex3_sp3_fd.out'

files_fd = [file_1_fd, file_2_fd]
labels_fd = ['SP1', 'SP3']
style_fd = ['-', '-', '-']

# Time-Domain
folder   = 'time_domain/'

looking_freq = 1.0
file_sta_1 = folder + '1D_UOX_FA_ex3_sp1_td.out'  
file_nos_1 = file_sta_1 + '.nos'  
file_sta_2 = folder + '1D_UOX_FA_ex3_sp3_td.out'  
file_nos_2 = file_sta_2 + '.nos'  

files_sta_td = [file_sta_1, file_sta_2]
files_nos_td = [file_nos_1, file_nos_2]
labels_td = ['SP1-TD', 'SP3-TD']
style_td = ['*', 'v', 's']
    
n_files_fd = len(files_fd)
n_files_td = len(files_nos_td)

assert(len(files_sta_td) == len(files_nos_td))
assert(len(files_sta_td) == len(labels_td))
assert(len(labels_fd) == n_files_fd)
ny= 1


#%% ---------------------------------------------------------------------------
# GET MESH

nx = 174
ny = 1

x = [0.08,  # water_strip
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #1
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #2
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #3
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #4
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #5
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #6
0.06430, 
0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 
0.16570, 0.16570,                                                       #7
0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 
0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
0.06430,
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #8
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #9
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #10
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #11
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #12
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #13
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #14
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #15
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #16
0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #17
0.08];  # water_strip

assert(len(x)== nx)
# assert(len(y)== ny)
line_idx = 53;

n_cells = nx * ny
x = np.cumsum(x)
x = np.insert(x, 0, 0.0) # Insert a 0 in the first position
x = (x[:-1] + x[1:]) / 2  #Compute the average of consecutive elements


#%% ---------------------------------------------------------------------------
# FREQUENCY DOMAIN

amp_g1_line_fd = []
amp_g2_line_fd = []
pha_g1_line_fd = []
pha_g2_line_fd = []
sta_g1_line_fd = []
sta_g2_line_fd = []

for i in range(len(files_fd)):
    print('FD', files_fd[i], '...')
    # Get From .OUT
    sta_g1_fd = parse_file(files_fd[i], 'Group 1 flux', n_max_lines=ny)
    sta_g2_fd = parse_file(files_fd[i], 'Group 2 flux', n_max_lines=ny)
    sta_g1_fd = np.array(sta_g1_fd)
    sta_g2_fd = np.array(sta_g2_fd)
    
    nos_g1_fd = parse_file_complex(files_fd[i], 'Flux Noise Group 1', n_max_lines=ny)
    nos_g2_fd = parse_file_complex(files_fd[i], 'Flux Noise Group 2', n_max_lines=ny)
    
    # Relative noise
    rel_g1_fd = nos_g1_fd / sta_g1_fd * 100 
    rel_g2_fd = nos_g2_fd / sta_g2_fd * 100  
    
    sta_g1_line_fd.append(sta_g1_fd)
    sta_g2_line_fd.append(sta_g2_fd)
    
    amp_g1_line_fd.append(np.abs(rel_g1_fd))
    amp_g2_line_fd.append(np.abs(rel_g2_fd))
    pha_g1_line_fd.append(np.angle(rel_g1_fd, deg=True))
    pha_g2_line_fd.append(np.angle(rel_g2_fd, deg=True))

#%% ---------------------------------------------------------------------------
#  TIME DOMAIN

amp_g1_line_td = []
amp_g2_line_td = []
pha_g1_line_td = []
pha_g2_line_td = []
sta_g1_line_td = []
sta_g2_line_td = []

for i in range(n_files_td):
    print('Time Domain', files_nos_td[i], '...')

    
    sta_g1 = parse_file(files_sta_td[i], 'Group 1 flux', n_max_lines=ny)
    sta_g2 = parse_file(files_sta_td[i], 'Group 2 flux', n_max_lines=ny)
    sta_g1 = np.array(sta_g1)
    sta_g2 = np.array(sta_g2)
        
    time = parse_file(files_sta_td[i], begin='Time vector', n_max_lines=1)
    time = time[:-1]
    n_steps = len(time)
    steps = range(0, n_steps)
    
    noise_g1_td= np.zeros([n_steps, n_cells])
    noise_g2_td= np.zeros([n_steps, n_cells])
    for st in steps:
        noise_g1_td[st] = parse_file(files_nos_td[i],
                                      'Noise of group 1 time step ' + str(st) ,
                                      n_max_lines=ny)
        noise_g2_td[st] = parse_file(files_nos_td[i],
                                      'Noise of group 2 time step ' + str(st),
                                       n_max_lines=ny)
    
    
    
    # Transpose and normalize
    noise_g1_td= np.transpose(noise_g1_td)
    noise_g2_td= np.transpose(noise_g2_td) 
    
    freq   = np.fft.rfftfreq(n_steps, d=time[1])
    fft_g1 = np.fft.rfft(noise_g1_td) * 2.0/ n_steps 
    fft_g2 = np.fft.rfft(noise_g2_td) * 2.0/ n_steps 
    
    
    # We cut at looking_freq Hz
    cut_freq = int (looking_freq * n_steps * time[1])
    assert(freq[cut_freq] == looking_freq)
    noise_g1 = np.zeros([n_cells], dtype='cfloat')
    noise_g2 = np.zeros([n_cells], dtype='cfloat')
    for node in range(n_cells):
        noise_g1[node] = fft_g1[node][cut_freq] 
        noise_g2[node] = fft_g2[node][cut_freq]
    
    amp_g1_td= np.abs(noise_g1) / sta_g1 * 100
    amp_g2_td= np.abs(noise_g2) / sta_g2 * 100

    pha_g1_td= np.angle(noise_g1, deg=True)  
    pha_g2_td= np.angle(noise_g2, deg=True)
    
    sta_g1_line_td.append(sta_g1)
    sta_g2_line_td.append(sta_g2)
    amp_g1_line_td.append(amp_g1_td)
    amp_g2_line_td.append(amp_g2_td)
    pha_g1_line_td.append(pha_g1_td)
    pha_g2_line_td.append(pha_g2_td)

#%% ---------------------------------------------------------------------------

# Static g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files_fd):
    ax1.plot(x, sta_g1_line_fd[i], style_fd[i], label=labels_fd[i], color=colors[i])
    ax1.plot(x, sta_g1_line_td[i], style_td[i], label=labels_td[i], color=colors[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Static Flux g1")
fig1.savefig(folder + problem + "_static_g1.pdf", format='pdf')


# Static g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files_fd):
    ax1.plot(x, sta_g2_line_fd[i], style_fd[i], label=labels_fd[i], color=colors[i])
    ax1.plot(x, sta_g2_line_td[i], style_td[i], label=labels_td[i], color=colors[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Static Flux g2")
fig1.savefig(folder + problem + "_static_g2.pdf", format='pdf')


# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files_fd):
    ax1.plot(x, amp_g1_line_fd[i], style_fd[i], label=labels_fd[i], color=colors[i])
    ax1.plot(x, amp_g1_line_td[i], style_td[i], label=labels_td[i], color=colors[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel(r"Relative Noise Magnitude (\%)")
fig1.savefig(folder + problem + "_noise_line_amp_g1.pdf", format='pdf')


# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files_fd):
    ax1.plot(x, amp_g2_line_fd[i], style_fd[i], label=labels_fd[i], color=colors[i])
    ax1.plot(x, amp_g2_line_td[i], style_td[i], label=labels_td[i], color=colors[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel(r"Relative Noise Magnitude (\%)")
fig1.savefig(folder + problem + "_noise_line_amp_g2.pdf", format='pdf')


# Print phase_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files_fd):
    ax1.plot(x, pha_g1_line_fd[i], style_fd[i], label=labels_fd[i], color=colors[i])
    ax1.plot(x, pha_g1_line_td[i], style_td[i], label=labels_td[i], color=colors[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Phase (deg)")
fig1.savefig(folder + problem + "_noise_line_pha_g1.pdf", format='pdf')


# Print phase_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files_fd):
    ax1.plot(x, pha_g2_line_fd[i], style_fd[i], label=labels_fd[i], color=colors[i])
    ax1.plot(x, pha_g2_line_td[i], style_td[i], label=labels_td[i], color=colors[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Phase (deg)")
fig1.savefig(folder + problem + "_noise_line_pha_g2.pdf", format='pdf')

