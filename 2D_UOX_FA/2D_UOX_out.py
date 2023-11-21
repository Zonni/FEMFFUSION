#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 13:03:14 2020

@author: zonni
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
          'lines.markersize': 5,
          'lines.markeredgewidth': 1,
          'legend.numpoints': 1, 
          'axes.formatter.useoffset': False,
          'figure.autolayout': True,
          }

rcParams.update(params)

problem = '2D_UOX_FA_ex2'
folder   = 'exercise_2/'
#
file_1 = folder + '2D_UOX_FA_sp1.out'
file_2 = folder + '2D_UOX_FA_sp3.out'
file_3 = folder + '2D_UOX_FA_sp5.out'
file_4 = folder + '2D_UOX_FA_sp7.out'

#
##files = [file_1, file_2, file_3]
##labels = ['Diffusion', 'SP3', 'FSP3']
##style = ['*-', '.--', '+']
#
files = [file_1, file_2, file_3, file_4]
labels = ['SP1', 'SP3', 'SP5', 'SP7']
style = ['*-', '.-', '-^', '--+']
n_files = len(files)
assert(len(files)==len(labels) and len(labels)==n_files)

# ====================================================================== #
# TIME DOMAIN
ny = 138
nx = 138
n_cells = nx * ny

folder_txt   = '../2D_UOX_FA/time_domain/'
file_txt_1 = folder_txt + '2D_UOX_diffusion_FE1'  
file_txt_2 = folder_txt + '2D_UOX_sp3_sint'  
#files_txt = [file_txt_1, file_txt_2]
#style_td = ['+C0', '*C1'] 
#labels_td = ['TD-Diff', 'TD-SP3']

files_txt = []
style_td = [] 
labels_td = []
n_files_txt = len(files_txt)


assert(len(files_txt)==n_files_txt and
       len(labels_td)==n_files_txt and 
       len(style_td)==n_files_txt)


cmap="viridis" 
#cmap="inferno"
#cmap="RdBu_r"


#%% ---------------------------------------------------------------------------
# Line Geometry
x = [0.08,
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
x = np.cumsum(np.array(x))
y = x.copy()
x_line = np.sqrt(2)*x

#%% ---------------------------------------------------------------------------
# GET FD DATA
noise_g1_line = []
noise_g2_line = []
phase_g1_line = []
phase_g2_line = []
static_g1_line= []
static_g2_line= []

for i in range(n_files):
    print(files[i], '...')
    # Get From .OUT
    sta_g1 = parse_file(files[i], 'Group 1 flux', n_max_lines=ny)
    sta_g2 = parse_file(files[i], 'Group 2 flux', n_max_lines=ny)
    sta_g1 = np.array(sta_g1)
    sta_g2 = np.array(sta_g2)
    
    nos_g1 = parse_file_complex(files[i], 'Flux Noise Group 1', n_max_lines=ny)
    nos_g2 = parse_file_complex(files[i], 'Flux Noise Group 2', n_max_lines=ny)
    
    # Relative noise
    relnos_g1 = nos_g1 / sta_g1 * 100 
    relnos_g2 = nos_g2 / sta_g2 * 100 
    
    # Get diagonal
    sta_line_g1 = np.zeros(ny)
    sta_line_g2 = np.zeros(ny)
    nos_line_g1 = np.zeros(ny, dtype=np.complex_)
    nos_line_g2 = np.zeros(ny, dtype=np.complex_)
    for i in range(ny):
        sta_line_g1[i] = sta_g1[i*(nx +1)]
        sta_line_g2[i] = sta_g2[i*(nx +1)]
        nos_line_g1[i] = relnos_g1[i*(nx +1)]
        nos_line_g2[i] = relnos_g2[i*(nx +1)]
    
    static_g1_line.append(sta_line_g1)
    static_g2_line.append(sta_line_g2)
    
    noise_g1_line.append(np.abs(nos_line_g1))
    noise_g2_line.append(np.abs(nos_line_g2))
    phase_g1_line.append(np.angle(nos_line_g1, deg=True))
    phase_g2_line.append(np.angle(nos_line_g2, deg=True))
  
#%% ---------------------------------------------------------------------------
# GET TXT DATA 
noise_g1_line_txt = []
noise_g2_line_txt = []
phase_g1_line_txt = []
phase_g2_line_txt = []
static_g1_line_txt = []
static_g2_line_txt = []


for i in range(n_files_txt):
    file_txt = files_txt[i]
    print(files_txt[i], '...')
    
    sta_g1_txt = np.loadtxt(file_txt + '_STATIC_FLX_G1.txt')
    sta_g2_txt = np.loadtxt(file_txt + '_STATIC_FLX_G2.txt')
    amp_g1_txt = np.loadtxt(file_txt + '_NOISE_AMP_G1.txt')
    amp_g2_txt = np.loadtxt(file_txt + '_NOISE_AMP_G2.txt')
    pha_g1_txt = np.loadtxt(file_txt + '_NOISE_PHASE_G1.txt')
    pha_g2_txt = np.loadtxt(file_txt + '_NOISE_PHASE_G2.txt')
    
    # Relative noise
    relnos_g1_txt = amp_g1_txt / sta_g1_txt * 100 
    relnos_g2_txt = amp_g2_txt / sta_g2_txt * 100 


    sta_line_g1 = sta_g1_txt.diagonal()
    sta_line_g2 = sta_g2_txt.diagonal()
    nos_line_g1 = relnos_g1_txt.diagonal()
    nos_line_g2 = relnos_g2_txt.diagonal()
    phs_line_g1 = pha_g1_txt.diagonal()
    phs_line_g2 = pha_g2_txt.diagonal()
       
    
    static_g1_line_txt.append(sta_line_g1)
    static_g2_line_txt.append(sta_line_g2)
    noise_g1_line_txt.append(nos_line_g1)
    noise_g2_line_txt.append(nos_line_g2)
    phase_g1_line_txt.append(phs_line_g1)
    phase_g2_line_txt.append(phs_line_g2)
    
#%% ---------------------------------------------------------------------------

# Static g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line, static_g1_line[i]/static_g1_line[i][0], style[i], label=labels[i])
for i in range(n_files_txt):
    ax1.plot(x_line, static_g1_line_txt[i]/static_g1_line_txt[i][0], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Static Flux g1 (AU)")
fig1.savefig(folder + problem + "_static_g1.pdf", format='pdf')


# Static g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_line_ftd, noise_g1_line_ftd, '-', label='Time-Domain')
for i in range(n_files):
    ax1.plot(x_line, static_g2_line[i]/static_g1_line[i][0], style[i], label=labels[i])
for i in range(n_files_txt):
    ax1.plot(x_line, static_g2_line_txt[i]/static_g1_line_txt[i][0], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Static Flux g2 (AU)")
fig1.savefig(folder + problem + "_static_g2.pdf", format='pdf')

# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line, noise_g1_line[i], style[i], label=labels[i])
for i in range(n_files_txt):
    ax1.plot(x_line, noise_g1_line_txt[i], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Relative Noise Magnitude (\%)")
fig1.savefig(folder + problem + "_noise_line_amp_g1.pdf", format='pdf')

# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line, noise_g2_line[i], style[i], label=labels[i])
for i in range(n_files_txt):
    ax1.plot(x_line, noise_g2_line_txt[i], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Relative Noise Magnitude (\%)")
fig1.savefig(folder + problem + "_noise_line_amp_g2.pdf", format='pdf')


# Print phase_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line, phase_g1_line[i], style[i], label=labels[i])
for i in range(n_files_txt):
    ax1.plot(x_line, phase_g1_line_txt[i], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Phase (deg)")
ax1.set_ylim([176.9, 177.3])
fig1.savefig(folder + problem + "_noise_line_pha_g1.pdf", format='pdf')


# Print phase_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_line_ftd, noise_g2_line_ftd, '-', label='Time-Domain')
for i in range(n_files):
    ax1.plot(x_line, phase_g2_line[i], style[i], label=labels[i])
for i in range(n_files_txt):
    ax1.plot(x_line, phase_g2_line_txt[i], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Phase (deg)")
ax1.set_ylim([176.9, 177.7])
fig1.savefig(folder + problem + "_noise_line_pha_g2.pdf", format='pdf')




#%% ---------------------------------------------------------------------------
# CON NORMALIZACIÓN

# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line, noise_g1_line[i]/noise_g1_line[i][0], style[i], label=labels[i])
for i in range(n_files_txt):
    ax1.plot(x_line, noise_g1_line_txt[i]/noise_g1_line_txt[i][0], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Relative Noise Magnitude (AU)")
fig1.savefig(folder + problem + "_noise_line_ampnorm_g1.pdf", format='pdf')

# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line, noise_g2_line[i]/noise_g2_line[i][0], style[i], label=labels[i])
for i in range(n_files_txt):
    ax1.plot(x_line, noise_g2_line_txt[i]/noise_g2_line_txt[i][0], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Relative Noise Magnitude (AU)")
fig1.savefig(folder + problem + "_noise_line_ampnorm_g2.pdf", format='pdf')


