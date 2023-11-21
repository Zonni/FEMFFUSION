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
from scipy.interpolate import interp1d

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
          'lines.markersize': 8,
          'lines.markeredgewidth': 1,
          'legend.numpoints': 1, 
          'axes.formatter.useoffset': False,
          'figure.autolayout': True,
          }

rcParams.update(params)

# FREQUENCY DOMAIN

problem = '2D_C5G7_noise'
folder   = '../2D_C5G7_noise/NN1.1/'
#
#file_1_sta = folder + 'static/' + 'c5g7_ref0_FE1_DSP1.out'
#file_1     = folder +             'c5g7_ref0_FE1_DSP1.out' 
#file_2_sta = folder + 'static/' + 'c5g7_ref0_FE1_DSP1.out'
#file_2     = folder + 'static/' + 'c5g7_ref0_FE1_DSP1.out' 
#file_3_sta = folder + 'static/' + 'c5g7_ref0_FE1_DSP3.out'
#file_3     = folder +             'c5g7_ref0_FE1_DSP3.out' 
#file_4_sta = folder + 'static/' + 'c5g7_ref0_FE1_DSP3.out'
#file_4     = folder + 'static/' + 'c5g7_ref0_FE1_DSP3.out' 
#
#files = [file_1, file_3]
#files_sta = [file_1_sta, file_3_sta]
#labels = ['FD-Diff', 'FD-DSP3']
#style = ['-', '--']


file_1_sta = folder + 'static/' + 'c5g7_ref1_FE2_DSP1.out'
file_1     = folder +             'c5g7_ref1_FE2_DSP1.out' 
file_2_sta = folder + 'static/' + 'c5g7_ref1_FE2_DSP3.out'
file_2     = folder +             'c5g7_ref1_FE2_DSP3.out' 

#
files = [file_1, file_2]
files_sta = [file_1_sta, file_2_sta]
labels = ['FD-DSP1', 'FD-DSP3']
style = ['-', '--']


n_files = len(files)
assert(len(files)==len(labels) and len(labels)==n_files)

# ====================================================================== #
# TIME DOMAIN

# 2D C5G7
ny = 51
nx = 51
n_cells = nx * ny
#
#
folder_txt   = '../2D_C5G7_noise/NN1.1/time_domain/'
file_txt_1 = folder_txt + '2D_C5G7_noise_diffusion_FE1'  
file_txt_2 = folder_txt + '2D_C5G7_noise_DSP3_ref1_FE2'  

files_txt = [file_txt_1, file_txt_2]
style_td = ['+C0', '*C1'] 
labels_td = ['TD-DSP1', 'TD-DSP3']



#folder_txt   = '../2D_C5G7_noise/NN1.1/time_domain/'
#file_txt_1 = folder_txt + '2D_C5G7_noise_diffusion_FE2'  
#file_txt_2 = folder_txt + '2D_C5G7_noise_DSP3_ref1_FE2'  

#files_txt = [file_txt_1, file_txt_2]
#style_td = ['+C0', '*C1'] 
#labels_td = ['TD-Diff', 'TD-SP3']

#files_txt = [file_txt_1]
#style_td = ['xC0'] 
#labels_td = ['TD-Diff']

#files_txt = []
#style_td = [] 
#labels_td = []


n_files_txt = len(files_txt)


assert(len(files_txt)==n_files_txt and
       len(labels_td)==n_files_txt and 
       len(style_td)==n_files_txt)


cmap="viridis" 
#cmap="inferno"
#cmap="RdBu_r"


#%% ---------------------------------------------------------------------------
# Line Geometry
x = [1.26,
     1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26,
     1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26,
     1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26,
     1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26,
     1.26, 1.26]
x = np.cumsum(np.array(x))
y = x.copy()
x_line = np.sqrt(2)*x


xs=np.linspace(0, 64.26*np.sqrt(2),51*10)

#%% ---------------------------------------------------------------------------
# GET FD DATA

noise_g1_line = []
#noise_g2_line = []
#noise_g3_line = []
#noise_g4_line = []
#noise_g5_line = []
#noise_g6_line = []
noise_g7_line = []

phase_g1_line = []
#phase_g2_line = []
#phase_g3_line = []
#phase_g4_line = []
#phase_g5_line = []
#phase_g6_line = []
phase_g7_line = []

static_g1_line= []
#static_g2_line= []
#static_g3_line= []
#static_g4_line= []
#static_g5_line= []
#static_g6_line= []
static_g7_line= []


for i in range(n_files):
    print(files[i], '...')
    # Get From .OUT
    sta_g1 = parse_file(files_sta[i], 'Group 1 flux', n_max_lines=ny)
#    sta_g2 = parse_file(files_sta[i], 'Group 2 flux', n_max_lines=ny)
#    sta_g3 = parse_file(files_sta[i], 'Group 3 flux', n_max_lines=ny)
#    sta_g4 = parse_file(files_sta[i], 'Group 4 flux', n_max_lines=ny)
#    sta_g5 = parse_file(files_sta[i], 'Group 5 flux', n_max_lines=ny)
#    sta_g6 = parse_file(files_sta[i], 'Group 6 flux', n_max_lines=ny)
    sta_g7 = parse_file(files_sta[i], 'Group 7 flux', n_max_lines=ny)
    
    sta_g1 = np.array(sta_g1)
#    sta_g2 = np.array(sta_g2)
#    sta_g3 = np.array(sta_g3)
#    sta_g4 = np.array(sta_g4)
#    sta_g5 = np.array(sta_g5)
#    sta_g6 = np.array(sta_g6)
    sta_g7 = np.array(sta_g7)

    
    nos_g1 = parse_file_complex(files[i], 'Flux Noise Group 1', n_max_lines=ny)
#    nos_g2 = parse_file_complex(files[i], 'Flux Noise Group 2', n_max_lines=ny)
#    nos_g3 = parse_file_complex(files[i], 'Flux Noise Group 3', n_max_lines=ny)        
#    nos_g4 = parse_file_complex(files[i], 'Flux Noise Group 4', n_max_lines=ny)
#    nos_g5 = parse_file_complex(files[i], 'Flux Noise Group 5', n_max_lines=ny)    
#    nos_g6 = parse_file_complex(files[i], 'Flux Noise Group 6', n_max_lines=ny)
    nos_g7 = parse_file_complex(files[i], 'Flux Noise Group 7', n_max_lines=ny)
    
    # Relative noise
    relnos_g1 = nos_g1 / sta_g1 * 100 
#    relnos_g2 = nos_g2 / sta_g2 * 100 
#    relnos_g3 = nos_g3 / sta_g3 * 100 
#    relnos_g4 = nos_g4 / sta_g4 * 100     
#    relnos_g5 = nos_g5 / sta_g5 * 100 
#    relnos_g6 = nos_g6 / sta_g6 * 100 
    relnos_g7 = nos_g7 / sta_g7 * 100 
    
    # Get diagonal
    sta_line_g1 = np.zeros(ny)
#    sta_line_g2 = np.zeros(ny)
#    sta_line_g3 = np.zeros(ny)
#    sta_line_g4 = np.zeros(ny)
#    sta_line_g5 = np.zeros(ny)
#    sta_line_g6 = np.zeros(ny)
    sta_line_g7 = np.zeros(ny)
    
    nos_line_g1 = np.zeros(ny, dtype=np.complex_)
#    nos_line_g2 = np.zeros(ny, dtype=np.complex_)
#    nos_line_g3 = np.zeros(ny, dtype=np.complex_)
#    nos_line_g4 = np.zeros(ny, dtype=np.complex_)
#    nos_line_g5 = np.zeros(ny, dtype=np.complex_)
#    nos_line_g6 = np.zeros(ny, dtype=np.complex_)
    nos_line_g7 = np.zeros(ny, dtype=np.complex_)
        
    for i in range(ny):
        sta_line_g1[i] = sta_g1[i*(nx +1)]
#        sta_line_g2[i] = sta_g2[i*(nx +1)]
#        sta_line_g3[i] = sta_g3[i*(nx +1)]
#        sta_line_g4[i] = sta_g4[i*(nx +1)]
#        sta_line_g5[i] = sta_g5[i*(nx +1)]
#        sta_line_g6[i] = sta_g6[i*(nx +1)]
        sta_line_g7[i] = sta_g7[i*(nx +1)]
        
        nos_line_g1[i] = relnos_g1[i*(nx +1)]
#        nos_line_g2[i] = relnos_g2[i*(nx +1)]
#        nos_line_g3[i] = relnos_g3[i*(nx +1)]
#        nos_line_g4[i] = relnos_g4[i*(nx +1)]
#        nos_line_g5[i] = relnos_g5[i*(nx +1)]
#        nos_line_g6[i] = relnos_g6[i*(nx +1)]
        nos_line_g7[i] = relnos_g7[i*(nx +1)]
        
        
    inter = interp1d(x_line, sta_line_g1, kind="cubic", fill_value="extrapolate")  
    sta_line_g1 = inter(xs)
    inter = interp1d(x_line, sta_line_g7, kind="cubic", fill_value="extrapolate")
    sta_line_g7 = inter(xs)
    
    static_g1_line.append(sta_line_g1)
#    static_g2_line.append(sta_line_g2)
#    static_g3_line.append(sta_line_g3)
#    static_g4_line.append(sta_line_g4)
#    static_g5_line.append(sta_line_g5)
#    static_g6_line.append(sta_line_g6)
    static_g7_line.append(sta_line_g7)
    
    noise_g1_line.append(np.abs(nos_line_g1))
#    noise_g2_line.append(np.abs(nos_line_g2))
#    noise_g3_line.append(np.abs(nos_line_g3))
#    noise_g4_line.append(np.abs(nos_line_g4))
#    noise_g5_line.append(np.abs(nos_line_g5))
#    noise_g6_line.append(np.abs(nos_line_g6))
    noise_g7_line.append(np.abs(nos_line_g7))
        
    phase_g1_line.append(np.angle(nos_line_g1, deg=True))
#    phase_g2_line.append(np.angle(nos_line_g2, deg=True))
#    phase_g3_line.append(np.angle(nos_line_g3, deg=True))
#    phase_g4_line.append(np.angle(nos_line_g4, deg=True))
#    phase_g5_line.append(np.angle(nos_line_g5, deg=True))
#    phase_g6_line.append(np.angle(nos_line_g6, deg=True))
    phase_g7_line.append(np.angle(nos_line_g7, deg=True))
  
#%% ---------------------------------------------------------------------------
# GET TXT DATA 
noise_g1_line_txt = []
#noise_g2_line_txt = []
#noise_g3_line_txt = []
#noise_g4_line_txt = []
#noise_g5_line_txt = []
#noise_g6_line_txt = []
noise_g7_line_txt = []

phase_g1_line_txt = []
#phase_g2_line_txt = []
#phase_g3_line_txt = []
#phase_g4_line_txt = []
#phase_g5_line_txt = []
#phase_g6_line_txt = []
phase_g7_line_txt = []

static_g1_line_txt = []
#static_g2_line_txt = []
#static_g3_line_txt = []
#static_g4_line_txt = []
#static_g5_line_txt = []
#static_g6_line_txt = []
static_g7_line_txt = []



for i in range(n_files_txt):
    file_txt = files_txt[i]
    print(files_txt[i], '...')
    
    sta_g1_txt = np.loadtxt(file_txt + '_STATIC_FLX_G1.txt')
#    sta_g2_txt = np.loadtxt(file_txt + '_STATIC_FLX_G2.txt')
#    sta_g3_txt = np.loadtxt(file_txt + '_STATIC_FLX_G3.txt')
#    sta_g4_txt = np.loadtxt(file_txt + '_STATIC_FLX_G4.txt')
#    sta_g5_txt = np.loadtxt(file_txt + '_STATIC_FLX_G5.txt')
#    sta_g6_txt = np.loadtxt(file_txt + '_STATIC_FLX_G6.txt')
    sta_g7_txt = np.loadtxt(file_txt + '_STATIC_FLX_G7.txt')
        
    amp_g1_txt = np.loadtxt(file_txt + '_NOISE_AMP_G1.txt')
#    amp_g2_txt = np.loadtxt(file_txt + '_NOISE_AMP_G2.txt')
#    amp_g3_txt = np.loadtxt(file_txt + '_NOISE_AMP_G3.txt')
#    amp_g4_txt = np.loadtxt(file_txt + '_NOISE_AMP_G4.txt')
#    amp_g5_txt = np.loadtxt(file_txt + '_NOISE_AMP_G5.txt')
#    amp_g6_txt = np.loadtxt(file_txt + '_NOISE_AMP_G6.txt')
    amp_g7_txt = np.loadtxt(file_txt + '_NOISE_AMP_G7.txt')

    pha_g1_txt = np.loadtxt(file_txt + '_NOISE_PHASE_G1.txt')
#    pha_g2_txt = np.loadtxt(file_txt + '_NOISE_PHASE_G2.txt')
#    pha_g3_txt = np.loadtxt(file_txt + '_NOISE_PHASE_G3.txt')
#    pha_g4_txt = np.loadtxt(file_txt + '_NOISE_PHASE_G4.txt')
#    pha_g5_txt = np.loadtxt(file_txt + '_NOISE_PHASE_G5.txt')
#    pha_g6_txt = np.loadtxt(file_txt + '_NOISE_PHASE_G6.txt')
    pha_g7_txt = np.loadtxt(file_txt + '_NOISE_PHASE_G7.txt')
    
    # Relative nois5
    relnos_g1_txt = amp_g1_txt / sta_g1_txt * 100 
#    relnos_g2_txt = amp_g2_txt / sta_g2_txt * 100 
#    relnos_g3_txt = amp_g3_txt / sta_g3_txt * 100 
#    relnos_g4_txt = amp_g4_txt / sta_g4_txt * 100 
#    relnos_g5_txt = amp_g5_txt / sta_g5_txt * 100 
#    relnos_g6_txt = amp_g6_txt / sta_g6_txt * 100 
    relnos_g7_txt = amp_g7_txt / sta_g7_txt * 100 


    sta_line_g1 = sta_g1_txt.diagonal()
#    sta_line_g2 = sta_g2_txt.diagonal()
#    sta_line_g3 = sta_g3_txt.diagonal()
#    sta_line_g4 = sta_g4_txt.diagonal()
#    sta_line_g5 = sta_g5_txt.diagonal()
#    sta_line_g6 = sta_g6_txt.diagonal()
    sta_line_g7 = sta_g7_txt.diagonal()
    
    nos_line_g1 = relnos_g1_txt.diagonal()
#    nos_line_g2 = relnos_g2_txt.diagonal()
#    nos_line_g3 = relnos_g3_txt.diagonal()
#    nos_line_g4 = relnos_g4_txt.diagonal()
#    nos_line_g5 = relnos_g5_txt.diagonal()
#    nos_line_g6 = relnos_g6_txt.diagonal()
    nos_line_g7 = relnos_g7_txt.diagonal()

    phs_line_g1 = pha_g1_txt.diagonal()
#    phs_line_g2 = pha_g2_txt.diagonal()
#    phs_line_g3 = pha_g3_txt.diagonal()
#    phs_line_g4 = pha_g4_txt.diagonal()
#    phs_line_g5 = pha_g5_txt.diagonal()
#    phs_line_g6 = pha_g6_txt.diagonal()
    phs_line_g7 = pha_g7_txt.diagonal()
       
    
#    inter = interp1d(x_line, sta_line_g7,
#                                         kind="cubic",
#                                         fill_value="extrapolate")
#    sta_line_g7 = inter(xs)
#    
    
    static_g1_line_txt.append(sta_line_g1)
#    static_g2_line_txt.append(sta_line_g2)
#    static_g3_line_txt.append(sta_line_g3)
#    static_g4_line_txt.append(sta_line_g4)
#    static_g5_line_txt.append(sta_line_g5)
#    static_g6_line_txt.append(sta_line_g6)
    static_g7_line_txt.append(sta_line_g7)
    
    noise_g1_line_txt.append(nos_line_g1)
#    noise_g2_line_txt.append(nos_line_g2)
#    noise_g3_line_txt.append(nos_line_g3)
#    noise_g4_line_txt.append(nos_line_g4)
#    noise_g5_line_txt.append(nos_line_g5)
#    noise_g6_line_txt.append(nos_line_g6)
    noise_g7_line_txt.append(nos_line_g7)
    
    phase_g1_line_txt.append(phs_line_g1)
#    phase_g2_line_txt.append(phs_line_g2)
#    phase_g3_line_txt.append(phs_line_g3)
#    phase_g4_line_txt.append(phs_line_g4)
#    phase_g5_line_txt.append(phs_line_g5)
#    phase_g6_line_txt.append(phs_line_g6)
    phase_g7_line_txt.append(phs_line_g7)
    
    
#%% ---------------------------------------------------------------------------

# Static g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(xs, static_g1_line[i], style[i], label=labels[i])
#for i in range(n_files_txt):
#    ax1.plot(xs, static_g7_line_txt[i], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Static Flux g1 (AU)")
ax1.set_xlim([0, 1.26*34 * np.sqrt(2)])
fig1.savefig(folder + problem + "_static_g1.pdf", format='pdf')


# Static g7
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(xs, static_g7_line[i], style[i], label=labels[i])
#for i in range(n_files_txt):
#    ax1.plot(xs, static_g7_line_txt[i], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Static Flux g7 (AU)")
ax1.set_xlim([0, 1.26*34 * np.sqrt(2)])
fig1.savefig(folder + problem + "_static_g7.pdf", format='pdf')


## rint noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line, noise_g1_line[i], style[i], label=labels[i])
for i in range(n_files_txt):
    ax1.plot(x_line, noise_g1_line_txt[i], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Relative Noise Magnitude g1 (\%)")
ax1.set_xlim([0, 1.26*34 * np.sqrt(2)])
fig1.savefig(folder + problem + "_noise_line_amp_g1.pdf", format='pdf')

# Print noise_g7
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line, noise_g7_line[i], style[i], label=labels[i])
for i in range(n_files_txt):
    ax1.plot(x_line, noise_g7_line_txt[i], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Relative Noise Magnitude g7 (\%)")
ax1.set_xlim([0, 1.26*34 * np.sqrt(2)])
fig1.savefig(folder + problem + "_noise_line_amp_g7.pdf", format='pdf')



# Print phase_g7
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line, phase_g1_line[1][0] * phase_g1_line[i]/phase_g1_line[i][0], style[i], label=labels[i])
for i in range(n_files_txt):
    ax1.plot(x_line,phase_g1_line[1][0]  * phase_g1_line_txt[i]/phase_g1_line_txt[i][0], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Phase g1 (deg)")
ax1.set_ylim([80, 90])
ax1.set_xlim([0, 1.26*34 * np.sqrt(2)])
fig1.savefig(folder + problem + "_noise_line_pha_g1.pdf", format='pdf')

# Print phase_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line, phase_g7_line[1][0] * phase_g7_line[i]/phase_g7_line[i][0], style[i], label=labels[i])
for i in range(n_files_txt):
    ax1.plot(x_line,phase_g7_line[1][0]  * phase_g7_line_txt[i]/phase_g7_line_txt[i][0], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Phase g7 (deg)")
ax1.set_ylim([80, 90])
ax1.set_xlim([0, 1.26*34 * np.sqrt(2)])
fig1.savefig(folder + problem + "_noise_line_pha_g7.pdf", format='pdf')



#%% ---------------------------------------------------------------------------
# CON NORMALIZACIÓN

# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line,  noise_g1_line[i]/noise_g1_line[i][0], style[i], label=labels[i])
for i in range(n_files_txt):
    ax1.plot(x_line,   noise_g1_line_txt[i]/noise_g1_line_txt[i][0], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Relative Noise Magnitude g1 (AU)")
ax1.set_xlim([0, 1.26*34 * np.sqrt(2)])
fig1.savefig(folder + problem + "_noise_line_ampnorm_g1.pdf", format='pdf')

# Print noise_g7
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line,  noise_g7_line[i]/noise_g1_line[i][0], style[i], label=labels[i])
for i in range(n_files_txt):
    ax1.plot(x_line,   noise_g7_line_txt[i]/noise_g1_line_txt[i][0], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Relative Noise Magnitude g7 (AU)")
fig1.savefig(folder + problem + "_noise_line_ampnorm_g7.pdf", format='pdf')




