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
from scipy.interpolate import interp1d


plt.close('all')

params = {'backend': 'pdf',
#          'font.family': 'serif',
          'font.size': 14,
          'axes.labelsize': 12,
          'legend.fontsize': 12,
          'xtick.labelsize': 13,
          'ytick.labelsize': 13,
          'text.usetex': False,
          'lines.linewidth': 2,
          'lines.markersize': 10,
          'lines.markeredgewidth': 1,
          'legend.numpoints': 1, 
          'axes.formatter.useoffset': False,
          'figure.autolayout': True,
          }

rcParams.update(params)

# FREQUENCY DOMAIN

problem = '2D_C5G7_noise'
folder   = 'NN1.2/'
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


file_1     = folder +             'c5g7_ref1_FE2_DSP1.out' 
file_2     = folder +             'c5g7_ref1_FE2_DSP3.out' 
file_3     = folder +             'c5g7_ref1_FE2_DSP5.out' 
file_4     = folder +             'c5g7_ref1_FE2_DSP7.out' 
#
##
files = [file_1, file_2, file_3, file_4]
labels = ['SP1', 'SP3', 'SP5', 'SP7']
style = ['-', '--', '-.', ':']

n_files = len(files)
assert(len(files)==len(labels) and len(labels)==n_files)

# ====================================================================== #
## TIME DOMAIN
#
# 2D C5G7
ny = 51
nx = 51
n_cells = nx * ny
##
##
#folder_txt   = 'NN1.2/time_domain/'
#file_txt_1 = folder_txt + '2D_C5G7_noise_diff_FE2'  
#file_txt_2 = folder_txt + '2D_C5G7_noise_dsp3_FE2'  
#
#files_txt = [file_txt_1, file_txt_2]
#style_td = ['+C0', '*C1'] 
#labels_td = ['TD-DSP1', 'TD-DSP3']
#
#
#n_files_txt = len(files_txt)
#
#
#assert(len(files_txt)==n_files_txt and
#       len(labels_td)==n_files_txt and 
#       len(style_td)==n_files_txt)


cmap="viridis" 
#cmap="inferno"
#cmap="RdBu_r"

def get_antidiagonal(mat):
    
    diag = np.zeros(34)
    for i in range(34):# Debido al bug de los components
        diag[i] = mat[i][33-i]
            
    return diag

def get_antidiagonal_comp(mat):
    
    diag = np.zeros(34, dtype=np.complex_)
    for i in range(34):# Debido al bug de los components
        diag[i] = mat[i][33-i]
            
    return diag
#%% ---------------------------------------------------------------------------
# Line Geometry
x =[]
for i in range(34):
    x.append((0.5+ i) *1.26)
    
x = np.array(x)

x_line = np.sqrt(2)*x
xs=np.linspace(0, (34*1.26)*np.sqrt(2),51*10)

#%% ---------------------------------------------------------------------------
# GET FD DATA

noise_g1_line = []
noise_g7_line = []

phase_g1_line = []
phase_g7_line = []

static_g1_line= []
static_g7_line= []

shape = [51, 51]


for i in range(n_files):
    print(files[i], '...')
    # Get From .OUT
    sta_g1 = parse_file(files[i], 'Group 1 flux', n_max_lines=ny)
    sta_g7 = parse_file(files[i], 'Group 7 flux', n_max_lines=ny)
    
    sta_g1 = np.reshape(np.array(sta_g1), shape)
    sta_g7 = np.reshape(np.array(sta_g7), shape)

    
    nos_g1 = parse_file_complex(files[i], 'Flux Noise Group 1', n_max_lines=ny)
    nos_g7 = parse_file_complex(files[i], 'Flux Noise Group 7', n_max_lines=ny)
    
    nos_g1 = np.reshape(np.array(nos_g1), shape)
    nos_g7 = np.reshape(np.array(nos_g7), shape)
    
    # Relative noise
    relnos_g1 = nos_g1 / sta_g1 * 100 
    relnos_g7 = nos_g7 / sta_g7 * 100 

    sta_line_g1 = get_antidiagonal(sta_g1)
    sta_line_g7 = get_antidiagonal(sta_g7)
    
    nos_line_g1 = get_antidiagonal_comp(relnos_g1)
    nos_line_g7 = get_antidiagonal_comp(relnos_g7)
        

    inter = interp1d(x_line, sta_line_g1, kind="cubic", fill_value="extrapolate")  
    sta_line_g1 = inter(xs)
    inter = interp1d(x_line, sta_line_g7, kind="cubic", fill_value="extrapolate")
    sta_line_g7 = inter(xs)
    
    static_g1_line.append(sta_line_g1)
    static_g7_line.append(sta_line_g7)
    
    noise_g1_line.append(np.abs(nos_line_g1))
    noise_g7_line.append(np.abs(nos_line_g7))
        
    phase_g1_line.append(np.angle(nos_line_g1, deg=True))
    phase_g7_line.append(np.angle(nos_line_g7, deg=True))
  
    

    
##%% ---------------------------------------------------------------------------
## GET TXT DATA  - TIME DOMAIN
#noise_g1_line_txt = []
#noise_g7_line_txt = []
#
#phase_g1_line_txt = []
#phase_g7_line_txt = []
#
#static_g1_line_txt = []
#static_g7_line_txt = []
#
#
#
#for i in range(n_files_txt):
#    file_txt = files_txt[i]
#    print(files_txt[i], '...')
#    
#    sta_g1_txt = np.loadtxt(file_txt + '_STATIC_FLX_G1.txt')
#    sta_g7_txt = np.loadtxt(file_txt + '_STATIC_FLX_G7.txt')
#        
#    amp_g1_txt = np.loadtxt(file_txt + '_NOISE_AMP_G1.txt')
#    amp_g7_txt = np.loadtxt(file_txt + '_NOISE_AMP_G7.txt')
#
#    pha_g1_txt = np.loadtxt(file_txt + '_NOISE_PHASE_G1.txt')
#    pha_g7_txt = np.loadtxt(file_txt + '_NOISE_PHASE_G7.txt')
#    
#    # Relative nois5
#    relnos_g1_txt = amp_g1_txt / sta_g1_txt * 100 
#    relnos_g7_txt = amp_g7_txt / sta_g7_txt * 100 
#
#
#    sta_line_g1 = get_antidiagonal(sta_g1_txt)
#    sta_line_g7 = get_antidiagonal(sta_g7_txt)
#    
#    nos_line_g1 = get_antidiagonal(relnos_g1_txt)
#    nos_line_g7 = get_antidiagonal(relnos_g7_txt)
#
#    phs_line_g1 = get_antidiagonal(pha_g1_txt)
#    phs_line_g7 = get_antidiagonal(pha_g7_txt)
#       
#    
#    
#    static_g1_line_txt.append(sta_line_g1)
#    static_g7_line_txt.append(sta_line_g7)
#    
#    noise_g1_line_txt.append(nos_line_g1)
#    noise_g7_line_txt.append(nos_line_g7)
#    
#    phase_g1_line_txt.append(phs_line_g1)
#    phase_g7_line_txt.append(phs_line_g7)
    
    
#%% ---------------------------------------------------------------------------
# Datos SN paolo
from data_2D_C5G7_SN import static_g1, static_g2, amp_g1, amp_g2, pha_g1, pha_g2
style_sn = '-oC2'
label_sn = 'S64 (Yi, 2021)'

l = int(len(static_g1)/2)
x_sta_g1 = np.zeros(l)
sta_g1 = np.zeros(l)
for i in range(l):
    x_sta_g1[i] = static_g1[2*i]
    sta_g1[i] = static_g1[2*i+1]
x_sta_g1 = x_sta_g1 * xs[-1] / x_sta_g1[-1];
    
l = int(len(static_g2)/2)
x_sta_g2 = np.zeros(l)
sta_g2 = np.zeros(l)
for i in range(l):
    x_sta_g2[i] = static_g2[2*i]
    sta_g2[i] = static_g2[2*i+1]
x_sta_g2 = x_sta_g2 * xs[-1] / x_sta_g2[-1];

l = int(len(amp_g1)/2)
x_amp_g1 = np.zeros(l)
amp_g1_sn = np.zeros(l)
for i in range(l):
    x_amp_g1[i] = amp_g1[2*i]
    amp_g1_sn[i] = amp_g1[2*i+1]
x_amp_g1 = x_amp_g1 * xs[-1] / x_amp_g1[-1];


amp_g1_sn /= sta_g1
#inter = interp1d(x_sta_g1, sta_g1, kind="cubic", fill_value="extrapolate")
#sta_g1_sn = inter(xs)
#inter = interp1d(x_amp_g1, amp_g1_sn, kind="cubic", fill_value="extrapolate")
#amp_g1_sn = inter(xs)
#amp_g1_sn /= sta_g1_sn

    
l = int(len(amp_g2)/2)
x_amp_g2 = np.zeros(l)
amp_g2_sn = np.zeros(l)
for i in range(l):
    x_amp_g2[i] = amp_g2[2*i]
    amp_g2_sn[i] = amp_g2[2*i+1]
x_amp_g2 = x_amp_g2 * xs[-1] / x_amp_g2[-1];  

amp_g2_sn /= sta_g2

#inter = interp1d(x_sta_g2, sta_g2, kind="linear", fill_value="extrapolate")
#sta_g2_sn = inter(xs)
#inter = interp1d(x_amp_g2, amp_g2_sn, kind="linear", fill_value="extrapolate")
#amp_g2_sn_int = inter(xs)
#amp_g2_sn_int /= sta_g2_sn

    
l = int(len(pha_g1)/2)
x_pha_g1 = np.zeros(l)
pha_g1_sn = np.zeros(l)
for i in range(l):
    x_pha_g1[i] = pha_g1[2*i]
    pha_g1_sn[i] = pha_g1[2*i+1]
x_pha_g1 = x_pha_g1 * xs[-1] / x_pha_g1[-1];  

l = int(len(pha_g2)/2)
x_pha_g2 = np.zeros(l)
pha_g2_sn = np.zeros(l)
for i in range(l):
    x_pha_g2[i] = pha_g2[2*i]
    pha_g2_sn[i] = pha_g2[2*i+1]
x_pha_g2 = x_pha_g2 * xs[-1] / x_pha_g2[-1];  
    
    
#%% ---------------------------------------------------------------------------

# Static g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(xs, static_g1_line[i], style[i], label=labels[i])

ax1.plot(x_sta_g1, sta_g1 * 1.04*max(static_g1_line[0])/max(sta_g1), style_sn, label=label_sn, markersize=3,linewidth= 1)
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Static Flux g1 (AU)")
#ax1.set_xlim([0, 1.26*34 * np.sqrt(2)])
fig1.savefig(folder + problem + "_static_g1.pdf", format='pdf')


# Static g7
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(xs, static_g7_line[i], style[i], label=labels[i])
ax1.plot(x_sta_g2, sta_g2 * 1.04*max(static_g7_line[0])/max(sta_g2), style_sn, label=label_sn, markersize=3, linewidth= 1)
#for i in range(n_files_txt):
#    ax1.plot(xs, static_g7_line_txt[i], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Static Flux g7 (AU)")
ax1.set_xlim([0, 1.26*34 * np.sqrt(2)])
fig1.savefig(folder + problem + "_static_g7.pdf", format='pdf')

# Print phase_g7
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line, phase_g1_line[0][0] * phase_g1_line[i]/phase_g1_line[i][0], style[i], label=labels[i])
#for i in range(n_files_txt):
#    ax1.plot(x_line, phase_g1_line_txt[i], style_td[i], label=labels_td[i])
ax1.plot(x_pha_g1, -90+pha_g1_sn, style_sn, label=label_sn, markersize=3, linewidth= 1)
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Phase g1 (deg)")
#ax1.set_ylim([80, 90])
ax1.set_xlim([0, 1.26*34 * np.sqrt(2)])
fig1.savefig(folder + problem + "_noise_line_pha_g1.pdf", format='pdf')

# Print phase_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line, phase_g7_line[0][0] * phase_g7_line[i]/phase_g7_line[i][0], style[i], label=labels[i])
#for i in range(n_files_txt):
#    ax1.plot(x_line,phase_g7_line[1][0]  * phase_g7_line_txt[i]/phase_g7_line_txt[i][0], style_td[i], label=labels_td[i])
ax1.plot(x_pha_g2, -90+pha_g2_sn, style_sn, label=label_sn, markersize=3, linewidth= 1)
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
#for i in range(n_files_txt):
#    ax1.plot(x_line,   noise_g1_line_txt[i]/noise_g1_line_txt[i][0], style_td[i], label=labels_td[i])
ax1.plot(x_amp_g1, amp_g1_sn /amp_g1_sn[20], style_sn, label=label_sn, markersize=3, linewidth= 1)
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Relative Noise Magnitude g1 (AU)")
ax1.set_xlim([0, 1.26*34 * np.sqrt(2)])
fig1.savefig(folder + problem + "_noise_line_ampnorm_g1.pdf", format='pdf')

## Print noise_g7
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x_line,  noise_g7_line[i]/noise_g7_line[i][0], style[i], label=labels[i])
#for i in range(n_files_txt):
#    ax1.plot(x_line,   noise_g7_line_txt[i]/noise_g7_line_txt[i][0], style_td[i], label=labels_td[i])
ax1.plot(x_amp_g2, amp_g2_sn /amp_g2_sn[20], style_sn, label=label_sn, markersize=3, linewidth= 1)
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("Diagonal length (cm)")
ax1.set_ylabel("Relative Noise Magnitude g7 (AU)")
ax1.set_xlim([0, 1.26*34 * np.sqrt(2)])
fig1.savefig(folder + problem + "_noise_line_ampnorm_g7.pdf", format='pdf')







