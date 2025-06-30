#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Antoni Vidal
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
#%% ---------------------------------------------------------------------------

problem = '2D_UOX_FA_ex3'
folder  = 'exercise_3/'
case    = '3' 
omega   = '2w0'

# problem = '2D_UOX_FA_ex6'
# folder  = 'exercise_6/'
# case    = '6' 
# omega   = 'w0'
 
#
file_1 = folder + 'FEMFFUSIONFD_SP1_' + omega
file_2 = folder + 'FEMFFUSIONFD_SP3_' + omega
file_3 = folder + 'FEMFFUSIONFD_SP5_' + omega
file_4 = folder + 'FEMFFUSIONFD_SP7_' + omega
file_5 = 'FEMFFUSION_TD/FEMFFUSION_CASE' + case
#
##files = [file_1, file_2, file_3]
##labels = ['Diffusion', 'SP3', 'FSP3']
##style = ['*-', '.--', '+']
#
# files = [file_1, file_5, file_2, file_3, file_4 ]
# labels = ['SP1', 'SP1-TD', 'SP3', 'SP5', 'SP7']
# style = ['-', '.', '.-', '-^', '--+' ]
files = [file_1, file_2, file_3, file_4]
labels = ['SP1', 'SP3', 'SP5', 'SP7']
style = ['-',  '.-', '-^', '--+']

n_files = len(files)
assert(len(files)==len(labels) and len(labels)==n_files)
problem
if case == '2':
    nx = 138
    ny = 138
    x = [0.08,  # water_strip
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #1
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #2
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #3
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #4
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #5
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #6
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #7
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
    0.08]  # water_strip
    y = [0.08,  # water_strip
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #1
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #2
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #3
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #4
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #5
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #6
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #7
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
    0.08]  # water_strip
    
    line_idx = 53;

    assert(len(x)== nx)
    assert(len(y)== ny)
    
if case == '3':
    nx = 174
    ny = 138
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

    y = [0.08,  # water_strip
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #1
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #2
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #3
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #4
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #5
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #6
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #7
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
    0.08]  # water_strip
    assert(len(x)== nx)
    assert(len(y)== ny)
    line_idx = 53;

if case == '4':
    nx = 210
    ny = 138
if case == '5':
    nx = 210
    ny = 138
if case == '6':
    nx = 246
    ny = 138
 
    x = [0.08,
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #1
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #2
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #3
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #4
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #5
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #6 
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #7
    0.06430,  
    0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 
    0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
    0.16570, 0.16570,                                                #  Pin 8
    0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 
    0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
    0.06430,                                                               
    0.06430,
    0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
    0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
    0.16570, 0.16570,                                                #  Pin 9
    0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 
    0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
    0.06430,        
    0.06430,
    0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 
    0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 
    0.16570, 0.16570,                                               #  Pin 10
    0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,   
    0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
    0.06430,  
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #11
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #12
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #13
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #14
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #15
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #16
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #17
    0.08]

    y = [0.08,  # water_strip
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #1
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #2
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #3
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #4
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #5
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #6
    0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215, #7
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
    0.08]  # water_strip
    line_idx = 1 + 8*8+4;
    
n_cells = nx * ny
x = np.cumsum(x)
x = np.array(x)
y = np.cumsum(y)
y = np.array(y)
X, Y = np.meshgrid(x, y)

#%% ---------------------------------------------------------------------------
# # TIME DOMAINomega
# folder_txt   = '../2D_UOX_FA/time_domain/'
# file_txt_1 = folder_txt + '2D_UOX_diffusion_FE1'  
# file_txt_2 = folder_txt + '2D_UOX_sp3_sint'  
# #files_txt = [file_txt_1, file_txt_2]
# #style_td = ['+C0', '*C1'] 
# #labels_td = ['TD-Diff', 'TD-SP3']

# files_txt = []
# style_td = [] 
# labels_td = []
# n_files_txt = len(files_txt)


#        len(labels_td)==n_files_txt and 
#        len(style_td)==n_files_txt)

#%% ---------------------------------------------------------------------------
# GET FD DATA FROM TXT
sta_g1_line = []
sta_g2_line = []
amp_g1_line = []
amp_g2_line = []
pha_g1_line = []
pha_g2_line = []
rel_g1_line = []
rel_g2_line = []
rea_g1_line = []
rea_g2_line = []
ima_g1_line = []
ima_g2_line = []


sta_g1 = []
sta_g2 = []
amp_g1 = []
amp_g2 = []
pha_g1 = []
pha_g2 = []
rel_g1 = []
rel_g2 = []
rea_g1 = []
rea_g2 = []
ima_g1 = []
ima_g2 = []

for i in range(n_files):
    file = files[i]
    print(files[i], '...')
    
    sta_g1.append(np.loadtxt(file + '_STATIC_FLX_G1.txt'))
    sta_g2.append(np.loadtxt(file + '_STATIC_FLX_G2.txt'))    
    amp_g1.append(np.loadtxt(file + '_NOISE_AMP_G1.txt'))
    amp_g2.append(np.loadtxt(file + '_NOISE_AMP_G2.txt'))
    pha_g1.append(np.loadtxt(file + '_NOISE_PHASE_G1.txt'))
    pha_g2.append(np.loadtxt(file + '_NOISE_PHASE_G2.txt'))
    rea_g1.append(np.loadtxt(file + '_NOISE_REAL_G1.txt'))
    rea_g2.append(np.loadtxt(file + '_NOISE_REAL_G2.txt'))
    ima_g1.append(np.loadtxt(file + '_NOISE_IMAG_G1.txt'))
    ima_g2.append(np.loadtxt(file + '_NOISE_IMAG_G2.txt'))
    
    # Relative noise
    rel_g1.append(amp_g1[i] / sta_g1[i] * 100)
    rel_g2.append(amp_g2[i] / sta_g2[i] * 100)

    sta_g1_line.append(sta_g1[i][line_idx])
    sta_g2_line.append(sta_g2[i][line_idx])
    amp_g1_line.append(amp_g1[i][line_idx])
    amp_g2_line.append(amp_g2[i][line_idx])
    pha_g1_line.append(pha_g1[i][line_idx])
    pha_g2_line.append(pha_g2[i][line_idx])
    rel_g1_line.append(rel_g1[i][line_idx])
    rel_g2_line.append(rel_g2[i][line_idx])
    rea_g1_line.append(rea_g1[i][line_idx])
    rea_g2_line.append(rea_g2[i][line_idx])
    ima_g1_line.append(ima_g1[i][line_idx])
    ima_g2_line.append(ima_g2[i][line_idx])

#%% ---------------------------------------------------------------------------

# Static g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x, sta_g1_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x [cm]")
ax1.set_ylabel("Static Flux g1 (AU)")
fig1.savefig(folder + problem + '_' + omega +"_static_g1.pdf", format='pdf')


# Static g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_ftd, noise_g1_line_ftd, '-', label='Time-Domain')
for i in range(n_files):
    ax1.plot(x, sta_g2_line[i], style[i], label=labels[i])
# for i in range(n_files_txt):
#     ax1.plot(x, static_g2_line_txt[i]/static_g1_line_txt[i][0], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x [cm]")
ax1.set_ylabel("Static Flux g2 (AU)")
fig1.savefig(folder + problem + '_' + omega +"_static_g2.pdf", format='pdf')


# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x, rel_g1_line[i], style[i], label=labels[i])
# for i in range(n_files_txt):
#     ax1.plot(x, noise_g1_line_txt[i], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x [cm]")
ax1.set_ylabel(r"Relative Noise Magnitude [\%]")
fig1.savefig(folder + problem + '_' + omega +"_noise_line_amp_g1.pdf", format='pdf')


# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x, rel_g2_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x [cm]")
ax1.set_ylabel(r"Relative Noise Magnitude [\%]")
fig1.savefig(folder + problem +  '_' + omega  + "_noise_line_amp_g2.pdf", format='pdf')

# Print phase_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x, pha_g1_line[i], style[i], label=labels[i])
# for i in range(n_files_txt):
#     ax1.plot(x, phase_g1_line_txt[i], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x [cm]")
ax1.set_ylabel("Phase [deg]")
# ax1.set_ylim([176.9, 177.3])
fig1.savefig(folder + problem + '_' + omega  +"_noise_line_pha_g1.pdf", format='pdf')

# Print phase_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x_ftd, noise_g2_line_ftd, '-', label='Time-Domain')
for i in range(n_files):
    ax1.plot(x, pha_g2_line[i], style[i], label=labels[i])
# for i in range(n_files_txt):
#     ax1.plot(x, phase_g2_line_txt[i], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x [AU]")
ax1.set_ylabel("Phase [deg]")
# ax1.set_ylim([176.9, 177.7])
fig1.savefig(folder + problem + '_' + omega  +"_noise_line_pha_g2.pdf", format='pdf')


# Print real_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x, rea_g1_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x [cm]")
ax1.set_ylabel("Real Noise Field  [AU]")
# ax1.set_ylim([176.9, 177.3])
fig1.savefig(folder + problem + '_' + omega  +"_noise_line_rea_g1.pdf", format='pdf')

# Print real_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x, rea_g2_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x [cm]")
ax1.set_ylabel("Real Noise Field  [AU]")
# ax1.set_ylim([176.9, 177.7])
fig1.savefig(folder + problem + '_' + omega  +"_noise_line_rea_g2.pdf", format='pdf')


# Print imag_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x, ima_g1_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x [cm]")
ax1.set_ylabel("Imaginary Noise Field  [AU]")
# ax1.set_ylim([176.9, 177.3])
fig1.savefig(folder + problem + '_' + omega  +"_noise_line_ima_g1.pdf", format='pdf')

# Print imag_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x, ima_g2_line[i], style[i], label=labels[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x [cm]")
ax1.set_ylabel("Imaginary Noise Field  [AU]")
# ax1.set_ylim([176.9, 177.7])
fig1.savefig(folder + problem + '_' + omega  +"_noise_line_ima_g2.pdf", format='pdf')

#%% ---------------------------------------------------------------------------
# CON NORMALIZACIÓN

# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x, rel_g1_line[i]/rel_g1_line[i][0], style[i], label=labels[i])
# for i in range(n_files_txt):
#     ax1.plot(x, noise_g1_line_txt[i]/noise_g1_line_txt[i][0], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x [cm]")
ax1.set_ylabel("Normalized Relative Noise Magnitude (AU)")
fig1.savefig(folder + problem + '_' + omega  +"_noise_line_ampnorm_g1.pdf", format='pdf')


# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files):
    ax1.plot(x, rel_g2_line[i]/rel_g2_line[i][0], style[i], label=labels[i])
# for i in range(n_files_txt):
#     ax1.plot(x, noise_g2_line_txt[i]/noise_g2_line_txt[i][0], style_td[i], label=labels_td[i])
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x [cm]")
ax1.set_ylabel("Normalized Relative Noise Magnitude (AU)")
fig1.savefig(folder + problem + '_' + omega  +"_noise_line_ampnorm_g2.pdf", format='pdf')

