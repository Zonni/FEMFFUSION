#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.insert(1, '../postprocess/')
import matplotlib.pyplot as plt
#import matplotlib.colors
from matplotlib import rcParams
from utils import parse_file, parse_file_complex, parse_file_opt
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
problem = '3D_ACP'

if (problem == '3D_ACP'):
    # file_1_fd = '1D_UOX_FA_ex2_diff.out'
    # files_fd = [file_1_fd ]
    # labels_fd = ['FD']
    # style_fd = ['-']
    
    # Time-Domain
    # file_sta_1 = 'ref_td/1D_UOX_FA_ex2_diffout'  
    # file_nos_1 = 'ref_td/1D_UOX_FA_ex2_diff.outnos'  
    # file_sta_1 = 'ref_td/1D_UOX_FA_ex2_diff_sintout'  
    # file_nos_1 = 'ref_td/1D_UOX_FA_ex2_diff_sint.outnos' 
    file_sta_1 = 'ACP100_FE2_td_FOM.out'  
    file_nos_1 = 'ACP100_FE2_td_FOM.out.nos'
    
    file_sta_2 = 'ACP100_FE2_td_ROM5.out'  
    file_nos_2 = 'ACP100_FE2_td_ROM5.out.nos'
    file_sta_3 = 'ACP100_FE2_td_ROM10.out'  
    file_nos_3 = 'ACP100_FE2_td_ROM10.out.nos'
    file_sta_4 = 'ACP100_FE2_td_ROM50.out'  
    file_nos_4 = 'ACP100_FE2_td_ROM50.out.nos'
    file_sta_5 = 'ACP100_FE2_td_ROM100.out'  
    file_nos_5 = 'ACP100_FE2_td_ROM100.out.nos'
    
    # files_sta_td = [file_sta_3]
    # files_nos_td = [file_nos_3]
    # labels_td = ['TD-ROM-10']
    # style_td = ['-', '--', '.-', '--', '.-']
    
    
    # files_sta_td = [file_sta_1, file_sta_3]
    # files_nos_td = [file_nos_1, file_nos_3]
    # labels_td = ['TD-FOM', 'TD-ROM-10']
    # style_td = ['-', '--', '.-', '--', '.-']
    
    
    files_sta_td = [file_sta_1,  file_sta_3, file_sta_4, file_sta_5]
    files_nos_td = [file_nos_1,  file_nos_3, file_nos_4, file_nos_5]
    labels_td = ['TD-FOM', 'TD-ROM-10',  'TD-ROM-50', 'TD-ROM-100']
    style_td = ['-', '--', '.-', '--', '.-']

# n_files_fd = len(files_fd)
n_files_td = len(files_nos_td)

assert(len(files_sta_td) == len(files_nos_td))
assert(len(files_sta_td) == len(labels_td))
# assert(len(labels_fd) == n_files_fd)


#%% ---------------------------------------------------------------------------
# GET MESH0.9099

n_cells = 1358
ny = 11 
nz = 14
# nx = 138
# ny= 1
# n_cells = nx * ny
# x = [0.08,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.13215, 0.13215, 0.18285, 0.18285, 0.18285, 0.18285, 0.13215, 0.13215,
#      0.08]
# assert(len(x)== nx)

# x = np.cumsum(np.array(x))
# x = np.insert(x, 0, 0.0) # Insert a 0 in the first position
# x = (x[:-1] + x[1:]) / 2  #Compute the average of consecutive elements



# #%% ---------------------------------------------------------------------------
# # FREQUENCY DOMAIN

# noise_g1_line_fd = []
# noise_g2_line_fd = []
# phase_g1_line_fd = []
# phase_g2_line_fd = []
# static_g1_line_fd = []
# static_g2_line_fd = []

# for i in range(len(files_fd)):
#     print('FD', files_fd[i], '...')
#     # Get From .OUT
#     sta_g1_fd = parse_file(files_fd[i], 'Group 1 flux', n_max_lines=ny)
#     sta_g2_fd = parse_file(files_fd[i], 'Group 2 flux', n_max_lines=ny)
#     sta_g1_fd = np.array(sta_g1_fd)
#     sta_g2_fd = np.array(sta_g2_fd)
    
#     nos_g1_fd = parse_file_complex(files_fd[i], 'Flux Noise Group 1', n_max_lines=ny)
#     nos_g2_fd = parse_file_complex(files_fd[i], 'Flux Noise Group 2', n_max_lines=ny)
    
#     # Relative noise
#     relnos_g1_fd = nos_g1_fd / sta_g1_fd * 100 
#     relnos_g2_fd = nos_g2_fd / sta_g2_fd * 100 
    
#     static_g1_line_fd.append(sta_g1_fd)
#     static_g2_line_fd.append(sta_g2_fd)
    
#     noise_g1_line_fd.append(np.abs(relnos_g1_fd))
#     noise_g2_line_fd.append(np.abs(relnos_g2_fd))
#     phase_g1_line_fd.append(np.angle(relnos_g1_fd, deg=True))
#     phase_g2_line_fd.append(np.angle(relnos_g2_fd, deg=True))


def interpolate(t_ref, y_ref, t1, y1):
    """
    Interpolates a signal to the first grid.
    """
    f1 = interp1d(t1, y1, kind='cubic', fill_value="extrapolate") 
    y1_interp = f1(t_ref)

    return y1_interp

#%% ---------------------------------------------------------------------------
#  TIME DOMAIN

noise_g1_line_td = []
noise_g2_line_td = []
phase_g1_line_td = []
phase_g2_line_td = []
static_g1_line_td = []
static_g2_line_td = []
noise_g1  = []
noise_g2  = []
sim_time = []

noise_interp_g1 = []
noise_interp_g2 = []

for i in range(n_files_td):
    print('Time Domain', files_nos_td[i], '...')
    
    sta_g1 = parse_file(files_sta_td[i], 'Group 1 flux', n_max_lines=ny*nz)
    sta_g2 = parse_file(files_sta_td[i], 'Group 2 flux', n_max_lines=ny*nz)
    sta_g1 = np.array(sta_g1)
    sta_g2 = np.array(sta_g2)
        
    time_fem = parse_file(files_sta_td[i], begin='Time vector', n_max_lines=1)
    time_fem = time_fem[:-1]
    n_steps = len(time_fem)
    
    current_noise_g1 = np.zeros([n_steps, n_cells])
    current_noise_g2 = np.zeros([n_steps, n_cells])
    # Open the file 
    with open(files_nos_td[i], 'r') as f:    
        for st in range(n_steps):
           
            #  Print something each processed 1e4 steps
            if (st % 1e4 == 0):
                print('Step: ', st)
                
            # Pass 'f' directly. It will start searching from wherever it left off!
            current_noise_g1[st] = parse_file_opt(f, 
                                      'Noise of group 1 time step ' + str(st), 
                                      n_max_lines=ny*nz)
                    
            current_noise_g2[st] = parse_file_opt(f, 
                                      'Noise of group 2 time step ' + str(st), 
                                      n_max_lines=ny*nz)
    
    current_noise_g1 = np.transpose(current_noise_g1)
    current_noise_g2 = np.transpose(current_noise_g2) 
    
    sim_time.append(time_fem)
    noise_g1.append(current_noise_g1)
    noise_g2.append(current_noise_g2)
    
    # if (i > 0):
    noise_interp_g1.append(interpolate(sim_time[0], noise_g1[0], sim_time[i], noise_g1[i]))
    noise_interp_g2.append(interpolate(sim_time[0], noise_g2[0], sim_time[i], noise_g2[i]))

 
    # freq   = np.fft.rfftfreq(n_steps, d=time_fem[1])
    # fft_g1 = np.fft.rfft(noise_g1) * 2.0/ n_steps 
    # fft_g2 = np.fft.rfft(noise_g2) * 2.0/ n_steps 
    
    # # We cut at looking_freq Hz
    # cut_freq = int (looking_freq * n_steps * time_fem[1])
    # assert(freq[cut_freq] == looking_freq)
    # noise_g1 = np.zeros([n_cells], dtype='cfloat')
    # noise_g2 = np.zeros([n_cells], dtype='cfloat')
    # for node in range(n_cells):
    #     noise_g1[node] = fft_g1[node][cut_freq] 
    #     noise_g2[node] = fft_g2[node][cut_freq]
    
    # amp_g1 = np.abs(noise_g1) / sta_g1 * 100
    # amp_g2 = np.abs(noise_g2) / sta_g2 * 100
    # # FEMFUSSION PERTURBATION is a SINE not a cosine
    # # This way we convert it
    # pha_g1 = np.angle(noise_g1, deg=True)
    # pha_g2 = np.angle(noise_g2, deg=True)  
    
    # static_g1_line_td.append(sta_g1)
    # static_g2_line_td.append(sta_g2)
    # noise_g1_line_td.append(amp_g1)
    # noise_g2_line_td.append(amp_g2)
    # phase_g1_line_td.append(pha_g1)
    # phase_g2_line_td.append(pha_g2)

#%% ---------------------------------------------------------------------------
#  ERRORS
for i in range(n_files_td):
    # Formula: RMS = sqrt(mean((y1 - y2)^2))
    errors_g1 = (noise_interp_g1[i] - noise_interp_g1[0]) /np.max(noise_interp_g1[0])
    rms_error_g1 = np.sqrt(np.mean(errors_g1**2)) * 100
    
    errors_g2 = (noise_interp_g2[i] - noise_interp_g2[0]) / np.max(noise_interp_g2[0])
    rms_error_g2 = np.sqrt(np.mean(errors_g2**2)) * 100
    
    # Using f-strings to format to 2 decimal places
    print(f"{i}  RMS ERROR G1: {rms_error_g1:.2f}  RMS ERROR G2: {rms_error_g2:.2f}")
    
#%% ---------------------------------------------------------------------------
#  PLOTS

# Plane 7
# plane = 8 # 7th plane 
# n_cells_per_plane = 97

middle_cell_idx = int((8*97+97/2) / 2)   
# 
# middle_cell_idx = int(n_cells/2) # cell 679 
#  Noise 1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(n_files_td):   
    ax1.plot(sim_time[i], noise_g1[i][middle_cell_idx], style_td[i], label=labels_td[i], color=colors[i])
ax1.grid(True)
ax1.legend(loc='upper right')
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Noise $\delta\phi_1$")
fig1.savefig(problem + "_td_g1.pdf", format='pdf')


# Noise 2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)

for i in range(n_files_td):   
    ax1.plot(sim_time[i], noise_g2[i][middle_cell_idx], style_td[i], label=labels_td[i], color=colors[i])
ax1.grid(True)
ax1.legend(loc='upper right')
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Noise $\delta\phi_2$")
fig1.savefig(problem + "_td_g2.pdf", format='pdf')



# middle_cell_idx = int(n_cells/2) # cell 679 
#  Noise 1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(sim_time[0], noise_g1[0][middle_cell_idx]/1000, style_td[0], label=labels_td[0], color=colors[0])
ax1.grid(True)
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Time domain Noise, $\delta\phi_1$")
fig1.savefig(problem + "_td_g1.png", format='png')


