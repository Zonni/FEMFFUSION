#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 13:03:14 2020

@author: zonni
"""

from utils import parse_vtk_file, parse_vtk_grid
import matplotlib.pyplot as plt
#import matplotlib.colors
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
folder   = '../2D_test_vibration/'
file_bpr = folder + '2D_test.out.vtk'
file_rfd = folder + '2D_test_ref.out.vtk'
file_prc = folder + '2D_test.out.parcs'

file_out = folder + 'timedomain/2D_test_ref3.out'  
static_file = folder + 'timedomain/2D_test_ref3.out.vtk'  
looking_freq = 1.0
problem = '2D_test'

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
    if (y[p] == mid_y):
        x_line_bpr.append(x[p])
        noise_g1_line_bpr.append(noise_g1[p])
        noise_g2_line_bpr.append(noise_g2[p])
        phase_g1_line_bpr.append(phase_g1[p])
        phase_g2_line_bpr.append(phase_g2[p])

#%% ---------------------------------------------------------------------------

# n_levels = 10
## Print noise_g1
#fig1 = plt.figure()
#ax1 = fig1.add_subplot(1, 1, 1)
#ax1.tricontour(x, y, noise_g1, levels=n_levels, linewidths=0.5, colors='k')
#cntr2 = ax1.tricontourf(x, y, noise_g1, levels=n_levels, cmap=cmap)
#ax1.set_aspect('equal')
#fig1.colorbar(cntr2, ax=ax1)
#ax1.set_ylabel("y (cm)")
#ax1.set_xlabel("x (cm)")
#ax1.set_title("Relative Fast Noise Magnitude (\%)")
#fig1.savefig(problem + "_noisevtk_g1_amp.pdf", format='pdf')
#
#
## Print noise_g2
#fig2 = plt.figure()
#ax2 = fig2.add_subplot(1, 1, 1)
#ax2.tricontour(x, y, noise_g2, levels=n_levels, linewidths=0.5, colors='k')
#cntr2 = ax2.tricontourf(x, y, noise_g2, levels=n_levels, cmap=cmap)
#fig2.colorbar(cntr2, ax=ax2)
#ax2.set_aspect('equal')
#ax2.set_ylabel("y (cm)")
#ax2.set_xlabel("x (cm)")
#ax2.set_title("Relative Thermal Noise Magnitude (\%)")
#fig2.savefig(problem + "_noise_g2_amp.pdf", format='pdf')
#
#
## Print noise_phase_g1
#norm = matplotlib.colors.Normalize(vmin=0, vmax=360)
#fig1 = plt.figure()
#ax1 = fig1.add_subplot(1, 1, 1)
#ax1.tricontour(x, y, phase_g1, levels=n_levels, linewidths=0.5, colors='k')
#cntr2 = ax1.tricontourf(x, y, phase_g1, levels=n_levels, cmap=cmap)
#fig1.colorbar(cntr2, ax=ax1)
#ax1.set_aspect('equal')
#ax1.set_ylabel("y (cm)")
#ax1.set_xlabel("x (cm)")
#ax1.set_title("Fast Noise Phase (deg)")
#fig1.savefig(problem + "_noisevtk_g1_phase.pdf", format='pdf')
#
## Print noise_phase_g2
#fig2 = plt.figure()
#ax2 = fig2.add_subplot(1, 1, 1)
#ax2.tricontour(x, y, phase_g2, levels=n_levels, linewidths=0.5, colors='k')
#cntr2 = ax2.tricontourf(x, y, phase_g2, levels=n_levels, cmap=cmap)
#fig2.colorbar(cntr2, ax=ax2)
#ax2.set_aspect('equal')
#ax2.set_ylabel("y (cm)")
#ax2.set_xlabel("x (cm)")
#ax2.set_title("Thermal Noise Phase (deg)")
#fig2.savefig(problem + "_noisevtk_g2_phase.pdf", format='pdf')

#%% ---------------------------------------------------------------------------
# Frequency domain reference  data
# Get From VTK
[x, y, z] = parse_vtk_grid(file_rfd)
stati_g1 = parse_vtk_file(file_rfd, "Static_Flux_g1")
stati_g2 = parse_vtk_file(file_rfd, "Static_Flux_g2")
noise_g1 = parse_vtk_file(file_rfd, "Noise_g1_Magnitude")
noise_g2 = parse_vtk_file(file_rfd, "Noise_g2_Magnitude")
phase_g1 = parse_vtk_file(file_rfd, "Noise_g1_Phase")
phase_g2 = parse_vtk_file(file_rfd, "Noise_g2_Phase")

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
x_line_rfd = []
noise_g1_line_rfd = []
noise_g2_line_rfd = []
phase_g1_line_rfd = []
phase_g2_line_rfd = []

for p in range(len(y)):
    if (y[p] == mid_y):
        x_line_rfd.append(x[p])
        noise_g1_line_rfd.append(noise_g1[p])
        noise_g2_line_rfd.append(noise_g2[p])
        phase_g1_line_rfd.append(phase_g1[p])
        phase_g2_line_rfd.append(phase_g2[p])
        
    
#%% ---------------------------------------------------------------------------
from utils import parse_fluxes_parcs_2D, parse_parcs_geometry
import numpy as np
delta_t_parcs = (1.0 / looking_freq) *  0.01

fl_parcs_fs, fl_parcs_th = parse_fluxes_parcs_2D(file_prc)
x_parcs, y_parcs = parse_parcs_geometry(file_prc)

# Remove last time step
fl_parcs_fs = fl_parcs_fs[0:-1]
fl_parcs_th = fl_parcs_th[0:-1]

n_steps_parcs = len(fl_parcs_th)
n_nodes_y = len(fl_parcs_th[0])
n_nodes_x = len(fl_parcs_th[0][0])
#
time_parcs = delta_t_parcs *  np.array(range(n_steps_parcs))

# Get Normalized static flux
st_flux_parcs_fs =  fl_parcs_fs[0].copy()
st_flux_parcs_th =  fl_parcs_th[0].copy()
##norm = np.max(static_flux_parcs_fs)

# Compute time and space dependent neutron noise
for st in range(n_steps_parcs) :
    for ny in range(n_nodes_y):
        for nx in range(n_nodes_x):
            fl_parcs_fs[st][ny][nx] -= st_flux_parcs_fs[ny][nx]
            fl_parcs_th[st][ny][nx] -= st_flux_parcs_th[ny][nx]

# Also, set the noise in time step 1 olso to 0
for ny in range(n_nodes_y):
    for nx in range(n_nodes_x):
        noise_fs_st1 = fl_parcs_fs[1][ny][nx]
        noise_th_st1 = fl_parcs_th[1][ny][nx]
        for st in range(1, n_steps_parcs) :
            fl_parcs_fs[st][ny][nx] -= noise_fs_st1
            fl_parcs_th[st][ny][nx] -= noise_th_st1
            

# Get noise at detectors
noise_parcs_fs = []
noise_parcs_th = []

idx = np.where((y_parcs == 30.5) + (y_parcs == 29.5))[0]

for nx in range(n_nodes_x):
    noise_parcs_fs.append(np.zeros(n_steps_parcs))
    noise_parcs_th.append(np.zeros(n_steps_parcs))
    
    for st in range(n_steps_parcs) : # Relative nosie
        noise_parcs_fs[nx][st] = 100 * (0.5 * 
                    fl_parcs_fs[st][idx[0]][nx] / st_flux_parcs_fs[idx[0]][nx] +
                    0.5 * 
                    fl_parcs_fs[st][idx[1]][nx] / st_flux_parcs_fs[idx[1]][nx])
        
        noise_parcs_th[nx][st] = 100 * (0.5 * 
                    fl_parcs_th[st][idx[0]][nx] / st_flux_parcs_th[idx[0]][nx] 
                    + 0.5 * 
                    fl_parcs_th[st][idx[1]][nx] / st_flux_parcs_th[idx[1]][nx])


# Compute Fast Fourier Transform
freq_prc = np.fft.rfftfreq(n_steps_parcs, d=delta_t_parcs)
noise_fft_fs = np.fft.rfft(noise_parcs_fs) * 2.0 / n_steps_parcs
noise_fft_th = np.fft.rfft(noise_parcs_th) * 2.0 / n_steps_parcs

# We cut at looking_freq Hz
cut_freq = int (looking_freq * n_steps_parcs * delta_t_parcs)
assert(abs(freq_prc[cut_freq] - looking_freq) < 1e-4)

# Get noise amplitude and phase
noise_g1_line_prc = np.zeros([n_nodes_x])
noise_g2_line_prc = np.zeros([n_nodes_x])
phase_g1_line_prc = np.zeros([n_nodes_x])
phase_g2_line_prc = np.zeros([n_nodes_x])
for node in range(n_nodes_x):
    noise_g1_line_prc[node] = abs(noise_fft_fs[node][cut_freq])
    noise_g2_line_prc[node] = abs(noise_fft_th[node][cut_freq])
    phase_g1_line_prc[node] = np.angle(noise_fft_fs[node][cut_freq], deg=True)
    phase_g2_line_prc[node] = np.angle(noise_fft_th[node][cut_freq], deg=True)
        
x_line_prc = x_parcs
#%% ---------------------------------------------------------------------------
# FEMFFUSION TIME DOMAIN
from utils import parse_file
import numpy as np


steps_fem = range(0, 300)
files =[]
for i, step in enumerate(steps_fem):
    files.append(file_out+ str(step) +  '.vtk')
    
time_fem = parse_file(file_out, begin='Time vector', n_max_lines=1)
if (len(time_fem) == (len(steps_fem) + 1)):
    time_fem = time_fem[:-1]
n_nodes_fem = len(parse_vtk_grid(files[0])[0])
n_steps_fem = len(steps_fem)
delta_t_fem = time_fem[1]

noise_fs = np.zeros([n_nodes_fem, n_steps_fem])
noise_th = np.zeros([n_nodes_fem, n_steps_fem])

for t, f in enumerate(files):
    ffs = parse_vtk_file(f, 'noise_g1')
    fth = parse_vtk_file(f, 'noise_g2')
    for node in range(n_nodes_fem):   
        noise_fs[node][t] = ffs[node] 
        noise_th[node][t] = fth[node] 

for t, f in enumerate(files):
    if (t==0.0):
        continue
    for node in range(n_nodes_fem):   
        noise_fs[node][t] -=  noise_fs[node][1]
        noise_th[node][t] -=  noise_th[node][1]

# Compute Fast Fourier Transform
freq_fem = np.fft.rfftfreq(n_steps_fem, d=delta_t_fem)
noise_fft_fs = np.fft.rfft(noise_fs) * 2.0 / n_steps_fem
noise_fft_th = np.fft.rfft(noise_th) * 2.0 / n_steps_fem

# We cut at looking_freq Hz
cut_freq = int (looking_freq * n_steps_fem * delta_t_fem)
assert(abs(freq_fem[cut_freq] - looking_freq) < 1e-4)

# Get noise amplitude and phase
noise_g1 = np.zeros([n_nodes_fem])
noise_g2 = np.zeros([n_nodes_fem])
phase_g1 = np.zeros([n_nodes_fem])
phase_g2 = np.zeros([n_nodes_fem])
for node in range(n_nodes_fem):
    noise_g1[node] = abs(noise_fft_fs[node][cut_freq])
    noise_g2[node] = abs(noise_fft_th[node][cut_freq])
    phase_g1[node] = np.angle(noise_fft_fs[node][cut_freq], deg=True)
    phase_g2[node] = np.angle(noise_fft_th[node][cut_freq], deg=True)
        
# Get Grid 
[x, y, z] = parse_vtk_grid(static_file)
# Static Fluxes
stati_g1 = parse_vtk_file(static_file, 'phi_g1_eig_1')
stati_g2 = parse_vtk_file(static_file, 'phi_g2_eig_1')

# Remove repeated data
stati_g1 = remove_repeated_data_point(x,y,z, stati_g1)
stati_g2 = remove_repeated_data_point(x,y,z, stati_g2)
noise_g1 = remove_repeated_data_point(x,y,z, noise_g1)
noise_g2 = remove_repeated_data_point(x,y,z, noise_g2)
phase_g1 = remove_repeated_data_point(x,y,z, phase_g1)
phase_g2 = remove_repeated_data_point(x,y,z, phase_g2)
x, y, z  = remove_repeated_point(x,y,z)

# Relative noise
noise_g1 = noise_g1 / stati_g1 * 100
noise_g2 = noise_g2 / stati_g2 * 100

# Midline
mid_y = max(y) / 2
x_line_ftd = []
noise_g1_line_ftd = []
noise_g2_line_ftd = []
phase_g1_line_ftd = []
phase_g2_line_ftd = []

for p in range(len(y)):
    if (y[p] == mid_y):
        x_line_ftd.append(x[p])
        noise_g1_line_ftd.append(noise_g1[p])
        noise_g2_line_ftd.append(noise_g2[p])
        phase_g1_line_ftd.append(phase_g1[p])
        phase_g2_line_ftd.append(phase_g2[p])
        

#%% ---------------------------------------------------------------------------

# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(x_line_ftd, noise_g1_line_ftd, '-', label='Time-Domain')
ax1.plot(x_line_rfd, noise_g1_line_rfd, 'X', label='Cell wise FD')
ax1.plot(x_line_bpr, noise_g1_line_bpr, 'o', label='Edge wise FD')
ax1.plot(x_line_prc, noise_g1_line_prc, '+', label='PARCS')
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Relative Noise Magnitude (\%)")
fig1.savefig(problem + "_noise_line_amp_g1.pdf", format='pdf')

# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(x_line_ftd, noise_g2_line_ftd, '-', label='Time-Domain')
ax1.plot(x_line_rfd, noise_g2_line_rfd, 'X', label='Cell wise FD')
ax1.plot(x_line_bpr, noise_g2_line_bpr, 'o', label='Edge wise FD')
ax1.plot(x_line_prc, noise_g2_line_prc, '+', label='PARCS')
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Relative Noise Magnitude (\%)")
fig1.savefig(problem + "_noise_line_amp_g2.pdf", format='pdf')

# Print phase_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(x_line_ftd, phase_g1_line_ftd, '-', label='Time-Domain')
ax1.plot(x_line_rfd, phase_g1_line_rfd, 'X',  label='Cell wise FD')
ax1.plot(x_line_bpr, phase_g1_line_bpr, 'o', label='Edge wise FD')
ax1.plot(x_line_prc, phase_g1_line_prc, '+', label='PARCS')
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Phase (deg)")
fig1.savefig(problem + "_noise_line_pha_g1.pdf", format='pdf')

# Print phase_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(x_line_ftd, phase_g2_line_ftd, '-', label='Time-Domain')
ax1.plot(x_line_rfd, phase_g2_line_rfd, 'X',  label='Cell wise FD')
ax1.plot(x_line_bpr, phase_g2_line_bpr, 'o', label='Edge wise FD')
ax1.plot(x_line_prc, phase_g2_line_prc, '+', label='PARCS')
ax1.grid(True)
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Phase (deg)")
ax1.legend(loc='best')
fig1.savefig(problem + "_noise_line_pha_g2.pdf", format='pdf')

#%% 
from scipy.interpolate import interp1d

f1 = interp1d(x_line_ftd, noise_g1_line_ftd, fill_value="extrapolate")
ftd_1 = f1(x_line_bpr)
f2 = interp1d(x_line_ftd, noise_g2_line_ftd, fill_value="extrapolate")
ftd_2 = f2(x_line_bpr)

f1 = interp1d(x_line_rfd, noise_g1_line_rfd, fill_value="extrapolate")
rfd_1 = f1(x_line_bpr)
f2 = interp1d(x_line_rfd, noise_g2_line_rfd, fill_value="extrapolate")
rfd_2 = f2(x_line_bpr)

f1 = interp1d(x_line_prc, noise_g1_line_prc, fill_value="extrapolate")
prc_1 = f1(x_line_bpr)
f2 = interp1d(x_line_prc, noise_g2_line_prc, fill_value="extrapolate")
prc_2 = f2(x_line_bpr)

bpr_1 = noise_g1_line_bpr
bpr_2 = noise_g2_line_bpr

def RAD(test, ref):
    err = list(abs(test-ref)/ref)
    err.remove(max(err))
    err = np.array(err)
    rmd = max(err)*100
    rad = np.mean(err)*100
    return rmd, rad


print(RAD(rfd_1, ftd_1), " ", RAD(rfd_2, ftd_2))
print(RAD(bpr_1, ftd_1), " ", RAD(bpr_2, ftd_2))
print(RAD(prc_1, ftd_1), " ", RAD(prc_2, ftd_2))


