#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Antoni Vidal """

import numpy as np
from utils import parse_line_fluxes_parcs_2D, parse_fluxes_parcs_static
from utils import parse_vtk_grid, parse_vtk_file, get_eigenvalues
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy.io

plt.close('all')
color = ['#1f77b4', '#d62728', '#ff7f0e', '#2ca02c', '#ff9896',
         '#c5b0d5', '#8c564b', '#e377c2', '#f7b6d2', '#7f7f7f',
         '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
lines = ['-', '--']
params = {'backend': 'pdf',
          'font.family': 'serif',
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
cmap="RdBu_r"

problem = '2D_test_vibration'

if problem == '2D_test_vibration':
    folder = "../2D_test_vibration/"
    
    # PARCS
    delta_t_parcs = 1e-2
    parcs_outfile = folder + '2D_test.out.parcs'
    n_line = 30
    looking_freq = 1.0
    keff_prc = 0.875889
    
    
    # FEMFFUSION - CELL-WISE
    y_line = 30.0
    tol = 1e-5
    fem_outfile = folder + '2D_test_ref.out.vtk'
    fem_out = folder + '2D_test_ref.out'
    
    # FEMFFUSION - BORDERS
    fem_FAv_outfile = folder + '2D_test.out.vtk'
    fem_FAv_out = folder + '2D_test.out'
    
    # CORESIM - FREQ
    DX = 0.1
    ny_line = 60
    cor_outfile = folder + 'output/RESULTS.mat'

# =============================================================================
print('PARSE PARCS .OUT FILES')
# =============================================================================
noise_fs_line = 0
noise_th_line = 0
x_line = 0
sta_flux_line_fs = 0
sta_flux_line_th = 0

static_fs, static_th = parse_fluxes_parcs_static(parcs_outfile)
norm = max(static_fs)
#x_prc, y_parcs = parse_parcs_geometry(parcs_outfile)
x_prc, flux_prc_fs, flux_prc_th = parse_line_fluxes_parcs_2D(parcs_outfile,
                                                             n_line)

# Remove last Parcs time step as we are counting 0
flux_prc_fs = flux_prc_fs[0:-1]
flux_prc_th = flux_prc_th[0:-1]

n_steps = len(flux_prc_th)
n_nodes = len(flux_prc_th[0])

# Normalize
for st in range(n_steps) :
    flux_prc_fs[st] /= norm
    flux_prc_th[st] /= norm

# Compute time and space dependent neutron noise
static_flux_prc_fs =  flux_prc_fs[0].copy()
static_flux_prc_th =  flux_prc_th[0].copy()
for st in range(n_steps) :
    flux_prc_fs[st] -= static_flux_prc_fs
    flux_prc_th[st] -= static_flux_prc_th

# Traspose
flux_prc_fs2 = np.zeros([n_nodes, n_steps]) 
flux_prc_th2 = np.zeros([n_nodes, n_steps])

for t in range(n_steps):
    for node in range(n_nodes):
        flux_prc_fs2[node][t] = flux_prc_fs[t][node]
        flux_prc_th2[node][t] = flux_prc_th[t][node]

freq = np.fft.rfftfreq(n_steps, d=delta_t_parcs)
ft_fs_parcs = np.fft.rfft(flux_prc_fs2) * 2.0/ n_steps 
ft_th_parcs = np.fft.rfft(flux_prc_th2) * 2.0/ n_steps

# We cut at looking_freq Hz
cut_freq = int (looking_freq * n_steps * delta_t_parcs)
assert(freq[cut_freq] == looking_freq)
esp_fs_parcs = np.zeros([n_nodes], dtype='cfloat')
esp_th_parcs = np.zeros([n_nodes], dtype='cfloat')
for node in range(n_nodes):
    esp_fs_parcs[node] = ft_fs_parcs[node][cut_freq]
    esp_th_parcs[node] = ft_th_parcs[node][cut_freq]
    
print(' Fast Noise Max PARCS: ',  max(abs(esp_fs_parcs)))
print(' Ther Noise Max PARCS: ',  max(abs(esp_th_parcs)))
print(' keff_parcs: ', round(keff_prc, 5))


## ============================================================================
print('FEMFFUSION FREQ')
## ============================================================================
# Get From VTK
[x, y, z] = parse_vtk_grid(fem_outfile)
n_nodes = len(x)
flux_g1 = parse_vtk_file(fem_outfile, "Static_Flux_g1")
flux_g2 = parse_vtk_file(fem_outfile, "Static_Flux_g2")
noise_g1 = parse_vtk_file(fem_outfile, "Noise_g1_Magnitude")
noise_g2 = parse_vtk_file(fem_outfile, "Noise_g2_Magnitude")
phase_g1 = parse_vtk_file(fem_outfile, "Noise_g1_Phase")
phase_g2 = parse_vtk_file(fem_outfile, "Noise_g2_Phase")

# GET LINE FEMFFUSION
x_line = []
sta_flux_line_fs = []
sta_flux_line_th = []
mag_noise_fs_line = []
mag_noise_th_line = []
pha_noise_fs_line = []
pha_noise_th_line = []

for node in range(n_nodes):
    if (abs(y[node] - y_line) < tol):
        x_line.append(x[node])
        sta_flux_line_fs.append(flux_g1[node])
        sta_flux_line_th.append(flux_g2[node])
        mag_noise_fs_line.append(noise_g1[node])
        mag_noise_th_line.append(noise_g2[node])
        pha_noise_fs_line.append(phase_g1[node]) 
        pha_noise_th_line.append(phase_g2[node])
   
sta_flux_line_fs = np.array(sta_flux_line_fs)  
sta_flux_line_th = np.array(sta_flux_line_th)  
mag_noise_fs_line = np.array(mag_noise_fs_line)    
mag_noise_th_line = np.array(mag_noise_th_line)
pha_noise_fs_line = np.array(pha_noise_fs_line)
pha_noise_th_line = np.array(pha_noise_th_line)

print(' Fast Noise Max FEMFFUSION: ',  max(mag_noise_fs_line))
print(' Ther Noise Max FEMFFUSION: ',  max(mag_noise_th_line))
#pha_noise_fs_line[302] = 85.5012
#pha_noise_fs_line[262] = 85.5012
#pha_noise_th_line[302] = 85.5012
#pha_noise_th_line[262] = 85.5012
[x_line,
 sta_flux_line_fs,
 sta_flux_line_th,
 mag_noise_fs_line,
 mag_noise_th_line,
 pha_noise_fs_line,
 pha_noise_th_line] = zip(*sorted(zip(
                   x_line,
                   sta_flux_line_fs,
                   sta_flux_line_th,
                   mag_noise_fs_line,
                   mag_noise_th_line,
                   pha_noise_fs_line,
                   pha_noise_th_line)))

keff_fem = round(get_eigenvalues(fem_out)[0], 5)
print(' keff_fem_cellwise: ', keff_fem)


## ============================================================================
print('FEMFFUSION BORDERS')
## ============================================================================
# Get From VTK
[x, y, z] = parse_vtk_grid(fem_FAv_outfile)
n_nodes = len(x)
flux_g1 = parse_vtk_file(fem_FAv_outfile, "Static_Flux_g1")
flux_g2 = parse_vtk_file(fem_FAv_outfile, "Static_Flux_g2")
noise_g1 = parse_vtk_file(fem_FAv_outfile, "Noise_g1_Magnitude")
noise_g2 = parse_vtk_file(fem_FAv_outfile, "Noise_g2_Magnitude")
phase_g1 = parse_vtk_file(fem_FAv_outfile, "Noise_g1_Phase")
phase_g2 = parse_vtk_file(fem_FAv_outfile, "Noise_g2_Phase")

# GET LINE FEMFFUSION
x_line = []
sta_flux_line_fs = []
sta_flux_line_th = []
mag_noise_fs_line = []
mag_noise_th_line = []
pha_noise_fs_line = []
pha_noise_th_line = []

for node in range(n_nodes):
    if (abs(y[node] - y_line) < tol):
        x_line.append(x[node])
        sta_flux_line_fs.append(flux_g1[node])
        sta_flux_line_th.append(flux_g2[node])
        mag_noise_fs_line.append(noise_g1[node])
        mag_noise_th_line.append(noise_g2[node])
        pha_noise_fs_line.append(phase_g1[node]) 
        pha_noise_th_line.append(phase_g2[node])
   
sta_flux_line_fs = np.array(sta_flux_line_fs)  
sta_flux_line_th = np.array(sta_flux_line_th)  
mag_noise_fs_line = np.array(mag_noise_fs_line)    
mag_noise_th_line = np.array(mag_noise_th_line)
pha_noise_fs_line = np.array(pha_noise_fs_line)
pha_noise_th_line = np.array(pha_noise_th_line)

print(' Fast Noise Max FEMFFUSION: ',  max(mag_noise_fs_line))
print(' Ther Noise Max FEMFFUSION: ',  max(mag_noise_th_line))
#pha_noise_fs_line[302] = 85.5012
#pha_noise_fs_line[262] = 85.5012
#pha_noise_th_line[302] = 85.5012
#pha_noise_th_line[262] = 85.5012
[x_line,
 sta_flux_line_fs,
 sta_flux_line_th,
 mag_noise_fs_line,
 mag_noise_th_line,
 pha_noise_fs_line,
 pha_noise_th_line] = zip(*sorted(zip(
                   x_line,
                   sta_flux_line_fs,
                   sta_flux_line_th,
                   mag_noise_fs_line,
                   mag_noise_th_line,
                   pha_noise_fs_line,
                   pha_noise_th_line)))

keff_fem = round(get_eigenvalues(fem_FAv_out)[0], 5)
print(' keff_fem_FAv: ', keff_fem)

## ============================================================================
print('CORESIM')
## ============================================================================


mat = scipy.io.loadmat(cor_outfile)
FLX1_line = mat['FLX1'][:, 59, 1];
FLX2_line = mat['FLX2'][:, 59, 1];
dFLX1_line = mat['dFLX1'][:, 59, 1];
dFLX2_line = mat['dFLX2'][:, 59, 1];
keff_cor = round(mat['keff'][0][0], 5)


n_nodes = len(FLX1_line)
x_cor = np.linspace(DX/2, (60.-DX/2), n_nodes)

print(' keff_cor: ', keff_cor)
print(' Fast Noise Max CORESIM: ',  max(abs(dFLX1_line)))
print(' Ther Noise Max CORESIM: ',  max(abs(dFLX2_line)))

## ============================================================================
print('PLOT')
## ============================================================================


lab = ["Fast Noise", "Thermal Noise"]
figsize = (9,6)

## Static Flux
fig0 = plt.figure(figsize=figsize)
ax0 = fig0.add_subplot(1, 1, 1)
ax0.plot(x_line, sta_flux_line_fs, marker='x', linestyle='-', c=color[0], label="FEM Fast Flux")
ax0.plot(x_line, sta_flux_line_th, marker='x', linestyle='-', c=color[1], label="FEM Thermal Flux")
ax0.grid(True)
ax0.set_ylabel("Static Neutron Flux")
ax0.set_xlabel("x (cm)")
ax0.plot(x_prc, static_flux_prc_fs, c=color[2], linestyle='', marker='o', label="PARCS Fast Flux")
ax0.plot(x_prc, static_flux_prc_th, c=color[3], linestyle='', marker='o', label="PARCS Thermal Flux")
ax0.plot(x_cor, FLX1_line, c=color[4], linestyle='--', marker='', label="CORESIM Fast Flux")
ax0.plot(x_cor, FLX2_line, c=color[5], linestyle='--', marker='', label="CORESIM Thermal Flux")
ax0.legend()
fig0.savefig(problem + "_staticflux.pdf", format='pdf')

#
fig2 = plt.figure(figsize=figsize)
ax2 = fig2.add_subplot(1, 1, 1)
ax2.plot(x_line, mag_noise_fs_line, marker='x', c=color[0], linestyle='-', label='FEM ' + lab[0])
ax2.plot(x_line, mag_noise_th_line, marker='x', c=color[1], linestyle='-', label='FEM ' + lab[1])

ax2.plot(x_prc, abs(esp_fs_parcs), c=color[2], linestyle='', marker='o', label='PARCS ' + lab[0])
ax2.plot(x_prc, abs(esp_th_parcs), c=color[3], linestyle='', marker='o', label='PARCS ' + lab[1])

ax2.plot(x_cor, abs(dFLX1_line), c=color[4], linestyle='--', marker='', label='CORESIM ' + lab[0])
ax2.plot(x_cor, abs(dFLX2_line), c=color[5], linestyle='--', marker='', label='CORESIM ' + lab[1])

ax2.grid(True)
ax2.set_ylabel("Noise Amplitude ")
ax2.set_xlabel("x (cm)")
ax2.legend(loc='best')
fig2.savefig(problem + "_fourier_amp_parcs.pdf", format='pdf')

#
fig3 = plt.figure(figsize=figsize)
ax3 = fig3.add_subplot(1, 1, 1)
ax3.plot(x_line, pha_noise_fs_line, c=color[0], linestyle='-', marker='x', label='FEM ' + lab[0])
ax3.plot(x_line, pha_noise_th_line, c=color[1], linestyle='-', marker='x', label='FEM ' + lab[1])

ax3.plot(x_prc, np.angle(esp_fs_parcs, deg=True), c=color[2], linestyle='',  marker='o', label='PARCS ' + lab[0])
ax3.plot(x_prc, np.angle(esp_th_parcs, deg=True), c=color[3], linestyle='',  marker='o', label='PARCS ' + lab[1])

ax3.plot(x_cor, np.angle(dFLX1_line, deg=True), c=color[4], linestyle='--', marker='', label='CORESIM ' + lab[0])
ax3.plot(x_cor, np.angle(dFLX2_line, deg=True), c=color[5], linestyle='--', marker='', label='CORESIM ' + lab[1])

ax3.grid(True)
ax3.set_ylabel("Noise Phase [deg]")
ax3.set_xlabel("x (cm)")
ax3.legend()
fig3.savefig(problem + "_fourier_phase_parcs.pdf", format='pdf')

#variables =  {'x_py': x, 'py_fs': esp_fs, 'py_th': esp_th}
#sio.savemat('CORESIM/output/time_1D_fft.mat', variables)
