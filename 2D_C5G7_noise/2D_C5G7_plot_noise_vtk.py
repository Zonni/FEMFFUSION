# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:33:28 2019

@author: zonni
"""
import sys
sys.path.append('../postprocess')
from utils import parse_vtk_file, parse_vtk_grid
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib import rcParams
from matplotlib import cm

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
          'lines.markersize': 5,
          'lines.markeredgewidth': 1,
          'legend.numpoints': 1, 
          'axes.formatter.useoffset': False,
          'figure.autolayout': True,
          }

rcParams.update(params)

file = 'NN1.2/c5g7_ref1_FE2_DSP3.out.vtk'
problem = 'NN1.2/c5g7_ref1_FE2_DSP3'
#cmap="viridis" 
#cmap="inferno"
cmap="viridis"

# Get From VTK
[x, y, z] = parse_vtk_grid(file)
stati_g1 = parse_vtk_file(file, "Static_Flux_g1")
stati_g2 = parse_vtk_file(file, "Static_Flux_g7")
noise_g1 = parse_vtk_file(file, "Noise_g1_Magnitude")
noise_g2 = parse_vtk_file(file, "Noise_g7_Magnitude")
phase_g1 = parse_vtk_file(file, "Noise_g1_Phase")
phase_g2 = parse_vtk_file(file, "Noise_g7_Phase")
n_levels = 150

#%% ---------------------------------------------------------------------------
# Remove repeated data
dict_stati_g1 = {}
dict_stati_g2 = {}
dict_noise_g1 = {}
dict_noise_g2 = {}
dict_phase_g1 = {}
dict_phase_g2 = {}
dict_x = {}
dict_y = {}
dict_z = {}
for p in range(len(noise_g1)):
        dict_stati_g1[(x[p], y[p], z[p])] = stati_g1[p]
        dict_stati_g2[(x[p], y[p], z[p])] = stati_g2[p]
        
        dict_noise_g1[(x[p], y[p], z[p])] = noise_g1[p]
        dict_noise_g2[(x[p], y[p], z[p])] = noise_g2[p]
        
        dict_phase_g1[(x[p], y[p], z[p])] = phase_g1[p]
        dict_phase_g2[(x[p], y[p], z[p])] = phase_g2[p]
        
        dict_x[(x[p], y[p], z[p])] = x[p]
        dict_y[(x[p], y[p], z[p])] = y[p]
        dict_z[(x[p], y[p], z[p])] = z[p]
        
stati_g1 = list(dict_stati_g1.values())
stati_g2 = list(dict_stati_g2.values())   
noise_g1 = list(dict_noise_g1.values())
noise_g2 = list(dict_noise_g2.values())
phase_g1 = list(dict_phase_g1.values())
phase_g2 = list(dict_phase_g2.values())
x = list(dict_x.values())
y = list(dict_y.values())
z = list(dict_z.values())

# -----------------------------------------------------------------------------
# Mean plane
dict_stati_g1 = {}
dict_stati_g2 = {}
dict_noise_g1 = {}
dict_noise_g2 = {}
dict_phase_g1 = {}
dict_phase_g2 = {}
dict_x = {}
dict_y = {}

times = {}
for p in range(len(noise_g1)):
    if (x[p], y[p]) in dict_noise_g1:
        dict_stati_g1[(x[p], y[p])] += stati_g1[p]
        dict_stati_g2[(x[p], y[p])] += stati_g2[p]
        
        dict_noise_g1[(x[p], y[p])] += noise_g1[p]
        dict_noise_g2[(x[p], y[p])] += noise_g2[p]
        dict_phase_g1[(x[p], y[p])] += phase_g1[p]
        dict_phase_g2[(x[p], y[p])] += phase_g2[p]
        times[(x[p], y[p])] += 1 
    else:
        dict_stati_g1[(x[p], y[p])] = stati_g1[p]
        dict_stati_g2[(x[p], y[p])] = stati_g2[p]
        
        dict_noise_g1[(x[p], y[p])] = noise_g1[p]
        dict_noise_g2[(x[p], y[p])] = noise_g2[p]
        dict_phase_g1[(x[p], y[p])] = phase_g1[p]
        dict_phase_g2[(x[p], y[p])] = phase_g2[p]
        
        dict_x[(x[p], y[p])] = x[p]
        dict_y[(x[p], y[p])] = y[p]
        times[(x[p], y[p])] = 1
        
for key in dict_noise_g1:
    dict_stati_g1[key] /= times[key]
    dict_stati_g2[key] /= times[key]
    dict_noise_g1[key] /= times[key]
    dict_noise_g2[key] /= times[key]
    dict_phase_g1[key] /= times[key]
    dict_phase_g2[key] /= times[key]
    
    # Relative noise
    dict_noise_g1[key] /= dict_stati_g1[key]
    dict_noise_g2[key] /= dict_stati_g2[key]
    
    # Percentage
    dict_noise_g1[key] *= 100
    dict_noise_g2[key] *= 100

 
static_g1 = list(dict_stati_g1.values()) 
static_g2 = list(dict_stati_g2.values())       
noise_g1 = list(dict_noise_g1.values()) 
noise_g2 = list(dict_noise_g2.values())
phase_g1 = list(dict_phase_g1.values())
phase_g2 = list(dict_phase_g2.values())
x = list(dict_x.values())
y = list(dict_y.values())



#%% ---------------------------------------------------------------------------
# Debido al bug del components
y2 = y
y = x
x = y2

# Print stati_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
cntr2 = ax1.tricontourf(x, y, stati_g1, levels=n_levels)
ax1.set_aspect('equal')
fig1.colorbar(cntr2, ax=ax1)
ax1.set_ylabel("y (cm)")
ax1.set_xlabel("x (cm)")
ax1.invert_yaxis()
ax1.xaxis.set_ticks_position('top') # the rest is the same
ax1.xaxis.set_label_position('top') 
fig1.tight_layout()
fig1.savefig(problem + "_staticvtk_g1.png", format='png',dpi=600)

#
## Print noise_g7
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
cntr2 = ax1.tricontourf(x, y, stati_g2, levels=n_levels)
ax1.set_aspect('equal')
fig1.colorbar(cntr2, ax=ax1)
ax1.set_ylabel("y (cm)")
ax1.set_xlabel("x (cm)")
ax1.invert_yaxis()
ax1.xaxis.set_ticks_position('top') # the rest is the same
ax1.xaxis.set_label_position('top') 
fig1.tight_layout()
fig1.savefig(problem + "_staticvtk_g7.png", format='png',dpi=600)

## Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
cntr2 = ax1.tricontourf(x, y, noise_g1, levels=n_levels, cmap=cmap)
ax1.set_aspect('equal')
fig1.colorbar(cntr2, ax=ax1)
ax1.set_ylabel("y (cm)")
ax1.set_xlabel("x (cm)")
ax1.invert_yaxis()
ax1.xaxis.set_ticks_position('top') # the rest is the same
ax1.xaxis.set_label_position('top') 
fig1.tight_layout()
fig1.savefig(problem + "_noisevtk_g1_amp.png", format='png',dpi=600)

#
## Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
cntr2 = ax1.tricontourf(x, y, noise_g2, levels=n_levels, cmap=cmap)
ax1.set_aspect('equal')
fig1.colorbar(cntr2, ax=ax1)
ax1.set_ylabel("y (cm)")
ax1.set_xlabel("x (cm)")
ax1.invert_yaxis()
ax1.xaxis.set_ticks_position('top') # the rest is the same
ax1.xaxis.set_label_position('top') 
fig1.tight_layout()
fig1.savefig(problem + "_noisevtk_g7_amp.png", format='png',dpi=600)
#
#
# Print noise_phase_g1
norm = matplotlib.colors.Normalize(vmin=0, vmax=360)
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
cntr2 = ax1.tricontourf(x, y, phase_g1, levels=n_levels, cmap=cmap)
fig1.colorbar(cntr2, ax=ax1)
ax1.set_aspect('equal')
ax1.set_ylabel("y (cm)")
ax1.set_xlabel("x (cm)")
ax1.invert_yaxis()
ax1.xaxis.set_ticks_position('top') # the rest is the same
ax1.xaxis.set_label_position('top') 
fig1.savefig(problem + "_noisevtk_g1_phase.png", format='png', dpi=600)
#
#
# Print noise_phase_g1
norm = matplotlib.colors.Normalize(vmin=0, vmax=360)
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
cntr2 = ax1.tricontourf(x, y, phase_g2, levels=n_levels, cmap=cmap)
fig1.colorbar(cntr2, ax=ax1)
ax1.set_aspect('equal')
ax1.set_ylabel("y (cm)")
ax1.set_xlabel("x (cm)")
ax1.invert_yaxis()
ax1.xaxis.set_ticks_position('top') # the rest is the same
ax1.xaxis.set_label_position('top') 
fig1.savefig(problem + "_noisevtk_g7_phase.png", format='png', dpi=600)
#