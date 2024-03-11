# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:33:28 2019

@author: zonni
"""
from utils import parse_vtk_file, parse_vtk_grid
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib import rcParams

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

file = '../3D_VVER1000_NOISE/3D_VVER1000_FE2.out.vtk'
problem = '3D_VVER1000_NOISE'
#cmap="viridis" 
#cmap="inferno"
cmap="RdBu_r"

# Get From VTK
[x, y, z] = parse_vtk_grid(file)
stati_g1 = parse_vtk_file(file, "Static_Flux_g1")
stati_g2 = parse_vtk_file(file, "Static_Flux_g2")
noise_g1 = parse_vtk_file(file, "Noise_g1_Magnitude")
noise_g2 = parse_vtk_file(file, "Noise_g2_Magnitude")
phase_g1 = parse_vtk_file(file, "Noise_g1_Phase")
phase_g2 = parse_vtk_file(file, "Noise_g2_Phase")
n_levels = 10

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

        
noise_g1 = list(dict_noise_g1.values()) 
noise_g2 = list(dict_noise_g2.values())
phase_g1 = list(dict_phase_g1.values())
phase_g2 = list(dict_phase_g2.values())
x = list(dict_x.values())
y = list(dict_y.values())


#%% ---------------------------------------------------------------------------

# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.tricontour(x, y, noise_g1, levels=n_levels, linewidths=0.5, colors='k')
cntr2 = ax1.tricontourf(x, y, noise_g1, levels=n_levels, cmap=cmap)
ax1.set_aspect('equal')
fig1.colorbar(cntr2, ax=ax1)
ax1.set_ylabel("y (cm)")
ax1.set_xlabel("x (cm)")
ax1.set_title("Relative Fast Noise Magnitude (\%)")
fig1.savefig(problem + "_noisevtk_g1_amp.pdf", format='pdf')


# Print noise_g2
fig2 = plt.figure()
ax2 = fig2.add_subplot(1, 1, 1)
ax2.tricontour(x, y, noise_g2, levels=n_levels, linewidths=0.5, colors='k')
cntr2 = ax2.tricontourf(x, y, noise_g2, levels=n_levels, cmap=cmap)
fig2.colorbar(cntr2, ax=ax2)
ax2.set_aspect('equal')
ax2.set_ylabel("y (cm)")
ax2.set_xlabel("x (cm)")
ax2.set_title("Relative Thermal Noise Magnitude (\%)")
fig2.savefig(problem + "_noisevtk_g2_amp.pdf", format='pdf')


# Print noise_phase_g1
norm = matplotlib.colors.Normalize(vmin=0, vmax=360)
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.tricontour(x, y, phase_g1, levels=n_levels, linewidths=0.5, colors='k')
cntr2 = ax1.tricontourf(x, y, phase_g1, levels=n_levels, cmap=cmap)
fig1.colorbar(cntr2, ax=ax1)
ax1.set_aspect('equal')
ax1.set_ylabel("y (cm)")
ax1.set_xlabel("x (cm)")
ax1.set_title("Fast Noise Phase (deg)")
fig1.savefig(problem + "_noisevtk_g1_phase.pdf", format='pdf')



# Print noise_phase_g2
fig2 = plt.figure()
ax2 = fig2.add_subplot(1, 1, 1)
ax2.tricontour(x, y, phase_g2, levels=n_levels, linewidths=0.5, colors='k')
cntr2 = ax2.tricontourf(x, y, phase_g2, levels=n_levels, cmap=cmap)
fig2.colorbar(cntr2, ax=ax2)
ax2.set_aspect('equal')
ax2.set_ylabel("y (cm)")
ax2.set_xlabel("x (cm)")
ax2.set_title("Thermal Noise Phase (deg)")
fig2.savefig(problem + "_noisevtk_g2_phase.pdf", format='pdf')