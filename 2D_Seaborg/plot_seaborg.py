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

file1 = '../2D_Seaborg/seaborg_diff.out.vtk'
file2 = '../2D_Seaborg/seaborg_sp3.out.vtk'
file3 = '../2D_Seaborg/seaborg_sp5.out.vtk'

files = [file1, file2, file3]
labels = ['Diffusion', 'SP3', 'SP5']

problem = '2D_Seaborg'
#cmap="viridis" 
#cmap="inferno"
cmap="RdBu_r"



xlist     = []
static_g1 = []
static_g2 = []
for i, file in enumerate(files): 
    # Get From VTK
    [x, y, z] = parse_vtk_grid(file)
    stati_g1 = parse_vtk_file(file, "phi_g1_eig_1")
    stati_g2 = parse_vtk_file(file, "phi_g2_eig_1")
    n_levels = 125
    
    #%% ---------------------------------------------------------------------------
    # Remove repeated data
    dict_stati_g1 = {}
    dict_stati_g2 = {}
    dict_x = {}
    dict_y = {}
    dict_z = {}
    for p in range(len(stati_g1)):
            dict_stati_g1[(x[p], y[p], z[p])] = stati_g1[p]
            dict_stati_g2[(x[p], y[p], z[p])] = stati_g2[p]
            
            dict_x[(x[p], y[p], z[p])] = x[p]
            dict_y[(x[p], y[p], z[p])] = y[p]
            dict_z[(x[p], y[p], z[p])] = z[p]
            
    stati_g1 = list(dict_stati_g1.values())
    stati_g2 = list(dict_stati_g2.values())   
    x = list(dict_x.values())
    y = list(dict_y.values())
    z = list(dict_z.values())

    #%% ---------------------------------------------------------------------------
    # GET 
    dict_stati_g1 = {}
    dict_stati_g2 = {}
    dict_x = {}
    dict_y = {}
    dict_z = {}
    for p in range(len(stati_g1)):
        if (abs(y[p]) < 0.1):
            dict_stati_g1[(x[p], y[p], z[p])] = stati_g1[p]
            dict_stati_g2[(x[p], y[p], z[p])] = stati_g2[p]
            
            dict_x[(x[p], y[p], z[p])] = x[p]
            dict_y[(x[p], y[p], z[p])] = y[p]
            dict_z[(x[p], y[p], z[p])] = z[p]
            
    stati_g1 = list(dict_stati_g1.values())
    stati_g2 = list(dict_stati_g2.values())   
    x = list(dict_x.values())
    y = list(dict_y.values())
    z = list(dict_z.values())
    
    x, stati_g1, stati_g2  = zip(*sorted(zip(x, stati_g1, stati_g2)))
    
    static_g1.append(stati_g1)
    static_g2.append(stati_g2)
    xlist.append(x)
    
    
#norm = min(stati_g1)
#stati_g1 /= norm
#stati_g2 /= norm

# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(len(files)):
    ax1.plot(xlist[i], static_g1[i], label=labels[i])
ax1.legend()
ax1.set_ylabel("Fast Flux")
ax1.set_xlabel("x (cm)")
ax1.grid(True)
fig1.savefig(problem + "_static_2.pdf", format='pdf')

# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(len(files)):
    ax1.plot(xlist[i], static_g2[i], label=labels[i])
ax1.legend()
ax1.set_ylabel("Thermal Flux")
ax1.set_xlabel("x (cm)")
ax1.grid(True)
fig1.savefig(problem + "_static_1.pdf", format='pdf')
