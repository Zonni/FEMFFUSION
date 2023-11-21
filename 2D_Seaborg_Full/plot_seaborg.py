# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:33:28 2019

@author: zonni
"""
import sys
sys.path.insert(1, '../postprocess')
from  utils import parse_vtk_file, parse_vtk_grid
import matplotlib.pyplot as plt
# import matplotlib.colors
from matplotlib import rcParams
import numpy as np

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

file1 = 'seaborg_SP1.out.vtk'
file2 = 'seaborg_SP3.out.vtk'
file3 = 'seaborg_SP3.out.vtk'

files = [file1, file2, file3]
labels = ['Diffusion', 'SP3', 'SP5']
styles = ['None', 'None', 'None']

# files = [file1]
# labels = ['Diffusion']



problem = '2D_Seaborg'
#cmap="viridis" 
#cmap="inferno"
cmap="RdBu_r"



xlist     = []
ylist     = []
static_g1 = []
static_g2 = []
static_g1x = []
static_g2x = []
static_g1y = []
static_g2y = []
for i, file in enumerate(files): 
    # Get From VTK
    [x, y, z] = parse_vtk_grid(file)
    stati_g1 = parse_vtk_file(file, "phi_g1_eig_1")
    stati_g2 = parse_vtk_file(file, "phi_g8_eig_1")

    
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
            
    stati_g1l = list(dict_stati_g1.values())
    stati_g2l = list(dict_stati_g2.values())   
    xl = list(dict_x.values())
    yl = list(dict_y.values())
    zl = list(dict_z.values())
    
    xl, stati_g1l, stati_g2l  = zip(*sorted(zip(xl, stati_g1l, stati_g2l)))
    
    static_g1x.append(stati_g1l)
    static_g2x.append(stati_g2l)
    xlist.append(xl)
    
    #%% ---------------------------------------------------------------------------
    # GET 
    dict_stati_g1 = {}
    dict_stati_g2 = {}
    dict_x = {}
    dict_y = {}
    dict_z = {}
    for p in range(len(stati_g1)):
        if (abs(x[p]) < 0.1):
            dict_stati_g1[(x[p], y[p], z[p])] = stati_g1[p]
            dict_stati_g2[(x[p], y[p], z[p])] = stati_g2[p]
            
            dict_x[(x[p], y[p], z[p])] = x[p]
            dict_y[(x[p], y[p], z[p])] = y[p]
            dict_z[(x[p], y[p], z[p])] = z[p]
            
    stati_g1l = list(dict_stati_g1.values())
    stati_g2l = list(dict_stati_g2.values())   
    xl = list(dict_x.values())
    yl = list(dict_y.values())
    zl = list(dict_z.values())
    
    yl, stati_g1l, stati_g2l  = zip(*sorted(zip(yl, stati_g1l, stati_g2l)))
    
    static_g1y.append(stati_g1l)
    static_g2y.append(stati_g2l)
    ylist.append(yl)
    
    
#norm = min(stati_g1)
#stati_g1 /= norm
#stati_g2 /= norm

# # Print noise_g1
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(1, 1, 1)
# for i in range(len(files)):
#     ax1.plot(xlist[i], static_g1[i], label=labels[i])
# ax1.legend()
# ax1.set_ylabel("Fast Flux")
# ax1.set_xlabel("x (cm)")
# ax1.grid(True)
# fig1.savefig(problem + "_static_2.pdf", format='pdf')

# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(len(files)):
    ax1.plot(xlist[i], static_g2x[i], label=labels[i], linestyle=styles[i],  marker='.')
ax1.legend()
ax1.set_ylabel("Thermal Flux g8")
ax1.set_xlabel("x (cm)")
ax1.grid(True)
fig1.savefig(problem + "_static_7x.pdf", format='pdf')


# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for i in range(len(files)):
    ax1.plot(ylist[i], static_g2y[i], label=labels[i],  marker='.', linestyle=styles[i])
ax1.legend()
ax1.set_ylabel("Thermal Flux g8")
ax1.set_xlabel("y (cm)")
ax1.grid(True)
fig1.savefig(problem + "_static_7y.pdf", format='pdf')


np.savetxt('y.txt', ylist[0])
np.savetxt('dif_y.txt', static_g2y[0])
np.savetxt('sp3_y.txt', static_g2y[1])
np.savetxt('sp5_y.txt', static_g2y[2])

np.savetxt('x.txt', xlist[0])
np.savetxt('dif_x.txt', static_g2x[0])
np.savetxt('sp3_x.txt', static_g2x[1])
np.savetxt('sp5_x.txt', static_g2x[2])



