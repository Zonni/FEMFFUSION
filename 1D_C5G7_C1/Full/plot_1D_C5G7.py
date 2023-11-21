# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:33:28 2019

@author: zonni
"""
import sys
sys.path.append('../../postprocess')
from utils import parse_vtk_file, parse_vtk_grid
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
from matplotlib import rcParams

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


files = ['het_SDP1.out.vtk',
         # 'het_SDP2.out.vtk',
         # 'het_SDP3.out.vtk',
         # 'het_SP1.out.vtk',
         'het_SP3.out.vtk',
         # 'het_SP5.out.vtk',
         'het_SP7.out.vtk']
         
# labels = ['SDP1', 'SDP2', 'SDP3', 'SP1', 'SP3', 'SP5', 'SP7' ]
# markers =['o', '^', 's', '.', '>', '*', '.' ]
lines = []
problem = '1D_C5G7'
ref_g7 = [2.9454309e-01, 2.9513030e-01, 2.9630649e-01, 2.9807568e-01, 3.0044526e-01, 3.0342843e-01, 3.0704895e-01, 3.1135077e-01, 3.1641784e-01, 3.2241608e-01, 3.2968477e-01, 3.3894130e-01, 3.5175201e-01, 3.7163838e-01, 4.0672630e-01, 4.7626775e-01, 6.2873747e-01, 9.2668584e-01, 1.1171538e+00, 1.2446732e+00, 1.3344554e+00, 1.4004951e+00, 1.4511081e+00, 1.4913137e+00, 1.5241669e+00, 1.5515378e+00, 1.5745743e+00, 1.5939816e+00, 1.6101915e+00, 1.6234668e+00, 1.6339656e+00, 1.6417818e+00, 1.6469691e+00, 1.6495564e+00]

x_ref = np.zeros(34)
for i in range (34):
    x_ref[i] = 1.26 * i+ 1.26/2

# Get From VTK
tol = 1e-3

labels = ['SDP1', 'SP3', 'SP7' ]
markers =[]

x_lines = []
for f, file in enumerate(files): 
    [x, y, z] = parse_vtk_grid(file)
    stati_g1 = parse_vtk_file(file, "phi_g7_eig_1")
    
    #%% ---------------------------------------------------------------------------
    # Remove repeated data
    dict_stati_g1 = {}
    dict_x = {}
    dict_y = {}
    dict_z = {}
    for p in range(len(stati_g1)):
            dict_stati_g1[(x[p], y[p], z[p])] = stati_g1[p]
            
            dict_x[(x[p], y[p], z[p])] = x[p]
            dict_y[(x[p], y[p], z[p])] = y[p]
            dict_z[(x[p], y[p], z[p])] = z[p]
            
    stati_g1 = list(dict_stati_g1.values())
    x = list(dict_x.values())
    y = list(dict_y.values())
    z = list(dict_z.values())
    
    # -----------------------------------------------------------------------------
    # Extract Line

          

    lines.append(np.array(list(stati_g1)))
    x_lines.append(list(x))

#%% ---------------------------------------------------------------------------

# Print phi_g1 line
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for f in range(len((files))): 
    ax1.plot(x_lines[f], lines[f], '-', label=labels[f])
ax1.plot(x_ref, ref_g7, '-',  label='S96')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Scalar Flux (AU)")
# ax1.set_title("y = 9.0 cm")
ax1.grid(True)
ax1.legend()
fig1.savefig(problem + "nazari_45.pdf", format='pdf')



