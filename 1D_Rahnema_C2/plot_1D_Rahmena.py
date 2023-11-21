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
import numpy as np
from matplotlib import rcParams
import openmoc

# OPENMOC REFERENCE
# out = openmoc.process.restore_simulation_state(filename='2D_Ranehma_ref30.pkl')
# out =  out['10-27-2023']['17:17:41']
# flxref_g1 = out['FSR scalar fluxes'][:,0]
# flxref_g2 = out['FSR scalar fluxes'][:,1]
# keff_ref = out['keff']

npzfile = np.load('2D_RahnemaC2.npz')
x_ref= npzfile['x']
fluxref_g1 = npzfile['flux_g1']
fluxref_g2 = npzfile['flux_g2']
keff_ref = npzfile['keff']
print(keff_ref)

plt.close('all')

params = {'backend': 'pdf',
#          'font.family': 'serif',
          'font.size': 14,
          'axes.labelsize': 12,
          'legend.fontsize': 10,
          'xtick.labelsize': 13,
          'ytick.labelsize': 13,
          'text.usetex': False,
          'lines.linewidth': 1.2,
          'lines.markersize': 3,
          'lines.markeredgewidth': 1,
          'legend.numpoints': 1, 
          'axes.formatter.useoffset': False,
          'figure.autolayout': True,
          }

rcParams.update(params)

problem = '1D_RahnemaC2'
files = [
     '1D_bwr_conf1_SDP1.out.vtk',
     '1D_bwr_conf1_SDP2.out.vtk',
     '1D_bwr_conf1_SDP3.out.vtk',
     '1D_bwr_conf1_SP1.out.vtk',
     '1D_bwr_conf1_SP3.out.vtk',
     '1D_bwr_conf1_SP5.out.vtk',
     '1D_bwr_conf1_SP7.out.vtk'
    ]
labels = ['SDP1', 'SDP2', 'SDP3', 'SP1', 'SP3', 'SP5', 'SP7' ]

markers =['o', '^', 's', '.', '>', '*', '.' ]
color = ['tab:blue', 'tab:olive', 'tab:purple', 
         'tab:orange', 'tab:green','tab:brown',
         'tab:gray', 'tab:red']
# labels = ['SDP2', 'SP5']

# ref_g7 = [2.9454309e-01, 2.9513030e-01, 2.9630649e-01, 2.9807568e-01, 3.0044526e-01, 3.0342843e-01, 3.0704895e-01, 3.1135077e-01, 3.1641784e-01, 3.2241608e-01, 3.2968477e-01, 3.3894130e-01, 3.5175201e-01, 3.7163838e-01, 4.0672630e-01, 4.7626775e-01, 6.2873747e-01, 9.2668584e-01, 1.1171538e+00, 1.2446732e+00, 1.3344554e+00, 1.4004951e+00, 1.4511081e+00, 1.4913137e+00, 1.5241669e+00, 1.5515378e+00, 1.5745743e+00, 1.5939816e+00, 1.6101915e+00, 1.6234668e+00, 1.6339656e+00, 1.6417818e+00, 1.6469691e+00, 1.6495564e+00]


# Get From VTK
tol = 1e-3



x_lines = []
lines_g1 = []
lines_g2 = []
for f, file in enumerate(files): 
    [x, y, z] = parse_vtk_grid(file)
    stati_g1 = parse_vtk_file(file, "phi_g1_eig_1")
    stati_g2 = parse_vtk_file(file, "phi_g2_eig_1")
    
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
    
    # --------------------------------------------------------------------------
    # Extract Line

    lines_g1.append(np.array(list(stati_g1)))
    lines_g2.append(np.array(list(stati_g2)))
    x_lines.append(list(x))
    
#%% ---------------------------------------------------------------------------   
#  REFERENCE
norm_g1 =  np.mean(lines_g1[-1]) / np.mean(fluxref_g1)
norm_g2 =  np.mean(lines_g2[-1]) / np.mean(fluxref_g2)
x_ref= x_ref - x_ref[0]

#%% ---------------------------------------------------------------------------


# SP3
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(x_lines[0], lines_g1[0], '-', marker=markers[0], c=color[0], label=labels[0])
ax1.plot(x_lines[4], lines_g1[4], '-', marker=markers[4], c=color[4], label=labels[4])
ax1.plot(x_ref, norm_g1 * fluxref_g1, '-',  c=color[7], label='OpenMOC')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Scalar Fast Flux (AU)")
ax1.grid(True)
ax1.legend(loc='upper center')
fig1.savefig(problem + "SP3_g1.pdf", format='pdf')


# SP5
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(x_lines[1], lines_g1[1], '-', marker=markers[1], c=color[1], label=labels[1])
ax1.plot(x_lines[5], lines_g1[5], '-', marker=markers[5], c=color[5], label=labels[5])
ax1.plot(x_ref, norm_g1 * fluxref_g1, '-',  c=color[7], label='OpenMOC')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Scalar Fast Flux (AU)")
ax1.grid(True)
ax1.legend(loc='upper center')
fig1.savefig(problem + "SP5_g1.pdf", format='pdf')


# SP7
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(x_lines[2], lines_g1[2], '-', marker=markers[2], c=color[2], label=labels[2])
ax1.plot(x_lines[6], lines_g1[6], '-', marker=markers[6], c=color[6], label=labels[6])
ax1.plot(x_ref, norm_g1 * fluxref_g1, '-',  c=color[7], label='OpenMOC')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Scalar Fast Flux (AU)")
ax1.grid(True)
ax1.legend(loc='upper center')
fig1.savefig(problem + "SP7_g1.pdf", format='pdf')

# SP3 G2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(x_lines[0], lines_g2[0], '-', marker=markers[0], c=color[0], label=labels[0])
ax1.plot(x_lines[4], lines_g2[4], '-', marker=markers[4], c=color[4], label=labels[4])
ax1.plot(x_ref, norm_g2 * fluxref_g2, '-',  c=color[7], label='OpenMOC')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Scalar Thermal Flux (AU)")
ax1.grid(True)
ax1.legend(loc='upper center')
fig1.savefig(problem + "SP3_g2.pdf", format='pdf')

# SP5 G2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(x_lines[1], lines_g2[1], '-', marker=markers[1], c=color[1], label=labels[1])
ax1.plot(x_lines[5], lines_g2[5], '-', marker=markers[5], c=color[5], label=labels[5])
ax1.plot(x_ref, norm_g2 * fluxref_g2, '-',  c=color[7], label='OpenMOC')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Scalar Thermal Flux (AU)")
ax1.grid(True)
ax1.legend(loc='upper center')
fig1.savefig(problem + "SP5_g2.pdf", format='pdf')


# SP7 G2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(x_lines[2], lines_g1[2], '-', marker=markers[2], c=color[2], label=labels[2])
ax1.plot(x_lines[6], lines_g1[6], '-', marker=markers[6], c=color[6], label=labels[6])
ax1.plot(x_ref, norm_g2 * fluxref_g2, '-',  c=color[7], label='OpenMOC')
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Scalar Thermal Flux (AU)")
ax1.grid(True)
ax1.legend(loc='upper center')
fig1.savefig(problem + "SP7_g2.pdf", format='pdf')






