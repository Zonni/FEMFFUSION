# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:33:28 2019

@author: zonni
"""
import sys
sys.path.append('../postprocess')
from utils import parse_vtk_over_line
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
from matplotlib import rcParams
import s16

plt.close('all')

npzfile = np.load('P3.npz')
x_ref= npzfile['x']
fluxref_g1 = npzfile['flux_g1']
keff_ref = npzfile['keff']
print(keff_ref)


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


files = ['P3_SDP3.out.vtk',
         'P3_SDP3_Nazari.out.vtk',
         'P3_SP7.out.vtk']
         
labels =  ['SDP3', 'SDP3 (Nazari et al., 2022)', 'SP7', 'OpenMOC']
markers =['o', '^',  '.', '-']
problem = 'Nazari_P3_SDP3'

# Get From VTK
y_line = 4.5
norm = 10. * 0.3735/0.35


point_a = (0, y_line, 0)
point_b = (10.0, y_line, 0)



lines = []
x_lines = []
for  file in files: 
    x, y = parse_vtk_over_line(file, "phi_g1_eig_1", point_a, point_b, resolution=50)
    x_lines.append(x)
    lines.append(y)
    
# x, y = parse_vtk_over_line(file, "phi_g1_eig_1", point_a, point_b, resolution=500)
# x_lines[-1] = x
# lines[-1] = y    
#%% ---------------------------------------------------------------------------
norm_g1 =  np.mean(lines[0]) / np.mean(fluxref_g1) / norm
x_ref= x_ref - x_ref[0]

#%% ---------------------------------------------------------------------------
# Print phi_g1 line
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for f in range(len((files))): 
    ax1.plot(x_lines[f], lines[f]/norm, markers[f], label=labels[f])

# ax1.plot(s16.x_ref, s16.y_ref, markers[f+1], label=labels[f+1])
ax1.plot(x_ref, norm_g1 * fluxref_g1, markers[f+1], label=labels[f+1])

ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Scalar Flux (AU)")
ax1.set_title("y=" + str(y_line) + " cm", fontsize=10)
ax1.grid(True)
ax1.legend()
fig1.savefig(problem +".pdf", format='pdf')







