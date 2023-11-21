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


files = ['P3_SDP1_Nazari.out.vtk',
         'P3_SDP2_Nazari.out.vtk',
         'P3_SDP3_Nazari.out.vtk',
         'P3_SP1.out.vtk',
         'P3_SP3.out.vtk',
         'P3_SP5.out.vtk',
         'P3_SP7.out.vtk',
         'P3_SDP1.out.vtk',
         'P3_SDP2.out.vtk',
         'P3_SDP3.out.vtk' ]
         
labels = ['SDP1_Naz', 'SDP2_Naz', 'SDP3_Naz', 'SP1', 'SP3', 'SP5', 'SP7', 'SDP1', 'SDP2', 'SDP3']
markers =['o', '^', 's', '.', '>', '*', '.-', '+', 'v', '<']
problem = 'Nazari_P3'

# Get From VTK
y_line = 8.967
norm = 10.


point_a = (0, y_line, 0)
point_b = (10.0, y_line, 0)


lines = []
x_lines = []
for f, file in enumerate(files): 
    x, y = parse_vtk_over_line(file, "phi_g1_eig_1", point_a, point_b, resolution=50)
    x_lines.append(x)
    lines.append(y)
#%% --------------------
#%% ---------------------------------------------------------------------------

# Print phi_g1 line
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for f in range(len((files))): 
    ax1.plot(x_lines[f], lines[f], markers[f], label=labels[f])

ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Scalar Flux (AU)")
ax1.set_title("y=" + str(y_line) + " cm")
ax1.grid(True)
ax1.legend()
fig1.savefig(problem + "nazari.pdf", format='pdf')





