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
          'lines.markersize': 0,
          'lines.markeredgewidth': 0,
          'legend.numpoints': 1, 
          'axes.formatter.useoffset': False,
          'figure.autolayout': True,
          }

rcParams.update(params)

g = 7
files = ['c5g7_SDP1.out.vtk',
         'c5g7_SDP2.out.vtk',
         'c5g7_SDP3.out.vtk' ]
         
labels =  ['SDP1', 'SDP2', 'SDP3']
markers =['-', '--', ':']
problem = 'C5G7_SDP1'

# Get From VTK
norm = 10. * 0.3735/0.35


point_a = (0, 0, 0)
point_b = (17*1.26, 17*1.26, 0)

# point_a = (0, 2*17*1.26, 0)
# point_b = (17*1.26, 17*1.26, 0)


lines = []
x_lines = []
for f, file in enumerate(files): 
    x, y = parse_vtk_over_line(file, "phi_g"+str(g)+ "_eig_1", point_a, point_b, resolution=1500)
    x_lines.append(x)
    lines.append(y)
#%% ---------------------------------------------------------------------------

# Print phi_g1 line
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for f in range(len((files))): 
    ax1.plot(x_lines[f], lines[f]/norm, markers[f], label=labels[f])

ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Scalar Flux (AU)")
ax1.set_title("Energy group " + str(g), fontsize=10)
ax1.grid(True)
ax1.legend()
fig1.savefig(problem +"_g" + str(g) +".pdf", format='pdf')





