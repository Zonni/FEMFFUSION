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
from pandas import read_csv
 

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

problem = 'Seaborg'
 
file_1 = "sp1_y.csv"
# using genfromtxt()

d = read_csv("sp1_y.csv")
sp1_y = np.array(d['phi_g8_eig_1'].values)
y = d['Points:1'].values
sp1_y /= max(sp1_y)

d = read_csv("sp3_y.csv")
sp3_y = np.array(d['phi_g8_eig_1'].values)
y = d['Points:1'].values
sp3_y /= max(sp3_y)

d = read_csv("sp5_y.csv")
sp5_y = np.array(d['phi_g8_eig_1'].values)
y = d['Points:1'].values
sp5_y /= max(sp5_y)


d = read_csv("sp1_x.csv")
sp1_x = np.array(d['phi_g8_eig_1'].values)
x = d['Points:0'].values
sp1_x /= max(sp1_x)

d = read_csv("sp3_x.csv")
sp3_x = np.array(d['phi_g8_eig_1'].values)
x = d['Points:0'].values
sp3_x /= max(sp3_x)

d = read_csv("sp5_x.csv")
sp5_x = np.array(d['phi_g8_eig_1'].values)
x = d['Points:0'].values
sp5_x /= max(sp5_x)




# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(x, sp1_x)
ax1.plot(x, sp3_x)
ax1.plot(x, sp5_x)
ax1.set_ylabel("Thermal Flux g8")
ax1.set_xlabel("x (cm)")
ax1.grid(True)
fig1.savefig(problem + "_static_7x.pdf", format='pdf')


# Print noise_g2
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(y, sp1_y)
ax1.plot(y, sp3_y)
ax1.plot(y, sp5_y)
ax1.set_ylabel("Thermal Flux g8")
ax1.set_xlabel("y (cm)")
ax1.grid(True)
fig1.savefig(problem + "_static_7y.pdf", format='pdf')





np.savetxt('y.txt', y)
np.savetxt('sp1_y.txt', sp1_y)
np.savetxt('sp3_y.txt', sp3_y)
np.savetxt('sp5_y.txt', sp5_y)

np.savetxt('x.txt', x)
np.savetxt('sp1_x.txt', sp1_x)
np.savetxt('sp3_x.txt', sp3_x)
np.savetxt('sp5_x.txt', sp5_x)



