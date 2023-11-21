# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 10:31:25 2018

@author: zonni
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
plt.close('all')
from matplotlib import style

style.use('seaborn-darkgrid')
 


cmap = mpl.cm.get_cmap('viridis')
fig = plt.figure()
ax1 = fig.add_subplot('111', projection = '3d')

n_rows = 2
n_cols = 2
n_bars = 4

xpos = np.array([1, 1, 2, 2]) - 0.5
ypos = np.array([1, 2, 1, 2]) - 0.5
zpos = np.zeros(n_rows * n_cols);

x = np.ones(n_bars)
y = np.ones(n_bars)
z = [0.6, 0.8, 0.3, 0.4]
transparency = 0.8

color_norm = mpl.colors.Normalize(vmin=min(z), vmax=max(z))
cl = [list(cmap(color_norm(val))[:3]) + [transparency] for val in z]
ax1.bar3d(xpos, ypos, zpos, x, y, z, color=cl, zsort='average')

plt.xticks(np.arange(1, n_rows+1))
plt.yticks(np.arange(1, n_cols+1))

# fake up the array of the scalar mappable
sm = plt.cm.ScalarMappable(cmap=cmap, norm=color_norm)
sm._A = []
plt.colorbar(sm)

ax1.set_xlabel('rows')
ax1.set_ylabel('cols')
ax1.set_zlabel('z')

plt.show()