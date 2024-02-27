#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 12:37:49 2019

@author: zonni
"""

###############################################################################
###############################################################################
###############################################################################
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
plt.close('all')

color = ['#1f77b4', '#d62728', '#2ca02c', '#ff7f0e', '#ff9896',
         '#c5b0d5', '#8c564b', '#e377c2', '#f7b6d2', '#7f7f7f',
         '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
lines = ['-', ':', '--', '-'] 

plt.style.use('default')
params = {'backend': 'pgf',
          'pgf.rcfonts': False,
          'font.family': 'serif',
          'font.size': 14,
          'axes.labelsize': 14,
          'legend.fontsize': 12,
          'xtick.labelsize': 11,
          'ytick.labelsize': 11,
          'text.usetex': True,
          'lines.linewidth': 3,
          'lines.markersize': 5,
          'lines.markeredgewidth': 1,
          'legend.numpoints': 1, 
          'figure.autolayout': True,
          }
rcParams.update(params)


power_25_parcs = [0.0000, 1.5787, 1.5868, 1.1539, 1.4966, 1.6361, 1.0398,
                  1.4047, 1.0398, 1.6361, 1.4966, 1.1539, 1.5868, 1.5787,
                  0.000000]
power_25 = [0.000000e+00, 1.574501e+00, 1.585469e+00,  1.152534e+00,
            1.509132e+00, 1.673297e+00, 1.047551e+00,  1.444753e+00,
            1.047551e+00, 1.673297e+00, 1.509132e+00,  1.152533e+00,
            1.585468e+00, 1.574501e+00, 0.000000e+00 ]
power_25_parcs_adfs = [0.0000, 1.4137, 1.5696, 1.2269, 1.6621, 1.8458, 1.2186,
                       1.6203, 1.2186, 1.8458, 1.6621, 1.2269, 1.5696, 1.4137,
                       0.0000]
power_25_ref = [0.000000e+00, 5.728300e+01, 6.236700e+01, 4.779900e+01,
             6.490700e+01, 7.249600e+01, 4.692900e+01, 6.369600e+01,
             4.692900e+01, 7.249600e+01, 6.490700e+01, 4.779900e+01,
             6.236700e+01, 5.728300e+01, 0.000000e+00]

power_25_ref = np.array(power_25_ref)
power_25_ref /= sum(power_25_ref) / sum(power_25_parcs)
fig4 = plt.figure()
ax4 = fig4.add_subplot(1, 1, 1)
ax4.plot(range(1, len(power_25)+1),
         power_25, 
         c=color[0],
         linestyle=lines[0],
         label='FEMFFUSION $p=3$')
ax4.plot(range(1, len(power_25_parcs)+1),
         power_25_parcs,
         c=color[1],
         linestyle=lines[1],
         marker='o',
         label='PARCS')
ax4.plot(range(1, len(power_25_parcs_adfs)+1),
         power_25_parcs_adfs,
         c=color[2],
         linestyle=lines[2],
         marker='d',
         label='PARCS ADF')
ax4.plot(range(1, len(power_25_ref)+1),
         power_25_ref,
         c=color[3],
         linestyle=lines[3],
         marker='X',
         label='REF')
ax4.set_ylabel("Power")
ax4.set_xlabel("Radial Position")
ax4.legend()
ax4.grid(True)
fig4.savefig("VVER1000CORTEX_radial_power.pdf", format='pdf')
plt.show()
