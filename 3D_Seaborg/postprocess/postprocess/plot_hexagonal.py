# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 13:23:02 2019

@author: zonni
"""
import matplotlib.pyplot as plt
from utils import plot_hexagonal_assemblies, get_flux, get_power
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
          'text.usetex': True,
          'lines.linewidth': 2,
          'lines.markersize': 5,
          'lines.markeredgewidth': 1,
          'legend.numpoints': 1, 
          'axes.formatter.useoffset': False,
          'figure.autolayout': True,
          }

rcParams.update(params)

problem = '2D_IAEA_noise'

##
if problem == '2D_IAEA_noise':
    file = '../2D_IAEA_noise/2D_IAEA_FE3.out'
    rows = [8, 9, 10, 11, 12, 13, 14, 15, 14, 13, 12, 11, 10, 9, 8]
    pitch = 20.0
if problem == '3D_VVER440_NOISE':
    file = '../3D_VVER440_NOISE/VVER440_EFPD_007.out'
    rows = [ 5, 10, 13, 14, 17, 18, 19, 20, 21, 20, 21, 22, 21,
            22, 21, 20, 21, 20, 19, 18, 17, 14, 13, 10, 5]
    pitch = 14.7
    
    
"-----------------------------------------------------------------------------"
" GET FLUXES "
fl1 = np.array(get_flux(file, 1))
norm = max(fl1)
fl1 = (fl1/norm)
fl2 = np.array(get_flux(file, 2))/norm

pwr = get_power(file)


"-----------------------------------------------------------------------------"
" PLOTS "
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, pwr, pitch, rows)
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Normalized Power')
fig.savefig(problem + "_power.pdf", format='pdf')


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, fl1, pitch, rows)
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Fast Flux')
fig.savefig(problem + "_flux_g1.pdf", format='pdf')


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, fl2, pitch, rows)
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Thermal Flux')
fig.savefig(problem + "_flux_g2.pdf", format='pdf')



