# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 13:23:02 2019

@author: zonni
"""
import matplotlib.pyplot as plt
from utils import plot_hexagonal_assemblies, get_flux, get_power
from utils import get_delta_flux
from matplotlib import rcParams
import numpy as np
import matplotlib as mpl
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
          'axes.titlesize': 14,
          'figure.autolayout': True,
          }

rcParams.update(params)
problem = '3D_VVER440_NOISE'

##
if problem == '3D_VVER440_NOISE':
    file = '../3D_VVER440_NOISE/VVER440_EFPD_007.out'
    rows = [ 5, 10, 13, 14, 17, 18, 19, 20, 21, 20, 21, 22, 21,
            22, 21, 20, 21, 20, 19, 18, 17, 14, 13, 10, 5]
    pitch = 14.7
    
    

#%%---------------------------------------------------------------------------
" GET FLUXES "
fl1 = np.array(get_flux(file, 1))
fl2 = np.array(get_flux(file, 2))
pwr = get_power(file)

n_ass = len(fl1)
n_ass_per_plane = sum(rows)
n_planes = int(n_ass/n_ass_per_plane)

assert(len(fl2) == n_ass)
assert(len(pwr) == n_ass)

# Average plane
fl1_meanplane = np.zeros(n_ass_per_plane)
fl2_meanplane = np.zeros(n_ass_per_plane)
pwr_meanplane = np.zeros(n_ass_per_plane)
for a in range(n_ass):
    fl1_meanplane[a%n_ass_per_plane] += fl1[a]
    fl2_meanplane[a%n_ass_per_plane] += fl2[a]
    pwr_meanplane[a%n_ass_per_plane] += pwr[a]

fl1_meanplane /= n_planes
fl2_meanplane /= n_planes
pwr_meanplane /= n_planes


#%%---------------------------------------------------------------------------
#" PLOTS "
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#plot_hexagonal_assemblies(fig, ax, pwr_meanplane, pitch, rows)
#ax.set_xlabel('x (cm)')
#ax.set_ylabel('y (cm)')
#ax.set_title('Neutron Power')
#fig.savefig(problem + "_power.pdf", format='pdf')
#
#
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#plot_hexagonal_assemblies(fig, ax, fl1_meanplane, pitch, rows)
#ax.set_xlabel('x (cm)')
#ax.set_ylabel('y (cm)')
#ax.set_title('Fast Flux')
#fig.savefig(problem + "_flux_g1.pdf", format='pdf')
#
#
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#plot_hexagonal_assemblies(fig, ax, fl2_meanplane, pitch, rows)
#ax.set_xlabel('x (cm)')
#ax.set_ylabel('y (cm)')
#ax.set_title('Thermal Flux')
#fig.savefig(problem + "_flux_g2.pdf", format='pdf')


#%%---------------------------------------------------------------------------
" GET NOISE - AVERAGE PLANE "
dfl1 = np.array(get_delta_flux(file, 1))
dfl2 = np.array(get_delta_flux(file, 2))

assert(len(dfl1) == n_ass)
assert(len(dfl2) == n_ass)

# Average plane
dfl1_meanplane = np.zeros(n_ass_per_plane)
dfl2_meanplane = np.zeros(n_ass_per_plane)
for a in range(n_ass):
    dfl1_meanplane[a%n_ass_per_plane] += abs(dfl1[a])/fl1[a] * 100
    dfl2_meanplane[a%n_ass_per_plane] += abs(dfl2[a])/fl2[a] * 100

dfl1_meanplane /= n_planes
dfl2_meanplane /= n_planes

# Plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, dfl1_meanplane, pitch, rows,
                          norm=mpl.colors.Normalize(vmin=0, vmax=20))
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Relative Fast Noise Magnitude  (\%)')
fig.savefig(problem + "_dflux_g1.pdf", format='pdf')


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, dfl2_meanplane, pitch, rows,
                          norm=mpl.colors.Normalize(vmin=0, vmax=20))
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Relative Thermal Noise Magnitude (\%)')
fig.savefig(problem + "_dflux_g2.pdf", format='pdf')

#%%---------------------------------------------------------------------------
" GET NOISE - AVERAGE AXIAL LINE "

# Average axial line
dfl1_meanline = np.zeros(n_planes, dtype=np.complex)
dfl2_meanline = np.zeros(n_planes, dtype=np.complex)
for a in range(n_ass):
    dfl1_meanline[int(a//n_ass_per_plane)] += dfl1[a] / fl1[a] * 100
    dfl2_meanline[int(a//n_ass_per_plane)] += dfl2[a] / fl2[a] * 100

dfl1_meanline /= n_ass_per_plane
dfl2_meanline /= n_ass_per_plane

z = range(n_planes)
# Plots
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(abs(dfl1_meanline), z)
ax.plot(abs(dfl2_meanline), z)
ax.grid(True)
ax.set_ylabel('z(cm)')
ax.set_xlabel('Relative Noise Magnitude  (\%)')

ax.legend(['Fast Noise' , 'Thermal Noise'])
fig.savefig(problem + "_dflux_avg_line_magnitude.pdf", format='pdf')

z = range(n_planes)
# Plots
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(np.angle(dfl1_meanline, deg=True), z)
ax.plot(np.angle(dfl2_meanline, deg=True), z)
ax.grid(True)
ax.set_ylabel('z(cm)')
ax.set_xlabel('Relative Noise Phase  (\%)')

ax.legend(['Fast Noise' , 'Thermal Noise'], )
fig.savefig(problem + "_dflux_avg_line_phase.pdf", format='pdf')

