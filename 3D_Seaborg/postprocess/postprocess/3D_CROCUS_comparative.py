#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modified June 2018

@author: Antoni Vidal
"""

#from utils import compareDistributions
import numpy as np
from utils import compareEig, readEig, latexRow
from utils import get_n_dofs, get_fe_degree, get_n_cells
from utils import get_power, get_n_refinements
from utils import remove_zeros, compareDistributions
files = []

## SELECT PROBLEM
problem = '3D_CROCUS'
#problem = '3D_CUBE'
#problem = '3D_CUBE_vacuum'
if problem == '3D_CUBE':
    folder = '../test/3D_cube/'
    n_blocks = 2
    ass_per_plane = 4
    n_planes = 20
    files.append(folder + 'cube.out')
    files_ref = folder + 'cube_ref.ref'
    
if problem == '3D_CUBE_vacuum':
    folder = '../test/3D_cube/'
    n_blocks = 2
    ass_per_plane = 4
    n_planes = 20
    files.append(folder + 'cube_vacuum.out')
    files_ref = folder + 'cube_ref_vacuum.ref'
    
    
if problem == '3D_CROCUS':
    folder = '../3D_CROCUS/'
    n_blocks = 2
    ass_per_plane = 322
    n_planes = 54
#    files.append(folder + 'CROCUS_FE2_ref0.out')
#    files.append(folder + 'CROCUS_FE3_ref0.out')
    files.append(folder + 'CROCUS2_FE1_ref0.out')
    files.append(folder + 'CROCUS2_FE2_ref0.out')
    files.append(folder + 'CROCUS2_FE2_ref1.out')
    files.append(folder + 'CROCUS2_FE3_ref0.out')
    
    files_ref = folder + 'CROCUS_PARCS.ref'

eig_ref = round(readEig(files_ref)[0], 5)


##
## GET GENERAL PARAMETERS
n_dofs = []
fe_degree = []
n_refs = []
eig = []
eig_err = []
n_cells = []
avg = []
mx = []
for f in files:
    n_dofs.append(get_n_dofs(f))
    fe_degree.append(get_fe_degree(f))
    n_refs.append(get_n_refinements(f))
    n_cells.append(get_n_cells(f))
    eigen = round(readEig(f)[0], 5)
    eig.append(eigen)
    eig_err.append(int(round(compareEig(eigen, eig_ref))))


##
## POWER DISTRIBUTIONS
power_ref = np.array(get_power(files_ref))
power_ref = power_ref / sum(power_ref) * len(power_ref) # Normalize
power = []
avg = []
tall = []
rms = []
mre = []
for i, f in enumerate(files):
    p = np.array(get_power(f))
    p = remove_zeros(p)     # Remove zeros
    p = p / sum(p) * len(p) # Normalize
    power.append(p)
    [avg_, mx_, rms_, _] = compareDistributions(p, power_ref)
    avg.append(avg_)
    mx.append(mx_)
    
## Print Values
print latexRow(['$r$', '$p$', 'n_cells', 'n_dofs', r'\lambda',
  r'\Delta\lambda','AVG', 'MAX'])
for i, f in enumerate(files):
    print latexRow([n_refs[i], fe_degree[i],
      n_cells[i], n_blocks*n_dofs[i],
      eig[i], eig_err[i], avg[i], mx[i] ])


#avg_planes = []
#max_planes =  []
#pow_planes = []
#pow_planes_ref = []
#for pl in range(n_planes):
#    pref_plane = power_ref[pl*ass_per_plane:(pl+1)*ass_per_plane]
#    pow_planes_ref.append(sum(pref_plane))
#pow_planes_ref = np.array(pow_planes_ref)
#pow_planes_ref /= ass_per_plane
#    
## Error Per Plane
#for i, f in enumerate(files):
#    pow_planes.append([])
#    avg_planes.append([])
#    max_planes.append([])
#    
#    for pl in range(n_planes):
#        p_plane= power[i][pl*ass_per_plane:(pl+1)*ass_per_plane]
#        pref_plane = power_ref[pl*ass_per_plane:(pl+1)*ass_per_plane]
#        [avg_, mx_, _, _] = compareDistributions(p_plane, pref_plane)
#        pow_planes[i].append(sum(p_plane))
#        avg_planes[i].append(avg_)
#        max_planes[i].append(mx_)
#
#        
#    pow_planes[i] = np.array(pow_planes[i])
#    pow_planes[i] /= ass_per_plane
#
#
## ----------------------------------------------------------------------- #
##### PLOT ERROR PER PLANE ####
## ----------------------------------------------------------------------- #
#import matplotlib.pyplot as plt
#from matplotlib import rcParams
#
#plt.close('all')
#
#color = ['#1f77b4', '#d62728', '#2ca02c', '#ff7f0e', '#ff9896',
#         '#c5b0d5', '#8c564b', '#e377c2', '#f7b6d2', '#7f7f7f',
#         '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
#lines = ['-', ':', '--', '-'] 
#
#plt.style.use('default')
#params = {'backend': 'pgf',
#          'pgf.rcfonts': False,
#          'font.family': 'serif',
#          'font.size': 14,
#          'axes.labelsize': 14,
#          'legend.fontsize': 12,
#          'xtick.labelsize': 11,
#          'ytick.labelsize': 11,
#          'text.usetex': True,
#          'lines.linewidth': 3,
#          'lines.markersize': 5,
#          'lines.markeredgewidth': 1,
#          'legend.numpoints': 1, 
#          'figure.autolayout': True,
#          }
#rcParams.update(params)
#
#
#fig1 = plt.figure()
#ax1 = fig1.add_subplot(1, 1, 1)
#for i, f in enumerate(files):
#    ax1.plot(range(1, n_planes+1), avg_planes[i],
#             c=color[i],
#             linestyle=lines[i],
#             marker='o',
#             label='$p=$' + str(fe_degree[i]) + ',  $r=$' + str(n_refs[i]))
#ax1.set_ylabel("AVG Error")
#ax1.set_xlabel("Planes")
#ax1.legend()
#ax1.grid(True)
#fig1.savefig(problem + "_avgerror_per_plane.pdf", format='pdf') 
#
## 
#fig2 = plt.figure()
#ax2 = fig2.add_subplot(1, 1, 1)
#for i, f in enumerate(files):
#    ax2.plot(range(1, n_planes+1), max_planes[i],
#             c=color[i],
#             linestyle=lines[i],
#             marker='o',
#             label='$p=$' + str(fe_degree[i]) + ',  $r=$' + str(n_refs[i]))
#ax2.set_ylabel("Max Error")
#ax2.set_xlabel("Planes")
#ax2.legend()
#ax2.grid(True)
#fig2.savefig(problem + "_maxerror_per_plane.pdf", format='pdf')
#
###############################################################################
#fig3 = plt.figure()
#ax3 = fig3.add_subplot(1, 1, 1)
#ax3.plot(pow_planes[2], range(1, n_planes+1),
#             c=color[2],
#             linestyle=lines[0],
#             label='FEMFFUSION $p=3$')
#ax3.plot(pow_planes_ref, range(1, n_planes+1),
#             c=color[1],
#             linestyle=lines[1],
#             marker='o',
#             label='PARCS')
#ax3.set_xlabel("Power")
#ax3.set_ylabel("Plane")
#ax3.legend()
#ax3.grid(True)
#fig3.savefig(problem + "_axial_power.pdf", format='pdf')
#
###############################################################################
#power_27 = [0.8502, 0.9670, 1.2057, 1.6714, 1.9317, 2.1906, 2.4156, 2.5932,
#            2.7157, 2.7782, 2.7781, 2.7155, 2.5928, 2.4148, 2.1894, 1.9299,
#            1.6688, 1.2028, 0.9638,  0.8463]
#power_27_parcs = [8.754350e-01, 9.738996e-01, 1.207365e+00, 1.680610e+00,
#                  1.936611e+00, 2.194763e+00, 2.419532e+00, 2.597071e+00,
#                  2.719468e+00, 2.781831e+00, 2.781748e+00, 2.719208e+00,
#                  2.596602e+00, 2.418772e+00, 2.193566e+00, 1.934796e+00,
#                  1.678022e+00, 1.204529e+00, 9.707886e-01, 8.714368e-01]
#
#fig4 = plt.figure()
#ax4 = fig4.add_subplot(1, 1, 1)
#ax4.plot(range(1, len(power_27)+1),
#         power_27, 
#         c=color[2],
#         linestyle=lines[0],
#         label='FEMFFUSION $p=3$')
#ax4.plot(range(1, len(power_27_parcs)+1),
#         power_27_parcs,
#         c=color[1],
#         linestyle=lines[1],
#         marker='o',
#         label='PARCS')
#ax4.set_ylabel("Power")
#ax4.set_xlabel("Radial Position")
#ax4.legend()
#ax4.grid(True)
#fig4.savefig(problem + "_radial_power.pdf", format='pdf')
#
#plt.show()
#
#
#
