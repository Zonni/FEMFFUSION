# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 13:18:12 2017

@author: zonni
"""
from utils import compare_distributions
from utils import compare_eig, get_eigenvalues, latex_row
from utils import get_n_dofs, get_fe_degree, get_n_cells, get_mesh_size
from utils import get_power, remove_zeros
#from utils import normalize_power
import numpy as np

#problem = 'Langenbuch'
problem = 'roseton'
#problem = '1D_Moving'
#problem = '3D_CROCUS'
#problem = '3D_VVER1000_CORTEX'
problem = '3D_VVER1000_NOISE'

if problem == 'Langenbuch':
    filename_1 = '../3D_Langenbuch/3D_Langenbuch_fe1.out'
    filename_2 = '../3D_Langenbuch/3D_Langenbuch_fe2.out'
    filename_3 = '../3D_Langenbuch/3D_Langenbuch_fe3.out'
    filename_ref = '../3D_Langenbuch/3D_Langenbuch_fe5.out'
    files = [filename_1, filename_2, filename_3]

elif problem == 'roseton':
    filename_1 = '../3D_roseton/roseton_FE1.out'
    filename_2 = '../3D_roseton/roseton_FE2.out'
    filename_3 = '../3D_roseton/roseton_FE3.out'
    filename_ref = '../3D_roseton/roseton_FE5.out'
    files = [filename_1, filename_2, filename_3]

elif problem == '1D_Moving':
    filename_1 = '../1D_MovingBar/1D_FE1.out'
    filename_2 = '../1D_MovingBar/1D_FE2.out'
    filename_3 = '../1D_MovingBar/1D_FE3.out'
    filename_ref = '../1D_MovingBar/1D_FE9.out'
    files = [filename_1, filename_2, filename_3]

elif problem == '3D_VVER1000_CORTEX':
    filename_1 = '../3D_VVER1000_CORTEX/3D_VVER1000_FE1.out'
    filename_2 = '../3D_VVER1000_CORTEX/3D_VVER1000_FE2.out'
    filename_3 = '../3D_VVER1000_CORTEX/3D_VVER1000_FE3.out'
    filename_4 = '../3D_VVER1000_CORTEX/VVER1000.parcs.ref'
    filename_5 = '../3D_VVER1000_CORTEX/VVER1000_adfs.parcs.ref'
    filename_ref = '../3D_VVER1000_CORTEX/VVER1000.out.ref'

   
    files = [filename_1, filename_2, filename_3, filename_4, filename_5]
    
elif problem == '3D_VVER1000_NOISE':
    filename_1 = '../3D_VVER1000_NOISE/3D_VVER1000_FE1.out'
    filename_2 = '../3D_VVER1000_NOISE/3D_VVER1000_FE2.out'
    filename_3 = '../3D_VVER1000_NOISE/3D_VVER1000_FE3.out'
    filename_ref = '../3D_VVER1000_NOISE/3D_VVER1000_FE5.out'
   
    files = [filename_1, filename_2, filename_3]


elif problem == '3D_CROCUS':
    folder = '../3D_CROCUS/'
    files = []
    files.append(folder + 'CROCUS2_FE1_ref0.out')
    files.append(folder + 'CROCUS2_FE2_ref0.out')
    files.append(folder + 'CROCUS2_FE2_ref1.out')
    files.append(folder + 'CROCUS2_FE3_ref0.out')
    
    filename_ref = folder + 'CROCUS_PARCS.ref'



# FE DEGREE
fe_degs = []
for f in files:
    fe_degs.append(get_fe_degree(f))
    
# CELLS
cells = []
for f in files:
    cells.append(get_n_cells(f))

# DOFS
dofs = []
for f in files:
    dofs.append(2*get_n_dofs(f))

# EIGENVALUES
eig_r = round(get_eigenvalues(filename_ref)[0], 5)
eigs = [] 
for f in files:
    eigs.append(round(get_eigenvalues(f)[0], 5))


err_eigs = [] 
for eig in eigs:
    err_eigs.append(int(round(compare_eig(eig, eig_r))))

## POWER DISTRIBUTIONS
means = []
mxs = []
power_ref = get_power(filename_ref)
power_ref = remove_zeros(power_ref)
#power_ref = normalize_power(power_ref)

for f in files:
    power = get_power(f)
    power = remove_zeros(power)
    #power = normalize_power(power)
    mean, mx, _, _ = compare_distributions(power, power_ref)
    means.append(mean)
    mxs.append(mx)
    
for i in range(len(files)):
    print(latex_row([fe_degs[i], dofs[i], cells[i], eigs[i], err_eigs[i], 
                    means[i], mxs[i]]))
    
###############################################################################
####    3D PART PER PARTS                                                  ####
###############################################################################
#    
# Ref
pow_planes_ref = []
[_, _, n_planes] = get_mesh_size(filename_ref)
power_ref = get_power(filename_ref)
ass_per_plane = len(power_ref) // n_planes
assert(len(power_ref) % n_planes == 0)
for pl in range(n_planes):
    pref_plane = power_ref[pl*ass_per_plane:ass_per_plane*(pl+1)]
    pow_planes_ref.append(sum(pref_plane))
pow_planes_ref = np.array(pow_planes_ref)
pow_planes_ref /= sum(pow_planes_ref) / len(pow_planes_ref)


pow_planes = []
for i, f in enumerate(files):
    pow_planes.append([])
    power = get_power(f)

    [_, _, n_planes] = get_mesh_size(f)
    ass_per_plane = len(power) // n_planes
    assert(len(power) % n_planes == 0)
    for pl in range(n_planes):
        p_plane= power[pl*ass_per_plane:(pl+1)*ass_per_plane]
        pow_planes[i].append(sum(p_plane))
    pow_planes[i] = np.array(pow_planes[i])
    pow_planes[i] /= sum(pow_planes[i]) / len(pow_planes[i])
    
    
avg_planes = len(files) * [np.zeros(n_planes)]
for i, f in enumerate(files):
    for pl in range(n_planes):
        if (pow_planes_ref[pl] != 0.0):
            num = abs(pow_planes[i][pl] - pow_planes_ref[pl])
            avg_planes[i][pl] = num

   
###############################################################################
#####         PLOT ERROR PER PLANE                                         ####
###############################################################################
import matplotlib.pyplot as plt
from matplotlib import rcParams

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

fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(pow_planes[1], range(1, n_planes+1),
             c=color[2],
             linestyle=lines[0],
             label='FEMFFUSION $p=2$')
ax1.plot(pow_planes_ref, range(1, n_planes+1),
             c=color[1],
             linestyle=lines[1],
             marker='o',
             label='PARCS')
ax1.set_xlabel("Power")
ax1.set_ylabel("Plane")
ax1.legend()
ax1.grid(True)
fig1.savefig(problem + "_axial_power.pdf", format='pdf')
