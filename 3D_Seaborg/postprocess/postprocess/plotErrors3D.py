# -*- coding: utf-8 -*-
from utils import parse_file
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
plt.close('all')

#problem = '3D_IAEA'
problem = 'VVER440'
    
if problem == '3D_IAEA':
    filename_1 = '../3D_IAEA/IAEA_SP1.out'
    filename_3 = '../3D_IAEA/IAEA_SP3.out'
    filename_r = '../3D_IAEA/IAEA_SP5.out'
    mesh_size  = (17, 17, 19)
elif problem == 'VVER440':
    filename_1 = '../VVER440/VVER440_SP1.out'
    filename_3 = '../VVER440/VVER440_SP3.out'
    filename_r = '../VVER440/VVER440_SP5.out'
    mesh_size  = (22,25,12)
else:
    print('ERROR! Set a valid problem ')
    assert False
    

####################################################
# PLOT PARAMETERS
####################################################
params = {'backend': 'pdf',
          'font.family': 'serif',
          'font.size': 14,
          'axes.labelsize': 14,
          'legend.fontsize': 14,
          'xtick.labelsize': 14,
          'ytick.labelsize': 14,
          'text.usetex': True,
          'lines.linewidth': 3,
          'lines.markersize': 5,
          'lines.markeredgewidth': 1,
          'legend.numpoints': 1
          }
rcParams.update(params)   

####################################################
####################################################
# Get Power per cell
n_lines = mesh_size[1] * mesh_size[2]
power_1 = parse_file(filename_1, begin='Power', end='Axial')
power_3 = parse_file(filename_3, begin='Power', end='Axial')
power_r = parse_file(filename_r, begin='Power', end='Axial')


# Calculate the relative eror per cell
error1 = len(power_1)*[0.0]
for i in range(len(power_1)):
        if (abs(power_r[i]) > 1e-3):
            error1[i] = 100* abs((power_1[i] - power_r[i])/power_r[i])

error1 = np.array(error1)            
error1 = error1.reshape((mesh_size[2], mesh_size[1], mesh_size[0]))

error3 = len(power_3)*[0.0]
for i in range(len(power_3)):
        if (abs(power_r[i]) > 1e-3):
            error3[i] = 100* abs((power_3[i] - power_r[i])/power_r[i])
            


error3 = np.array(error3)            
error3 = error3.reshape((mesh_size[2], mesh_size[1], mesh_size[0]))

####################################################
#  RADIAL ERROR
####################################################
# Get Average Radial Error
radial_error = mesh_size[1] * [mesh_size[0] * [0.0]]
radial_error = np.array(radial_error)

for i in range(mesh_size[0]):
    for j in range(mesh_size[1]):
        for k in range(mesh_size[2]):
            radial_error[j][i] += error1[k][j][i]/mesh_size[2]

            
# Plot Radial Error
#fig, ax, cbar = plotCell(radial_error)
#plt.axis('off')
#cbar.set_label('Relative Error (\%)', size=16,  labelpad=12)
#plt.savefig(problem + '_radialError.pdf')


####################################################
#  Axial DISTRIBUTION
####################################################

axial_1 = parse_file(filename_1, begin='Axial', n_max_lines=mesh_size[2])
axial_3 = parse_file(filename_3, begin='Axial', n_max_lines=mesh_size[2])
axial_r = parse_file(filename_r, begin='Axial',n_max_lines=mesh_size[2])


color = ['#1f77b4', '#d62728', '#2ca02c', '#ff7f0e', '#ff9896',
         '#c5b0d5', '#8c564b', '#e377c2', '#f7b6d2', '#7f7f7f',
         '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
lines = ['-', '--', ':', '--'] 
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(range(1, mesh_size[2]+1), axial_1, '--', color=color[0], label='SP1')
ax.plot(range(1, mesh_size[2]+1), axial_3, '-.', color=color[1], label='SP3')
ax.plot(range(1, mesh_size[2]+1), axial_r, ':', color=color[2], label='SP5')
ax.tick_params(labelsize=14)
ax.set_xlim([1, mesh_size[2]])
ax.grid()
ax.legend()
ax.set_xlabel('Plane Number', size=14)
ax.set_ylabel('Axial Power', size=14)

plt.savefig(problem + '_axialDist.pdf')

####################################################
#  Axial ERROR
####################################################

axial_error1 = mesh_size[2] * [0.0]
axial_error1 = np.array(axial_error1)
axial_error3 = np.array(axial_error1)
for i in range(mesh_size[0]):
    for j in range(mesh_size[1]):
        for k in range(mesh_size[2]):
            axial_error1[k] += error1[k][j][i]/(mesh_size[0] * mesh_size[1])
            axial_error3[k] += error3[k][j][i]/(mesh_size[0] * mesh_size[1])



fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(range(1, mesh_size[2]+1), axial_error1, '--', color=color[0], linewidth=3, label='SP1')
ax.plot(range(1, mesh_size[2]+1), axial_error3, '-.', color=color[1], linewidth=3, label='SP3')
ax.tick_params(labelsize=14)
ax.set_xlim([1, mesh_size[2]])
ax.grid()
ax.legend(loc='best')
ax.set_xlabel('Plane Number', size=14)
ax.set_ylabel('Relative Error (\%)', size=14)

plt.savefig(problem + '_axialError.pdf')
            
