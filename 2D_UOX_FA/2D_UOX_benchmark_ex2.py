# -*- coding: utf-8 -*-
"""
@author: Antoni Vidal """

#from utils import  parse_file, parse_file_complex, compare_distributions
#from utils import   parse_vtk_file, parse_vtk_grid, parse_file
import numpy as np
import matplotlib.pyplot as plt
#import scipy.io as sio

#import scipy.io as sio
plt.close('all')



    
#%% ===========================================================================
problem = '2D_UOX_FA'
if problem == '2D_UOX_FA':
    
#    folder_cr = "../2D_UOX_FA/CORE_SIMplus/"
    folder_fd = "2D_UOX_FA/FD/"
    folder_sn = "2D_UOX_FA/SN/"
    codes = ['SN', 'FEMFFUSION']
    labels = codes  
    folders = [folder_sn, folder_fd]
    
    

#%% ===========================================================================

    


#%% ===========================================================================
# SN
sta_g1 = []
sta_g2 = []
amp_g1 = [] 
amp_g2 = []
pha_g1 = [] 
pha_g2 = [] 

for i in range(len(codes)):
    file_sta_g1 = folders[i] + codes[i] + '_CASE2_STATIC_FLX_G1.txt'
    file_sta_g2 = folders[i] + codes[i] + '_CASE2_STATIC_FLX_G2.txt'
    file_amp_g1 = folders[i] + codes[i] + '_CASE2_NOISE_AMP_G1.txt'
    file_amp_g2 = folders[i] + codes[i] + '_CASE2_NOISE_AMP_G2.txt'
    file_pha_g1 = folders[i] + codes[i] + '_CASE2_NOISE_PHASE_G1.txt'
    file_pha_g2 = folders[i] + codes[i] + '_CASE2_NOISE_PHASE_G2.txt'
    
    sta_g1.append(np.loadtxt(file_sta_g1))
    sta_g2.append(np.loadtxt(file_sta_g2))
    amp_g1.append(np.loadtxt(file_amp_g1))
    amp_g2.append(np.loadtxt(file_amp_g2))
    pha_g1.append(np.loadtxt(file_pha_g1))
    pha_g2.append(np.loadtxt(file_pha_g2))
    
    # NORM
    amp_g1[i] = amp_g1[i] / sta_g1[i]
    amp_g2[i] = amp_g2[i] / sta_g2[i]
    
    amp_g1[i] /= amp_g1[i][0][0]
    amp_g2[i] /= amp_g2[i][0][0]


#%% ===========================================================================
# FEMFFUSION-FD

plt.figure()
plt.grid(True)
for i in range(len(codes)): 
    plt.plot(np.diag(amp_g1[i]), label=labels[i], linewidth=2)
plt.xlabel('Index')
plt.ylabel('Amp 1')
plt.legend()


plt.figure()
plt.grid(True)
for i in range(len(codes)): 
    plt.plot(np.diag(amp_g2[i]), label=labels[i], linewidth=2)

plt.xlabel('Index')
plt.ylabel('Amp 2')
plt.legend()

plt.figure()
plt.grid(True)
for i in range(len(codes)): 
    plt.plot(np.diag(pha_g1[i]), label=labels[i], linewidth=2)
plt.xlabel('Index')
plt.ylabel('Pha 1')
plt.legend()


plt.figure()
plt.grid(True)
for i in range(len(codes)): 
    plt.plot(np.diag(pha_g2[i]), label=labels[i], linewidth=2)

plt.xlabel('Index')
plt.ylabel('Pha 2')
plt.legend()



#
#plt.figure()
#plt.plot((amp_g2_td.reshape(ny, nx) - amp_g2_fem.reshape(ny, nx))/ amp_g2_td.reshape(ny, nx) *100, origin='lower', vmax= 35.0);
#plt.title('Thermal Flux Noise Difference')
#cbar = plt.colorbar()
#cbar.set_label('%', rotation=0)
#plt.xlabel('X index')
#plt.ylabel('Y index')
#plt.savefig(folder_td + 'Thermal Flux Noise Difference.pdf')
##