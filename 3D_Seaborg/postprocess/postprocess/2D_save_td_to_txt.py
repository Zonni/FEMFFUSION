# -*- coding: utf-8 -*-
"""
@author: Antoni Vidal """

from utils import  parse_file
#from utils import   parse_vtk_file, parse_vtk_grid, parse_file
import numpy as np
import matplotlib.pyplot as plt
#import scipy.io as sio
plt.close('all')
#%% ===========================================================================

# Time Domain Files 
looking_freq = 1.0 # Hz
folder  = '../2D_UOX_FA/time_domain/'

#file_sta = folder + '2D_UOX_diffusionout'  
#file_nos = folder + '2D_UOX_diffusion.outnos'
#out_file = folder + '2D_UOX_diffusion'

#file_sta = folder + '2D_UOX_diffusion_1cycleout'  
#file_nos = folder + '2D_UOX_diffusion_1cycle.out.nos'
#out_file = folder + '2D_UOX_diffusion_1cycle'

#file_sta = folder + '2D_UOX_diffusion_1cycleout'  
#file_nos = folder + '2D_UOX_diffusion.out.nos'
#out_file = folder + '2D_UOX_diffusion_1cycle'


#file_sta = folder + '2D_UOX_diffusion_10cycleout'  
#file_nos = folder + '2D_UOX_diffusion_10cycle.out.nos'
#out_file = folder + '2D_UOX_diffusion_10cycle'

#file_sta = folder + '2D_UOX_sp3out'  
#file_nos = folder + '2D_UOX_sp3.outnos'
#out_file = folder + '2D_UOX_sp3'
#
#file_sta = folder + '2D_UOX_diffusion_sint.out'  
#file_nos = folder + '2D_UOX_diffusion_sint.out.nos'
#out_file = folder + '2D_UOX_diffusion_sint'

#file_sta = folder + '2D_UOX_sp3_sint.out'  
#file_nos = folder + '2D_UOX_sp3_sint.out.nos'
#out_file = folder + '2D_UOX_sp3_sint'


#file_sta = folder + '2D_UOX_sp1.out'  
#file_nos = folder + '2D_UOX_sp1.out.nos'
#out_file = folder + '2D_UOX_sp1'

file_sta = folder + '2D_UOX_diffusion.out'  
file_nos = folder + '2D_UOX_diffusion.out.nos'
out_file = folder + '2D_UOX_diffusion_FE1'



ny = 138
nx = 138
n_cells = nx * ny

cmap="viridis" 

   
#%% ===========================================================================
# GET TD DATA FROM FEMFFUSION AND POSTPROCESS IT 

time_fem = parse_file(file_sta, begin='Time vector', n_max_lines=1)
time_fem = time_fem[:-1]
n_steps = len(time_fem)
steps = range(0, n_steps)

sta_g1_fem = parse_file(file_sta, 'Group 1 flux', n_max_lines=ny)
sta_g2_fem = parse_file(file_sta, 'Group 2 flux', n_max_lines=ny)
sta_g1_fem = np.array(sta_g1_fem)
sta_g2_fem = np.array(sta_g2_fem)

noise_g1_fem = np.zeros([n_steps, n_cells])
noise_g2_fem = np.zeros([n_steps, n_cells])
for st in steps:
    noise_g1_fem[st] = parse_file(file_nos,
                              'Noise of group 1 time step ' + str(st) ,
                              n_max_lines=nx)
    noise_g2_fem[st] = parse_file(file_nos,
                              'Noise of group 2 time step ' + str(st),
                              n_max_lines=nx)

# Transpose and normalize
noise_g1_fem = np.transpose(noise_g1_fem)
noise_g2_fem = np.transpose(noise_g2_fem) 

freq = np.fft.rfftfreq(n_steps, d=time_fem[1])
fft_g1 = np.fft.rfft(noise_g1_fem) * 2.0/ n_steps 
fft_g2 = np.fft.rfft(noise_g2_fem) * 2.0/ n_steps


# We cut at looking_freq Hz
cut_freq = int (looking_freq * n_steps * time_fem[1])
assert(freq[cut_freq] == looking_freq)
noise_g1 = np.zeros([n_cells], dtype='cfloat')
noise_g2 = np.zeros([n_cells], dtype='cfloat')
for node in range(n_cells):
    noise_g1[node] = fft_g1[node][cut_freq]
    noise_g2[node] = fft_g2[node][cut_freq]

amp_g1_fem = np.abs(noise_g1)
amp_g2_fem = np.abs(noise_g2)

# FEMFUSSION PERTURBATION is a SINE not a cosine
# This way we convert it
pha_g1_fem = np.angle(noise_g1, deg=True) + 90
pha_g2_fem = np.angle(noise_g2, deg=True) + 90
   

#%% ===========================================================================
# Reshape and save in the a seriers of TXT files

sta_g1_fem = sta_g1_fem.reshape(ny, nx)
sta_g2_fem = sta_g2_fem.reshape(ny, nx)
amp_g1_fem = amp_g1_fem.reshape(ny, nx)
amp_g2_fem = amp_g2_fem.reshape(ny, nx)
pha_g1_fem = pha_g1_fem.reshape(ny, nx)
pha_g2_fem = pha_g2_fem.reshape(ny, nx)


np.savetxt(out_file + '_STATIC_FLX_G1.txt', sta_g1_fem)
np.savetxt(out_file + '_STATIC_FLX_G2.txt', sta_g2_fem)
np.savetxt(out_file + '_NOISE_AMP_G1.txt', amp_g1_fem)
np.savetxt(out_file + '_NOISE_AMP_G2.txt', amp_g2_fem)
np.savetxt(out_file + '_NOISE_PHASE_G1.txt', pha_g1_fem)
np.savetxt(out_file + '_NOISE_PHASE_G2.txt', pha_g2_fem)

