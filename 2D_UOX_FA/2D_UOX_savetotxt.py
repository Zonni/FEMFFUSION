# -*- coding: utf-8 -*-
"""
@author: Antoni Vidal """

from utils import  parse_file, compare_distributions
#from utils import   parse_vtk_file, parse_vtk_grid, parse_file
import numpy as np
import matplotlib.pyplot as plt
#import scipy.io as sio
plt.close('all')
#%% ===========================================================================
problem = '2D_UOX_FA'
case = '2'

if problem == '2D_UOX_FA':
    folder_cor = "../2D_UOX_FA/CORE_SIMplus/"
    file_cor_sta_g1 = folder_cor + 'CORESIM_STATIC_FLX_G1.txt'
    file_cor_sta_g2 = folder_cor + 'CORESIM_STATIC_FLX_G2.txt'
    file_cor_amp_g1 = folder_cor + 'CORESIM_CASE2_NOISE_AMP_G1.txt'
    file_cor_amp_g2 = folder_cor + 'CORESIM_CASE2_NOISE_AMP_G2.txt'
    file_cor_pha_g1 = folder_cor + 'CORESIM_CASE2_NOISE_PHASE_G1.txt'
    file_cor_pha_g2 = folder_cor + 'CORESIM_CASE2_NOISE_PHASE_G2.txt'
    
    looking_freq = 1.0 # Hz
    folder = "../2D_UOX_FA/exercise_2/"
    file_fem_static = folder + '2D_UOX_FA_diff.out'
    file_fem_noise  = folder + '2D_UOX_FA_diff.out'

#%% ===========================================================================
# STATIC FLUX COMPARISON
static_flux_g1_cor = np.loadtxt(file_cor_sta_g1)
static_flux_g2_cor = np.loadtxt(file_cor_sta_g2)
n_nodes = static_flux_g1_cor.size
nx =  int(np.sqrt(n_nodes))

static_flux_g1_cor = static_flux_g1_cor.reshape(1, n_nodes)[0]
static_flux_g2_cor = static_flux_g2_cor.reshape(1, n_nodes)[0] 

static_flux_g1_fem = parse_file(file_fem_static, 'Group 1 flux', n_max_lines=nx)
static_flux_g2_fem = parse_file(file_fem_static, 'Group 2 flux', n_max_lines=nx)

static_flux_g1_fem = np.array(static_flux_g1_fem)
static_flux_g2_fem = np.array(static_flux_g2_fem)



#%% ===========================================================================
# GET DATA FROM CORESIM+

amp_g1_cor = np.loadtxt(file_cor_amp_g1)
amp_g2_cor = np.loadtxt(file_cor_amp_g1)

pha_g1_cor = np.loadtxt(file_cor_pha_g1)
pha_g2_cor = np.loadtxt(file_cor_pha_g2)

amp_g1_cor = amp_g1_cor.reshape(1, n_nodes)[0]
amp_g2_cor = amp_g2_cor.reshape(1, n_nodes)[0]
pha_g1_cor = pha_g1_cor.reshape(1, n_nodes)[0]
pha_g2_cor = pha_g2_cor.reshape(1, n_nodes)[0]
       
#%% ===========================================================================
# GET DATA FROM FEMFFUSION AND POSTPROCESS IT 

time_fem = parse_file(file_fem_static, begin='Time vector', n_max_lines=1)
time_fem = time_fem[:-1]
n_steps = len(time_fem)
steps = range(0, n_steps)

noise_g1_fem = np.zeros([n_steps, n_nodes])
noise_g2_fem = np.zeros([n_steps, n_nodes])
for st in steps:
    noise_g1_fem[st] = parse_file(file_fem_noise,
                              'Noise of group 1 time step ' + str(st) ,
                              n_max_lines=nx)
    noise_g2_fem[st] = parse_file(file_fem_noise,
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
noise_g1 = np.zeros([n_nodes], dtype='cfloat')
noise_g2 = np.zeros([n_nodes], dtype='cfloat')
for node in range(n_nodes):
    noise_g1[node] = fft_g1[node][cut_freq]
    noise_g2[node] = fft_g2[node][cut_freq]

amp_g1_fem = np.abs(noise_g1)
amp_g2_fem = np.abs(noise_g2)
# FEMFUSSION PERTURBATION is a SINE not a cosine
# This way we convert it
pha_g1_fem = np.angle(noise_g1, deg=True) + 90
pha_g2_fem = np.angle(noise_g2, deg=True) + 90


#%% ===========================================================================

np.savetxt(folder + 'FEMFFUSION_STATIC_FLX_G1.txt', static_flux_g1_fem.reshape(nx, nx))
np.savetxt(folder + 'FEMFFUSION_STATIC_FLX_G2.txt', static_flux_g2_fem.reshape(nx, nx))
np.savetxt(folder + 'FEMFFUSION_CASE2_NOISE_PHASE_G1.txt', pha_g1_fem.reshape(nx, nx))
np.savetxt(folder + 'FEMFFUSION_CASE2_NOISE_PHASE_G2.txt', pha_g2_fem.reshape(nx, nx))
np.savetxt(folder + 'FEMFFUSION_CASE2_NOISE_AMP_G2.txt', amp_g1_fem.reshape(nx, nx))
np.savetxt(folder + 'FEMFFUSION_CASE2_NOISE_AMP_G1.txt', amp_g2_fem.reshape(nx, nx))

#%% ===========================================================================
# Normalize 

# Normalize STATIC FLUX
norm  = sum(static_flux_g1_cor)  / sum(static_flux_g1_fem)
static_flux_g1_fem *= norm
static_flux_g2_fem *= norm 

mean_error1, _, _, _ = compare_distributions(
        static_flux_g1_fem,
        static_flux_g1_cor)
mean_error2, _, _, _ = compare_distributions(
        static_flux_g2_fem,
        static_flux_g2_cor)
print('STATIC FLX G1: ', mean_error1, '%')
print('STATIC FLX G2: ', mean_error2, '%')


# Normalize Noise Amplitude
amp_g1_fem /= amp_g1_fem[0]
amp_g2_fem /= amp_g2_fem[0]


#%% COMPUTE DIFFERENCES
mean_error1, _, _, _ = compare_distributions(
        amp_g1_fem,
        amp_g1_cor)
mean_error2, _, _, _ = compare_distributions(
        amp_g2_fem,
        amp_g2_cor )
print('MAGNITUDE FLX G1: ', mean_error1, '%')
print('MAGNITUDE FLX G2: ', mean_error2, '%')


mean_error1, _, _, _ = compare_distributions(
        pha_g1_fem,
        pha_g1_cor)
mean_error2, _, _, _ = compare_distributions(
        pha_g2_fem,
        pha_g2_cor )
print('PHASE FLX G1: ', mean_error1, '%')
print('PHASE FLX G2: ', mean_error2, '%')



import matplotlib.pyplot as plt

plt.figure()
plt.imshow(abs(amp_g2_fem.reshape(nx, nx) - amp_g2_cor.reshape(nx, nx))/ amp_g2_cor.reshape(nx, nx) *100);
plt.title('Thermal Flux Noise Difference')
cbar = plt.colorbar()
cbar.set_label('%', rotation=0)
plt.xlabel('X index')
plt.ylabel('Y index')

plt.figure()
plt.imshow(abs(amp_g1_fem.reshape(nx, nx) - amp_g1_cor.reshape(nx, nx))/ amp_g1_cor.reshape(nx, nx) *100);
cbar = plt.colorbar()
plt.title('Fast Flux Noise Difference')
plt.xlabel('X index')
plt.ylabel('Y index')
cbar.set_label('%', rotation=0)
plt.show()

