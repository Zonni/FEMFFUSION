# -*- coding: utf-8 -*-
# import sys, os
# sys.path.insert(1, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../postprocess')))
import sys
sys.path.insert(1, '../../postprocess/')
from utils import  parse_file, get_mesh_size
#from utils import   parse_vtk_file, parse_vtk_grid, parse_file
import numpy as np
import matplotlib.pyplot as plt
#import scipy.io as sio
plt.close('all')
#%% ===========================================================================
problem = '1D_mechanical_vibration'

looking_freq = 1.0 # Hz
omega   = 'w0'

file_static = '1D_5cells_zeroflux.out'
file_noise  = file_static + '.nos'

#%% ===========================================================================
# MESH
nx, ny, _ = get_mesh_size(file_static)
x = np.ones(nx) * 0.2
assert(len(x)== nx)
x = np.cumsum(x)
x = np.insert(x, 0, 0)    # Insert a 0 in the first position
x = (x[:-1] + x[1:]) / 2  # Compute the average of consecutive elements

n_cells = nx * ny

#%% ===========================================================================
# GET STATIC FLUX
static_flux_g1 = parse_file(file_static, 'Group 1 flux', n_max_lines=ny)
static_flux_g2 = parse_file(file_static, 'Group 2 flux', n_max_lines=ny)

static_flux_g1 = np.array(static_flux_g1)
static_flux_g2 = np.array(static_flux_g2)

       
#%% ===========================================================================
# GET DATA FROM FEMFFUSION AND POSTPROCESS IT 

time = parse_file(file_static, begin='Time vector', n_max_lines=1)
power = parse_file(file_static, begin='Total Power vector', n_max_lines=1)

# Remove the last element to have an even number of points
time = time[:-1]
power = power[:-1]

n_steps = len(time)
steps = range(0, n_steps)

noise_g1 = np.zeros([n_steps, n_cells])
noise_g2 = np.zeros([n_steps, n_cells])
for st in steps:
    noise_g1[st] = parse_file(file_noise,
                              'Noise of group 1 time step ' + str(st) ,
                              n_max_lines=ny)
    noise_g2[st] = parse_file(file_noise,
                              'Noise of group 2 time step ' + str(st),
                              n_max_lines=ny)

# Transpose to use the FFT directly
noise_g1 = noise_g1.transpose()
noise_g2 = noise_g2.transpose() 

# FAST FOURIER TRANSFORM 
freq = np.fft.rfftfreq(n_steps, d=time[1])
fft_g1 = np.fft.rfft(noise_g1) * 2.0/ n_steps 
fft_g2 = np.fft.rfft(noise_g2) * 2.0/ n_steps


# We cut at looking_freq Hz
cut_freq = int (looking_freq * n_steps * time[1])
assert(freq[cut_freq] == looking_freq)

noise_g1 = np.zeros([n_cells], dtype=complex)
noise_g2 = np.zeros([n_cells], dtype=complex)
for node in range(n_cells):
    noise_g1[node] = fft_g1[node][cut_freq]
    noise_g2[node] = fft_g2[node][cut_freq]
    
noise_g1 = noise_g1.reshape(ny, nx)
noise_g2 = noise_g2.reshape(ny, nx)

rea_g1 = np.real(noise_g1)
rea_g2 = np.real(noise_g2)

ima_g1 = np.imag(noise_g1)
ima_g2 = np.imag(noise_g2) 

amp_g1 = np.abs(noise_g1)
amp_g2 = np.abs(noise_g2)

pha_g1 = np.angle(noise_g1, deg=True) 
pha_g2 = np.angle(noise_g2, deg=True) 


#%% PLOT FEMFFUSSION  RESULTS

plt.figure()
plt.plot(time, power)
plt.title('')
plt.grid(True)
plt.xlabel('t (x)')
plt.ylabel('Neutron Power [AU]')
# cbar.set_label('Noise', rotation=90)
# plt.savefig(problem +'_total_power.png')

# # Static Results 
# plt.figure()
# plt.pcolormesh(X, Y, static_flux_g1, shading='auto')
# cbar = plt.colorbar()
# plt.title('Steady State Fast Flux [AU]')
# plt.xlabel('x (cm)')
# plt.ylabel('y (cm)')
# # cbar.set_label('Noise', rotation=90)
# plt.savefig(folder + code + '_' + omega +'_sta_g1.pdf')

# plt.figure()
# plt.pcolormesh(X, Y, static_flux_g2, shading='auto')
# plt.title('Steady State Thermal Flux [AU]')
# cbar = plt.colorbar()
# # cbar.set_label('Noise', rotation=90)
# plt.xlabel('x (cm)')
# plt.ylabel('y (cm)')
# plt.savefig(folder + code + '_' + omega +'_sta_g2.pdf')

# # Amplitude Results 
# plt.figure()
# plt.pcolormesh(X, Y, amp_g1, shading='auto')
# cbar = plt.colorbar()
# # cbar.set_label('Noise', rotation=90)
# plt.title('Fast Neutron Noise Amplitude [AU]')
# plt.xlabel('x (cm)')
# plt.ylabel('y (cm)')
# plt.savefig(folder + code + '_' + omega +'_amp_g2.pdf')

# plt.figure()
# plt.pcolormesh(X, Y, amp_g2, shading='auto')
# plt.title('Thermal Neutron Noise Amplitude [AU]')
# cbar = plt.colorbar()
# # cbar.set_label('Noise', rotation=90)
# plt.xlabel('x (cm)')
# plt.ylabel('y (cm)')
# plt.savefig(folder + code + '_' + omega + '_amp_g2.pdf')

# # Amplitude Results 
# plt.figure()
# plt.pcolormesh(X, Y, amp_g1, shading='auto')
# cbar = plt.colorbar()
# # cbar.set_label('Noise', rotation=90)
# plt.title('Fast Neutron Noise Amplitude [AU]')
# plt.xlabel('x (cm)')
# plt.ylabel('y (cm)')
# plt.savefig(folder + code + '_' + omega +'_amp_g2.pdf')

# plt.figure()
# plt.pcolormesh(X, Y, amp_g2, shading='auto')
# plt.title('Thermal Neutron Noise Amplitude [AU]')
# cbar = plt.colorbar()
# # cbar.set_label('Noise', rotation=90)
# plt.xlabel('x (cm)')
# plt.ylabel('y (cm)')
# plt.savefig(folder + code + '_' + omega + '_amp_g2.pdf')

# # Relative Amplitude Results 
# plt.figure()
# plt.pcolormesh(X, Y, 100*amp_g1/static_flux_g1, shading='auto')
# cbar = plt.colorbar()
# # cbar.set_label('Noise', rotation=90)
# plt.title(r'Fast Neutron Noise Relative Amplitude [%]')
# plt.xlabel('x (cm)')
# plt.ylabel('y (cm)')
# plt.savefig(folder + code +'_' + omega  +'_rel_amp_g2.pdf')

# plt.figure()
# plt.pcolormesh(X, Y,  100*amp_g2/static_flux_g2, shading='auto')
# plt.title(r'Thermal Neutron Noise Relative Amplitude [%]')
# cbar = plt.colorbar()
# # cbar.set_label('Noise', rotation=90)
# plt.xlabel('x (cm)')
# plt.ylabel('y (cm)')
# plt.savefig(folder + code + '_' + omega + '_rel_amp_g2.pdf')

# # Phase Results 
# plt.figure()
# plt.pcolormesh(X, Y, pha_g1, shading='auto')
# cbar = plt.colorbar()
# plt.title('Fast Flux Noise Phase [deg]')
# plt.xlabel('x (cm)')
# plt.ylabel('y (cm)')
# # cbar.set_label('Noise', rotation=90)
# plt.savefig(folder + code +'_pha_g1.pdf')

# plt.figure()
# plt.pcolormesh(X, Y, pha_g2, shading='auto')
# plt.title('Thermal Flux Noise Phase [deg]')
# cbar = plt.colorbar()
# # cbar.set_label('Noise', rotation=90)
# plt.xlabel('x (cm)')
# plt.ylabel('y (cm)')
# plt.savefig(folder + code + '_' + omega + '_pha_g2.pdf')


# # REAL Results 
# plt.figure()
# plt.pcolormesh(X, Y, rea_g1, shading='auto')
# cbar = plt.colorbar()
# plt.title('Fast Flux Real Part')
# plt.xlabel('x (cm)')
# plt.ylabel('y (cm)')
# # cbar.set_label('Noise', rotation=90)
# plt.savefig(folder + code +'_rea_g1.pdf')

# plt.figure()
# plt.pcolormesh(X, Y, rea_g2, shading='auto')
# plt.title('Thermal Flux Real Part')
# cbar = plt.colorbar()
# # cbar.set_label('Noise', rotation=90)
# plt.xlabel('x (cm)')
# plt.ylabel('y (cm)')
# plt.savefig(folder + code + '_' + omega + '_rea_g2.pdf')

# # IMAGINARY Results 
# plt.figure()
# plt.pcolormesh(X, Y, ima_g1, shading='auto')
# cbar = plt.colorbar()
# plt.title('Fast Flux Imaginary Part')
# plt.xlabel('x (cm)')
# plt.ylabel('y (cm)')
# # cbar.set_label('Noise', rotation=90)
# plt.savefig(folder + code +'_ima_g1.pdf')

# plt.figure()
# plt.pcolormesh(X, Y, ima_g2, shading='auto')
# plt.title('Thermal Flux Imaginary Part')
# cbar = plt.colorbar()
# # cbar.set_label('Noise', rotation=90)
# plt.xlabel('x (cm)')
# plt.ylabel('y (cm)')
# plt.savefig(folder + code + '_' + omega + '_ima_g2.pdf')