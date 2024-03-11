# -*- coding: utf-8 -*-
"""
@author: Antoni Vidal
"""

from utils import  parse_file, parse_file_complex, compare_distributions
from utils import plot_hexagonal_assemblies
#from utils import   parse_vtk_file, parse_vtk_grid, parse_file
import numpy as np
import matplotlib.pyplot as plt
#import scipy.io as sio
plt.close('all')
#%% ===========================================================================
problem = '3D_VVER1000_NOISE'
if problem == '3D_VVER1000_NOISE':
   
    # FEMFFUSION-fd
    folder = "../3D_VVER1000_NOISE/"
    file_fd_out = folder + '3D_VVER1000_FE3.out'
#    folder = '../../femffusion_vibration/3D_VVER1000/fe2/'
#    file_fd_out = folder + '3D_VVER1000_FE2.out'
#    file_fd_nos = folder + '3D_VVER1000_FE2.out.nos'
#    file_fd_static = folder + '3D_VVER1000_FE2.out'
#    
    # FEMFFUSION-td
    looking_freq = 1.0 # Hz
    folder_td = '../../femffusion_vibration/3D_VVER1000/fe3/'
    file_td_nos = folder_td + '3D_VVER1000_FE3.out.nos'
    file_td_static = folder_td + '3D_VVER1000_FE3.out'
    
    # 
    nx_rows = [ 9, 10, 11, 12, 13, 14, 15, 16, 17, 
               16, 15, 14, 13, 12, 11, 10, 9]
    ny = 17
    nz = 24
    assert(len(nx_rows) == ny)
    nxy = sum(nx_rows)
    n_lines = ny * nz
    n_cells = nxy * nz
    pitch = 23.6
#%% ===========================================================================
# STATIC FLUX COMPARISON

static_flux_g1_fd = parse_file(file_fd_out, 'Group 1 flux', n_max_lines=n_lines)
static_flux_g2_fd = parse_file(file_fd_out, 'Group 2 flux', n_max_lines=n_lines)

static_flux_g1_fd = np.array(static_flux_g1_fd)
static_flux_g2_fd = np.array(static_flux_g2_fd)

#%% ===========================================================================
# GET DATA FROM FEMFFUSION-FREQ

noise_g1_fd = parse_file_complex(file_fd_out,
                                  begin='Flux Noise Group 1',
                                  n_max_lines=n_lines)
noise_g2_fd = parse_file_complex(file_fd_out,
                                  begin='Flux Noise Group 2',
                                  n_max_lines=n_lines)

assert(n_cells == len(noise_g1_fd))
assert(n_cells == len(noise_g2_fd))

# Relative neutron noise
noise_g1_fd = noise_g1_fd / static_flux_g1_fd
noise_g2_fd = noise_g2_fd / static_flux_g2_fd

# AMPLITUDE
amp_g1_fd = np.abs(noise_g1_fd) 
amp_g2_fd = np.abs(noise_g2_fd)

# PHASE
pha_g1_fd  = np.angle(noise_g1_fd, deg=True) 
pha_g2_fd  = np.angle(noise_g2_fd, deg=True) 

### Normalize Noise Amplitude
#norm = np.mean(amp_g2_fem)
#amp_g1_fem /= norm
#amp_g2_fem /= norm


#norm  = 1.0 / max(static_flux_g1_fem)
#static_flux_g1_fem *= norm
#static_flux_g2_fem *= norm 

#%% ===========================================================================
# AVERAGED PLANE

static_flux_g1_fd_mean = np.zeros(nxy)
static_flux_g2_fd_mean = np.zeros(nxy)
for cell in range(n_cells):
    static_flux_g1_fd_mean[cell%nxy] += static_flux_g1_fd[cell]
    static_flux_g2_fd_mean[cell%nxy] += static_flux_g2_fd[cell]
    
static_flux_g1_fd_mean /= nz
static_flux_g2_fd_mean /= nz


amp_g1_fd_mean = np.zeros(nxy)
amp_g2_fd_mean = np.zeros(nxy)
pha_g1_fd_mean = np.zeros(nxy)
pha_g2_fd_mean = np.zeros(nxy)
for cell in range(n_cells):
    amp_g1_fd_mean[cell % nxy] += amp_g1_fd[cell]
    amp_g2_fd_mean[cell % nxy] += amp_g2_fd[cell]
    pha_g1_fd_mean[cell % nxy] += pha_g1_fd[cell]
    pha_g2_fd_mean[cell % nxy] += pha_g2_fd[cell]
    
amp_g1_fd_mean /= nz
amp_g2_fd_mean /= nz
pha_g1_fd_mean /= nz
pha_g2_fd_mean /= nz

#%% PLOT FEMFFUSSION  RESULTS

## Plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, static_flux_g1_fd_mean, pitch, nx_rows)
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Static Fast Flux')
fig.savefig(folder + problem + "_static_g1.pdf", format='pdf')
#
#
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, static_flux_g2_fd_mean, pitch, nx_rows)
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Static Thermal Flux')
fig.savefig(folder + problem + "_static_g2.pdf", format='pdf')


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, amp_g1_fd_mean * 100, pitch, nx_rows)
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Relative Fast Neutron Noise Amplitude FD (\%)')
fig.savefig(folder + problem + "_amp_g1_fd.pdf", format='pdf')


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, amp_g2_fd_mean * 100 , pitch, nx_rows)
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Relative Thermal Neutron Noise Amplitude FD (\%)')
fig.savefig(folder + problem + "_amp_g2_fd.pdf", format='pdf')


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, pha_g1_fd_mean, pitch, nx_rows)
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Fast Neutron Noise Phase FD (deg)')
fig.savefig(folder + problem + "_pha_g1_fd.pdf", format='pdf')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, pha_g2_fd_mean, pitch, nx_rows)
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Thermal Neutron Noise Phase FD (deg)')
fig.savefig(folder + problem + "_pha_g2_fd.pdf", format='pdf')
    
#%% ===========================================================================
#
## static Flux
#static_flux_g1_fd = parse_file(file_fd_static, 'Group 1 flux', n_max_lines=n_lines)
#static_flux_g2_fd = parse_file(file_fd_static, 'Group 2 flux', n_max_lines=n_lines)
#static_flux_g1_fd = np.array(static_flux_g1_fd)
#static_flux_g2_fd = np.array(static_flux_g2_fd)
#       
## static Flux
#time_fd = parse_file(file_fd_static, begin='Time vector', n_max_lines=1)
#time_fd = time_fd[:-1]
#n_steps = len(time_fd)
#steps = range(0, n_steps)
#
#noise_g1_fd = np.zeros([n_steps, n_cells])
#noise_g2_fd = np.zeros([n_steps, n_cells])
#for st in steps:
#    noise_g1_fd[st] = parse_file(file_fd_nos,
#                                 'Noise of group 1 time step ' + str(st),
#                                 n_max_lines=n_lines)
#    noise_g2_fd[st] = parse_file(file_fd_nos,
#                                 'Noise of group 2 time step ' + str(st),
#                                  n_max_lines=n_lines)
#    assert(n_cells == len(noise_g1_fd[st]))
#    assert(n_cells == len(noise_g2_fd[st]))
#
#
#
## Transpose and normalize
#noise_g1_fd = np.transpose(noise_g1_fd)
#noise_g2_fd = np.transpose(noise_g2_fd) 
#
#freq   = np.fft.rfftfreq(n_steps, d=time_fd[1])
#fft_g1 = np.fft.rfft(noise_g1_fd) * 2.0/ n_steps 
#fft_g2 = np.fft.rfft(noise_g2_fd) * 2.0/ n_steps
#
## We cut at looking_freq Hz
#cut_freq = int (looking_freq * n_steps * time_fd[1])
#assert(freq[cut_freq] == looking_freq)
#noise_g1 = np.zeros([n_cells], dtype='cfloat')
#noise_g2 = np.zeros([n_cells], dtype='cfloat')
#for cell in range(n_cells):
#    noise_g1[cell] = fft_g1[cell][cut_freq] / static_flux_g1_fd[cell]
#    noise_g2[cell] = fft_g2[cell][cut_freq] / static_flux_g2_fd[cell]
#
#amp_g1_fd = np.abs(noise_g1)
#amp_g2_fd = np.abs(noise_g2)
#
#pha_g1_fd = np.angle(noise_g1, deg=True) 
#pha_g2_fd = np.angle(noise_g2, deg=True)


#%% ===========================================================================
# GET FEMFFUSION-td DATA AND POSTPROCESS IT 

# static Flux
static_flux_g1_td = parse_file(file_td_static, 'Group 1 flux', n_max_lines=n_lines)
static_flux_g2_td = parse_file(file_td_static, 'Group 2 flux', n_max_lines=n_lines)
static_flux_g1_td = np.array(static_flux_g1_td)
static_flux_g2_td = np.array(static_flux_g2_td)
       
# static Flux
time_td = parse_file(file_td_static, begin='Time vector', n_max_lines=1)
time_td = time_td[:-1]
n_steps = len(time_td)
steps = range(0, n_steps)

noise_g1_td = np.zeros([n_steps, n_cells])
noise_g2_td = np.zeros([n_steps, n_cells])
for st in steps:
    noise_g1_td[st] = parse_file(file_td_nos,
                                 'Noise of group 1 time step ' + str(st),
                                 n_max_lines=n_lines)
    noise_g2_td[st] = parse_file(file_td_nos,
                                 'Noise of group 2 time step ' + str(st),
                                  n_max_lines=n_lines)
    assert(n_cells == len(noise_g1_td[st]))
    assert(n_cells == len(noise_g2_td[st]))



# Transpose and normalize
noise_g1_td = np.transpose(noise_g1_td)
noise_g2_td = np.transpose(noise_g2_td) 

freq   = np.fft.rfftfreq(n_steps, d=time_td[1])
fft_g1 = np.fft.rfft(noise_g1_td) * 2.0/ n_steps 
fft_g2 = np.fft.rfft(noise_g2_td) * 2.0/ n_steps

# We cut at looking_freq Hz
cut_freq = int (looking_freq * n_steps * time_td[1])
assert(freq[cut_freq] == looking_freq)
noise_g1 = np.zeros([n_cells], dtype='cfloat')
noise_g2 = np.zeros([n_cells], dtype='cfloat')
for cell in range(n_cells):
    noise_g1[cell] = fft_g1[cell][cut_freq] / static_flux_g1_td[cell]
    noise_g2[cell] = fft_g2[cell][cut_freq] / static_flux_g2_td[cell]

amp_g1_td = np.abs(noise_g1)
amp_g2_td = np.abs(noise_g2)

pha_g1_td = np.angle(noise_g1, deg=True) 
pha_g2_td = np.angle(noise_g2, deg=True)



#%% ===========================================================================
# AVERAGED PLANE

static_flux_g1_td_mean = np.zeros(nxy)
static_flux_g2_td_mean = np.zeros(nxy)
for cell in range(n_cells):
    static_flux_g1_td_mean[cell%nxy] += static_flux_g1_td[cell]
    static_flux_g2_td_mean[cell%nxy] += static_flux_g2_td[cell]
    
static_flux_g1_td_mean /= nz
static_flux_g2_td_mean /= nz


amp_g1_td_mean = np.zeros(nxy)
amp_g2_td_mean = np.zeros(nxy)
pha_g1_td_mean = np.zeros(nxy)
pha_g2_td_mean = np.zeros(nxy)
for cell in range(n_cells):
    amp_g1_td_mean[cell % nxy] += amp_g1_td[cell]
    amp_g2_td_mean[cell % nxy] += amp_g2_td[cell]
    pha_g1_td_mean[cell % nxy] += pha_g1_td[cell]
    pha_g2_td_mean[cell % nxy] += pha_g2_td[cell]
    
amp_g1_td_mean /= nz
amp_g2_td_mean /= nz
pha_g1_td_mean /= nz
pha_g2_td_mean /= nz

## Plot
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#plot_hexagonal_assemblies(fig, ax, static_flux_g1_td_mean, pitch, nx_rows)
#ax.set_xlabel('x (cm)')
#ax.set_ylabel('y (cm)')
#ax.set_title('Static Fast Flux')
#fig.savefig(folder_td + problem + "_static_g1.pdf", format='pdf')
#
#
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#plot_hexagonal_assemblies(fig, ax, static_flux_g2_td_mean, pitch, nx_rows)
#ax.set_xlabel('x (cm)')
#ax.set_ylabel('y (cm)')
#ax.set_title('Static Thermal Flux  (\%)')
#fig.savefig(folder_td + problem + "_static_g2.pdf", format='pdf')


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, amp_g1_td_mean * 100, pitch, nx_rows)
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Relative Fast Neutron Noise Amplitude TD (\%)')
fig.savefig(folder_td + problem + "_amp_g1_td.pdf", format='pdf')


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, amp_g2_td_mean * 100 , pitch, nx_rows)
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Relative Thermal Neutron Noise Amplitude TD (\%)')
fig.savefig(folder_td + problem + "_amp_g2_td.pdf", format='pdf')


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, pha_g1_td_mean, pitch, nx_rows)
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Fast Neutron Noise Phase TD (deg)')
fig.savefig(folder_td + problem + "_pha_g1_td.pdf", format='pdf')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_hexagonal_assemblies(fig, ax, pha_g2_td_mean, pitch, nx_rows)
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Thermal Neutron Noise Phase TD (deg)')
fig.savefig(folder_td + problem + "_pha_g2_td.pdf", format='pdf')
#%% ===========================================================================
# COMPARATIVE 


mean_error1, _, _, _ = compare_distributions(
                            static_flux_g1_fd,
                            static_flux_g1_td)
mean_error2, _, _, _ = compare_distributions(
                            static_flux_g2_fd,
                            static_flux_g2_td)
print('STATIC FLX G1: ', mean_error1, '%')
print('STATIC FLX G2: ', mean_error2, '%')


## Normalize Noise Amplitude
mean_error1, _, _, _ = compare_distributions(
                                             amp_g1_fd,
                                             amp_g1_td)
mean_error2, _, _, _ = compare_distributions(
                                             amp_g2_fd,
                                             amp_g2_td)
print('MAGNITUDE FLX G1: ', mean_error1, '%')
print('MAGNITUDE FLX G2: ', mean_error2, '%')


mean_error1, _, _, _ = compare_distributions(
                                                pha_g1_fd,
                                                pha_g1_td)
mean_error2, _, _, _ = compare_distributions(
                                                pha_g2_fd,
                                                pha_g2_td )
print('PHASE FLX G1: ', mean_error1, '%')
print('PHASE FLX G2: ', mean_error2, '%')
#
#
#
#plt.figure()
#plt.imshow((amp_g2_td.reshape(ny, nx) - amp_g2_td.reshape(ny, nx))/ amp_g2_td.reshape(ny, nx) *100, origin='lower', vmax= 35.0);
#plt.title('Thermal Flux Noise Difference')
#cbar = plt.colorbar()
#cbar.set_label('%', rotation=0)
#plt.xlabel('X index')
#plt.ylabel('Y index')
#plt.savefig(folder_td + 'Thermal Flux Noise Difference.pdf')
##
##plt.figure()
##plt.imshow((amp_g1_td.reshape(ny, nx) - amp_g1_td.reshape(ny, nx))/ amp_g1_td.reshape(ny, nx) *100, origin='lower');
##cbar = plt.colorbar()
##plt.title('Fast Flux Noise Difference')
##plt.xlabel('X index')
##plt.ylabel('Y index')
##cbar.set_label('%', rotation=0)
##plt.savefig(folder_td + 'Fast Flux Noise Difference.pdf')
##
##
#
##
##plt.figure()
##plt.plot(amp_g2_td.reshape(ny, nx).diagonal() / static_flux_g2_td.reshape(ny, nx).diagonal(),
##         label='FEMFFUSION-freq');
##plt.plot(amp_g2_td.reshape(ny, nx).diagonal() / static_flux_g2_td.reshape(ny, nx).diagonal(), 
##         label='FEMFFUSION-time', linewidth=2);
##plt.xlabel('Index')
##plt.ylabel('Relative Thermal Neutron Noise Amplitude (%)')
##plt.savefig('Relative_Thermal_Neutron_Noise_Amplitude.pdf')
##plt.legend()
##plt.grid()
#
#
#
#
#plt.figure()
#plt.imshow(amp_g1_td.reshape(ny, nx), origin='lower');
#cbar = plt.colorbar()
#plt.title('Fast Flux Noise TD')
#plt.xlabel('X index')
#plt.ylabel('Y index')
#cbar.set_label('Noise', rotation=90)
#plt.savefig(folder + 'Fast Flux Noise.pdf')
#
#plt.figure()
#plt.imshow(amp_g2_td.reshape(ny, nx), origin='lower');
#plt.title('Thermal Flux Noise TD')
#cbar = plt.colorbar()
#cbar.set_label('Noise', rotation=90)
#plt.xlabel('X index')
#plt.ylabel('Y index')
#plt.savefig(folder + 'Thermal Flux Noise.pdf')
#
#plt.show()