# -*- coding: utf-8 -*-
"""
@author: Antoni Vidal """

from utils import  parse_file, parse_file_complex, compare_distributions
#from utils import   parse_vtk_file, parse_vtk_grid, parse_file
import numpy as np
import matplotlib.pyplot as plt
#import scipy.io as sio
plt.close('all')
#%% ===========================================================================
problem = '2D_UOX_FA'
if problem == '2D_UOX_FA':
   
    # FEMFFUSION-fd
    folder = "../2D_UOX_FA/exercise_3_fe2/"
    file_fem_out = folder + '2D_UOX_FA.out'
    
    # FEMFFUSION-td
    folder_td = '../2D_UOX_FA/FEMFFUSION-td/ex3/'
    file_td_sta_g1 = folder_td + 'FEMFFUSION_STATIC_FLX_G1.txt'
    file_td_sta_g2 = folder_td + 'FEMFFUSION_STATIC_FLX_G2.txt'
    file_td_amp_g1 = folder_td + 'FEMFFUSION_EX3_NOISE_AMP_G1.txt'
    file_td_amp_g2 = folder_td + 'FEMFFUSION_EX3_NOISE_AMP_G2.txt'
    file_td_pha_g1 = folder_td + 'FEMFFUSION_EX3_NOISE_PHASE_G1.txt'
    file_td_pha_g2 = folder_td + 'FEMFFUSION_EX3_NOISE_PHASE_G2.txt'
    
    nx = 174
    ny = 138
    n_cells = nx * ny

#%% ===========================================================================
# STATIC FLUX COMPARISON

static_flux_g1_fem = parse_file(file_fem_out, 'Group 1 flux', n_max_lines=ny)
static_flux_g2_fem = parse_file(file_fem_out, 'Group 2 flux', n_max_lines=ny)

static_flux_g1_fem = np.array(static_flux_g1_fem)
static_flux_g2_fem = np.array(static_flux_g2_fem)

#%% ===========================================================================
# GET DATA FROM FEMFFUSION-FREQ

noise_g1_fem = parse_file_complex(file_fem_out,
                                  begin='Flux Noise Group 1',
                                  n_max_lines=ny)
noise_g2_fem = parse_file_complex(file_fem_out,
                                  begin='Flux Noise Group 2',
                                  n_max_lines=ny)

assert(n_cells == len(noise_g1_fem))
assert(n_cells == len(noise_g2_fem))


# AMPLITUDE
amp_g1_fem = np.abs(noise_g1_fem)
amp_g2_fem = np.abs(noise_g2_fem)

# PHASE
pha_g1_fem = np.angle(noise_g1_fem, deg=True)
pha_g2_fem = np.angle(noise_g2_fem, deg=True)



## Normalize Noise Amplitude
norm = np.mean(amp_g2_fem)
amp_g1_fem /= norm
amp_g2_fem /= norm


norm  = 1.0 / max(static_flux_g1_fem)
static_flux_g1_fem *= norm
static_flux_g2_fem *= norm 


#%% ===========================================================================
#
np.savetxt(folder + 'FEMFFUSION_EX3_STATIC_FLX_G1.txt', static_flux_g1_fem.reshape(ny, nx))
np.savetxt(folder + 'FEMFFUSION_EX3_STATIC_FLX_G2.txt', static_flux_g2_fem.reshape(ny, nx))
np.savetxt(folder + 'FEMFFUSION_EX3_NOISE_PHASE_G1.txt', pha_g1_fem.reshape(ny, nx))
np.savetxt(folder + 'FEMFFUSION_EX3_NOISE_PHASE_G2.txt', pha_g2_fem.reshape(ny, nx))
np.savetxt(folder + 'FEMFFUSION_EX3_NOISE_AMP_G1.txt', amp_g1_fem.reshape(ny, nx))
np.savetxt(folder + 'FEMFFUSION_EX3_NOISE_AMP_G2.txt', amp_g2_fem.reshape(ny, nx))

#%% PLOT FEMFFUSSION  RESULTS

plt.figure()
plt.imshow(amp_g1_fem.reshape(ny, nx), origin='lower');
cbar = plt.colorbar()
plt.title('Fast Flux Noise FD')
plt.xlabel('X index')
plt.ylabel('Y index')
cbar.set_label('Noise', rotation=90)
plt.savefig(folder + 'Fast Flux Noise.pdf')

plt.figure()
plt.imshow(amp_g2_fem.reshape(ny, nx), origin='lower');
plt.title('Thermal Flux Noise FD')
cbar = plt.colorbar()
cbar.set_label('Noise', rotation=90)
plt.xlabel('X index')
plt.ylabel('Y index')
plt.savefig(folder + 'Thermal Flux Noise.pdf')


#%% ===========================================================================
# GET AND COMPARE FEMFFUSION-td DATA 

static_flux_g1_td = np.loadtxt(file_td_sta_g1)
static_flux_g2_td = np.loadtxt(file_td_sta_g2)
amp_g1_td = np.loadtxt(file_td_amp_g1)
amp_g2_td = np.loadtxt(file_td_amp_g2)
pha_g1_td = np.loadtxt(file_td_pha_g1)
pha_g2_td = np.loadtxt(file_td_pha_g2)

# Reshape 
static_flux_g1_td = static_flux_g1_td.reshape(1, n_cells)[0]
static_flux_g2_td = static_flux_g2_td.reshape(1, n_cells)[0]
amp_g1_td = amp_g1_td.reshape(1, n_cells)[0]
amp_g2_td = amp_g2_td.reshape(1, n_cells)[0]
pha_g1_td = pha_g1_td.reshape(1, n_cells)[0]
pha_g2_td = pha_g2_td.reshape(1, n_cells)[0]

norm = np.mean(amp_g2_td)
amp_g1_td /= norm
amp_g2_td /= norm

# Normalize STATIC FLUX
norm  = 1.0 / max(static_flux_g1_td)
static_flux_g1_td *= norm
static_flux_g2_td *= norm 

mean_error1, _, _, _ = compare_distributions(
                            static_flux_g1_td,
                            static_flux_g1_fem)
mean_error2, _, _, _ = compare_distributions(
                            static_flux_g2_td,
                            static_flux_g2_fem)
print('STATIC FLX G1: ', mean_error1, '%')
print('STATIC FLX G2: ', mean_error2, '%')


## Normalize Noise Amplitude
mean_error1, _, _, _ = compare_distributions(
                                                amp_g1_td,
                                                amp_g1_fem)
mean_error2, _, _, _ = compare_distributions(
                                                amp_g2_td,
                                                amp_g2_fem)
print('MAGNITUDE FLX G1: ', mean_error1, '%')
print('MAGNITUDE FLX G2: ', mean_error2, '%')


mean_error1, _, _, _ = compare_distributions(
                                                pha_g1_td,
                                                pha_g1_fem)
mean_error2, _, _, _ = compare_distributions(
                                                pha_g2_td,
                                                pha_g2_fem )
print('PHASE FLX G1: ', mean_error1, '%')
print('PHASE FLX G2: ', mean_error2, '%')



plt.figure()
plt.imshow((amp_g2_td.reshape(ny, nx) - amp_g2_fem.reshape(ny, nx))/ amp_g2_td.reshape(ny, nx) *100, origin='lower', vmax= 35.0);
plt.title('Thermal Flux Noise Difference')
cbar = plt.colorbar()
cbar.set_label('%', rotation=0)
plt.xlabel('X index')
plt.ylabel('Y index')
plt.savefig(folder_td + 'Thermal Flux Noise Difference.pdf')
#
#plt.figure()
#plt.imshow((amp_g1_td.reshape(ny, nx) - amp_g1_fem.reshape(ny, nx))/ amp_g1_fem.reshape(ny, nx) *100, origin='lower');
#cbar = plt.colorbar()
#plt.title('Fast Flux Noise Difference')
#plt.xlabel('X index')
#plt.ylabel('Y index')
#cbar.set_label('%', rotation=0)
#plt.savefig(folder_td + 'Fast Flux Noise Difference.pdf')
#
#

#
#plt.figure()
#plt.plot(amp_g2_fem.reshape(ny, nx).diagonal() / static_flux_g2_fem.reshape(ny, nx).diagonal(),
#         label='FEMFFUSION-freq');
#plt.plot(amp_g2_td.reshape(ny, nx).diagonal() / static_flux_g2_td.reshape(ny, nx).diagonal(), 
#         label='FEMFFUSION-time', linewidth=2);
#plt.xlabel('Index')
#plt.ylabel('Relative Thermal Neutron Noise Amplitude (%)')
#plt.savefig('Relative_Thermal_Neutron_Noise_Amplitude.pdf')
#plt.legend()
#plt.grid()




plt.figure()
plt.imshow(amp_g1_td.reshape(ny, nx), origin='lower');
cbar = plt.colorbar()
plt.title('Fast Flux Noise TD')
plt.xlabel('X index')
plt.ylabel('Y index')
cbar.set_label('Noise', rotation=90)
plt.savefig(folder + 'Fast Flux Noise.pdf')

plt.figure()
plt.imshow(amp_g2_td.reshape(ny, nx), origin='lower');
plt.title('Thermal Flux Noise TD')
cbar = plt.colorbar()
cbar.set_label('Noise', rotation=90)
plt.xlabel('X index')
plt.ylabel('Y index')
plt.savefig(folder + 'Thermal Flux Noise.pdf')

plt.show()