# -*- coding: utf-8 -*-
"""
@author: Antoni Vidal
"""
import sys
sys.path.insert(1, '../postprocess/')

from utils import  parse_file, parse_file_complex, compare_distributions
#from utils import   parse_vtk_file, parse_vtk_grid, parse_file
import numpy as np
import matplotlib.pyplot as plt
#import scipy.io as sio
plt.close('all')
#%% ===========================================================================
problem = '2D_UOX_FA'
if problem == '2D_UOX_FA':
   
    # FEMFFUSION-FD
    code = "FEMFUSIONFDSP4"
    folder = "exercise_3/"
    file_fem_out = folder + '2D_UOX_FA_ex3_sp3.out'
    
    # FEMFFUSION-TD
    folder_td = 'FEMFFUSION-td/ex3/'
    file_td_sta_g1 = folder_td + 'FEMFFUSION_CASE3_STATIC_FLX_G2.txt'
    file_td_sta_g2 = folder_td + 'FEMFFUSION_CASE3_STATIC_FLX_G2.txt'
    file_td_amp_g1 = folder_td + 'FEMFFUSION_CASE3_NOISE_AMP_G1.txt'
    file_td_amp_g2 = folder_td + 'FEMFFUSION_CASE3_NOISE_AMP_G2.txt'
    file_td_pha_g1 = folder_td + 'FEMFFUSION_CASE3_NOISE_PHASE_G1.txt'
    file_td_pha_g2 = folder_td + 'FEMFFUSION_CASE3_NOISE_PHASE_G2.txt'
    
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



# ## Normalize Noise Amplitude
# norm = np.mean(amp_g2_fem)
# amp_g1_fem /= norm
# amp_g2_fem /= norm


# norm  = 1.0 / max(static_flux_g1_fem)
# static_flux_g1_fem *= norm
# static_flux_g2_fem *= norm 


#%% ===========================================================================
# #
# np.savetxt(folder + code + '_EX3_STATIC_FLX_G1.txt', static_flux_g1_fem.reshape(ny, nx))
# np.savetxt(folder + code + '_EX3_STATIC_FLX_G2.txt', static_flux_g2_fem.reshape(ny, nx))
# np.savetxt(folder + code + '_EX3_NOISE_PHASE_G1.txt', pha_g1_fem.reshape(ny, nx))
# np.savetxt(folder + code + '_EX3_NOISE_PHASE_G2.txt', pha_g2_fem.reshape(ny, nx))
# np.savetxt(folder + code + '_EX3_NOISE_AMP_G1.txt', amp_g1_fem.reshape(ny, nx))
# np.savetxt(folder + code + '_EX3_NOISE_AMP_G2.txt', amp_g2_fem.reshape(ny, nx))

