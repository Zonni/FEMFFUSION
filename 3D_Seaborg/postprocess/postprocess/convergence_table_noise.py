# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 13:18:12 2017

@author: zonni
"""

from utils import compareDistributions
from utils import compareEig, readEig, latexRow
from utils import get_n_dofs, get_flux, get_delta_flux
import numpy as np

#problem = 'Langenbuch'
#problem = 'roseton'
#problem = '1D_Moving'
problem = '2D_IAEA_noise'

if problem == 'Langenbuch':
    filename_1 = '../3D_Langenbuch/3D_Langenbuch_fe1.out'
    filename_2 = '../3D_Langenbuch/3D_Langenbuch_fe2.out'
    filename_3 = '../3D_Langenbuch/3D_Langenbuch_fe3.out'
    filename_r = '../3D_Langenbuch/3D_Langenbuch_fe5.out'
    n_lines = 11*10
elif problem == 'roseton':
    filename_1 = '../3D_roseton/roseton_FE1.out'
    filename_2 = '../3D_roseton/roseton_FE2.out'
    filename_3 = '../3D_roseton/roseton_FE3.out'
    filename_r = '../3D_roseton/roseton_FE5.out'
    n_lines = 5*12
elif problem == '1D_Moving':
    filename_1 = '../1D_MovingBar/1D_FE1.out'
    filename_2 = '../1D_MovingBar/1D_FE2.out'
    filename_3 = '../1D_MovingBar/1D_FE3.out'
    filename_r = '../1D_MovingBar/1D_FE9.out'
    n_lines = 1
elif problem == '2D_IAEA_noise':
    file_1 = '../2D_IAEA_noise/2D_IAEA_FE1.out'
    file_2 = '../2D_IAEA_noise/2D_IAEA_FE2.out'
    file_3 = '../2D_IAEA_noise/2D_IAEA_FE3.out'
    file_5 = '../2D_IAEA_noise/2D_IAEA_FE5.out'
#    file_p = '../2D_IAEA_noise/2D_IAEA.parcs'
    file_ref = '../2D_IAEA_noise/2D_IAEA_FE5_ref2.out'
    #file_ref = '../2D_IAEA_noise/2D_IAEA_ref.ref'
    files = [file_1, file_2, file_3, file_5]
else:
    raise Exception('Not valid problem: ' + problem )


"-----------------------------------------------------------------------------"
# DOFS
dofs = []
for f in files:
    dofs.append(get_n_dofs(f))


# EIGENVALUES
eigs = []
for f in files:
    eigs.append(round(readEig(f)[0], 5))
eig_r = round(readEig(file_ref)[0], 5)

# DELTA Eigenvalue
err_eigs = []
for i in range(len(files)):
    err_eigs.append(int(round(compareEig(eigs[i], eig_r))))

"-----------------------------------------------------------------------------"
# FLUX DISTRIBUTIONS

fl1 = []
fl2 = []
for f in files:
    fl = np.array(get_flux(f, 1))
    norm = max(fl)
    fl1.append(fl/norm)
    fl2.append(np.array(get_flux(f, 2))/norm)
    

fl1_ref = np.array(get_flux(file_ref, 1))
norm = max(fl1_ref)
fl1_ref =  fl1_ref / norm
fl2_ref = np.array(get_flux(file_ref, 2))/norm

# ERROR IN POWER
fl1_mean = []
fl2_mean = []
for i in range(len(files)):
    m, _, _, _ = compareDistributions(fl1[i], fl1_ref)
    fl1_mean.append(m)
    m, _, _, _ = compareDistributions(fl2[i], fl2_ref)
    fl2_mean.append(m)
    
"-----------------------------------------------------------------------------"
# NOISE DISTRIBUTIONS - MAGNITUDE

delta_fl1 = []
delta_fl2 = []
for f in files:
    delta_fl1.append(abs(np.array(get_delta_flux(f, 1))))
    delta_fl2.append(abs(np.array(get_delta_flux(f, 2))))
    
delta_fl1_ref = abs(np.array(get_delta_flux(file_ref, 1)))
delta_fl2_ref = abs(np.array(get_delta_flux(file_ref, 2)))


# ERROR IN POWER
delta_fl1_mean = []
delta_fl2_mean = []
for i in range(len(files)):
    m, _, _, _ = compareDistributions(delta_fl1[i], delta_fl1_ref)
    delta_fl1_mean.append(m)
    m, _, _, _ = compareDistributions(delta_fl2[i], delta_fl2_ref)
    delta_fl2_mean.append(m)


"-----------------------------------------------------------------------------"
# LATEX

for i in range(len(files)):
    print(latexRow([i, dofs[i], err_eigs[i],
                    fl1_mean[i], fl2_mean[i],
                    delta_fl1_mean[i], delta_fl2_mean[i]]))
print(eig_r)
