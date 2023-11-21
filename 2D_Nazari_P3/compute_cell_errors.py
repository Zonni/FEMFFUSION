# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:33:28 2019

@author: zonni
"""
import sys
sys.path.append('../postprocess')
from utils import compare_distributions, compare_eig
from utils import get_eigenvalues, get_flux
from utils import latex_row, get_cpu_time
import numpy as np
import itertools

files = ['P3_SDP1.out',
         'P3_SDP2.out',
         'P3_SDP3.out',
         'P3_SDP1_Nazari.out',
         'P3_SDP2_Nazari.out',
         'P3_SDP3_Nazari.out',
         'P3_SP1.out',
         'P3_SP3.out',
         'P3_SP5.out',
         'P3_SP7.out',
]

npzfile = np.load('P3.npz')
keff_ref = float(npzfile['keff'])


npzfile = np.load('P3_cells.npz')
flux_g1 = npzfile['flux_g1']
flat_list = [item for sublist in flux_g1 for item in sublist]
flux_ref_g1 = np.array(flat_list)
flux_ref_g1 = sum(get_flux(files[2], 1)) / sum(flux_ref_g1) * flux_ref_g1




# EIGENVALUES
eig_r = keff_ref
eigs = [] 
for f in files:
    eigs.append(round(get_eigenvalues(f)[0], 5))

err_eigs = [] 
for eig in eigs:
    err_eigs.append(int(round(compare_eig(eig, eig_r))))

## POWER DISTRIBUTIONS
means = []
rms = []
times = []
flux_ref = flux_ref_g1



for f in files:
    flux = get_flux(f, 1)

    _mean, _, _rms, _ = compare_distributions(flux, flux_ref)
    means.append(_mean)
    rms.append(_rms)
    times.append(get_cpu_time(f))

    
for i in range(len(files)):
    print(latex_row([files[i], eigs[i], err_eigs[i], 
                    rms[i], times[i]]))
    





