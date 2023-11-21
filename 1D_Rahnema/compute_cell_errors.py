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

files = [
     '1D_bwr_conf1_SDP1.out',
     '1D_bwr_conf1_SDP2.out',
     '1D_bwr_conf1_SDP3.out',
     '1D_bwr_conf1_SP1.out',
     '1D_bwr_conf1_SP3.out',
     '1D_bwr_conf1_SP5.out',
     '1D_bwr_conf1_SP7.out'
    ]


npzfile = np.load('2D_Rahnema.npz')
x_ref= npzfile['x']
fluxref_g1 = npzfile['flux_g1']
fluxref_g2 = npzfile['flux_g2']
keff_ref = npzfile['keff']

x=[1.158, 3.231, 3.231, 3.231, 3.231, 1.158,
   1.158, 3.231, 3.231, 3.231, 3.231, 1.158,
   1.158, 3.231, 3.231, 3.231, 3.231, 1.158,
   1.158, 3.231, 3.231, 3.231, 3.231, 1.158,
   1.158, 3.231, 3.231, 3.231, 3.231, 1.158,
   1.158, 3.231, 3.231, 3.231, 3.231, 1.158,
   1.158, 3.231, 3.231, 3.231, 3.231, 1.158]

x = np.cumsum(x) + x_ref[0]
j = 0
fref_g1 = []
fref_g1.append(0.0)
fref_g2 = []
fref_g2.append(0.0)
n =  []
n.append(0)
for i in range(len(x_ref)):
    if x_ref[i] < x[j]:
        fref_g1[j] += fluxref_g1[i]
        fref_g2[j] += fluxref_g2[i]
        n[j] += 1 
    else:
        j += 1
        print(j, x_ref[i])
        fref_g1.append(fluxref_g1[i])
        fref_g2.append(fluxref_g2[i])
        n.append(1)

n = np.array(n)
fref_g1 =  np.array(fref_g1)
flux_ref_g1 = fref_g1 / n
flux_ref_g1 = 3710.840520 / sum(flux_ref_g1) * flux_ref_g1
fref_g2 =  np.array(fref_g2)
flux_ref_g2 = fref_g2 / n
flux_ref_g2 = 1311.656740000000 / sum(flux_ref_g2) * flux_ref_g2


# EIGENVALUES
eig_r = round(float(keff_ref), 5)
eigs = [] 
for f in files:
    eigs.append(round(get_eigenvalues(f)[0], 5))

err_eigs = [] 
for eig in eigs:
    err_eigs.append(int(round(compare_eig(eig, eig_r))))

## POWER DISTRIBUTIONS
means = []
rms1 = []
rms2 = []
times = []
# flux_ref = get_flux(filename_ref, 1)


for f in files:
    fluxg1= get_flux(f, 1)
    fluxg2 = get_flux(f, 2)
    _mean, _, _rms, _ = compare_distributions(fluxg1, flux_ref_g1)
    rms1.append(_rms)
    _mean, _, _rms, _ = compare_distributions(fluxg2, flux_ref_g2)
    rms2.append(_rms)
    times.append(get_cpu_time(f))
    
for i in range(len(files)):
    print(latex_row([files[i], eigs[i], err_eigs[i], 
                    rms1[i],  rms2[i], round(times[i], 2)]))
    





