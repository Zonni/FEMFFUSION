# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:33:28 2019

@author: zonni
"""
import sys
sys.path.append('../postprocess')
from utils import compare_distributions, compare_eig
from utils import get_n_total_dofs
from utils import get_eigenvalues, get_power, comparePowerC5G7, get_cpu_time 

from utils import latex_row
import numpy as np


files = ['c5g7_SDP1.out',
         'c5g7_SDP2.out',
         'c5g7_SDP3.out',
         'c5g7_SDP1_Nazari.out',
         'c5g7_SDP2_Nazari.out',
         # 'c5g7_SDP3_Nazari.out',
         'c5g7_SP1.out',
         'c5g7_SP3.out',
         'c5g7_SP5.out',
         'c5g7_SP7.out',
]

filename_ref = 'c5g7_transport.ref'


# EIGENVALUES
eig_r = round(get_eigenvalues(filename_ref)[0], 5)
eigs = [] 
for f in files:
    eigs.append(round(get_eigenvalues(f)[0], 5))

err_eigs = [] 
for eig in eigs:
    err_eigs.append(int(round(compare_eig(eig, eig_r))))

## POWER DISTRIBUTIONS
means = []
rms = []
power_ref = get_power(filename_ref)
times = []
dofs = []


for f in files:
    power = get_power(f)
    norm = sum(power_ref)/sum(power)
    power = np.array(power)
    power *= norm
    _mean, _, _rms, _ = compare_distributions(power, power_ref)
    
    cpu = get_cpu_time(f)
    
    means.append(_mean)
    rms.append(_rms)
    times.append(cpu)
    
    dofs.append(get_n_total_dofs(f))
    
    
for i in range(len(files)):
    print(latex_row([files[i], dofs[i], eigs[i], err_eigs[i], 
                   means[i], rms[i], times[i]]))
    





