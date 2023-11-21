#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modified June 2018

@author: Antoni Vidal
"""

#from utils import compareDistributions
import numpy as np
from utils import compareEig, readEig, latexRow
from utils import get_n_dofs, get_fe_degree, get_n_cells
from utils import get_axial_refinements, get_radial_refinements
from utils import get_power, comparePowerC5G7

files = []

#problem = '2D_SP1'
#problem = '2D_SP3'
problem = '3D_SP1'
#problem = '3D_SP3'

if problem == '2D_SP1':
    folder = '../2D_C5G7/'
    n_blocks = 7
    files.append(folder + 'c5g7_ref0_FE1.out')
    files.append(folder + 'c5g7_ref0_FE2.out')
    files.append(folder + 'c5g7_ref0_FE3.out')
    files.append(folder + 'c5g7_ref1_FE1.out')
    files.append(folder + 'c5g7_ref1_FE2.out')
    files.append(folder + 'c5g7_ref1_FE3.out')
    files.append(folder + 'c5g7_ref2_FE1.out')
    files.append(folder + 'c5g7_ref2_FE2.out')
    files.append(folder + 'c5g7_ref2_FE3.out')
    files_ref = folder + 'c5g7_transport.ref'

if problem == '2D_SP3':
    folder = '../2D_C5G7/'
    n_blocks = 14
    files.append(folder + 'c5g7_SP3_ref0_FE1.out')
    files.append(folder + 'c5g7_SP3_ref0_FE2.out')
    files.append(folder + 'c5g7_SP3_ref0_FE3.out')
    files.append(folder + 'c5g7_SP3_ref1_FE1.out')
    files.append(folder + 'c5g7_SP3_ref1_FE2.out')
    files.append(folder + 'c5g7_SP3_ref2_FE1.out')
    files_ref = folder + 'c5g7_transport.ref'

if problem == '3D_SP1':
    folder = '../3D_C5G7/'
    n_blocks = 7

    files.append(folder + 'c5g7_FE2_rad1_ax0.out')
    files.append(folder + 'c5g7_FE2_rad1_ax1.out')
    files.append(folder + 'c5g7_FE2_rad1_ax2.out')
    files.append(folder + 'c5g7_FE2_rad1_ax3.out')
    files_ref = folder + 'c5g7_transport.ref'


if problem == '3D_SP3':
    folder = '../3D_C5G7/'
    n_blocks = 14
    files.append(folder + 'c5g7_SP3_FE2_rad1_ax0.out')
    files.append(folder + 'c5g7_SP3_FE2_rad1_ax1.out')
    files.append(folder + 'c5g7_SP3_FE2_rad1_ax2.out')
    files.append(folder + 'c5g7_SP3_FE2_rad1_ax3.out')
    files_ref = folder + 'c5g7_transport.ref'

eig_ref = round(readEig(files_ref)[0], 5)
# Get the values
n_dofs = []
fe_degree = []
n_ref_axi = []
n_ref_rad = []
eig = []
eig_err = []
n_cells = []
for f in files:
    n_dofs.append(get_n_dofs(f))
    fe_degree.append(get_fe_degree(f))
    n_ref_rad.append(get_radial_refinements(f))
    n_ref_axi.append(get_axial_refinements(f))
    n_cells.append(get_n_cells(f))
    eigen = round(readEig(f)[0], 5)
    eig.append(eigen)
    eig_err.append(int(round(compareEig(eigen, eig_ref))))


# POWER DISTRIBUTIONS
power_ref = get_power(files_ref)
power = []
avg = []
tall = []
rms = []
mre = []
for i, f in enumerate(files):
    pow_ = np.array(get_power(f))
    pow_norm = pow_/sum(pow_)*sum(power_ref)
    power.append(pow_norm.tolist())
    avg_, tall_, rms_, mre_ = comparePowerC5G7(power[i], power_ref)
    avg.append(avg_)
    tall.append(tall_)
    rms.append(rms_)
    mre.append(mre_)



# Print Values
print(latexRow(['$r_a$', '$r_r$', '$p$', 'n_cells', 'n_dofs', r'\lambda',
                r'\Delta\lambda','AVG', 'MRE']))
for i, f in enumerate(files):
    print(latexRow([n_ref_axi[i], n_ref_rad[i], fe_degree[i],
                    n_cells[i], n_blocks*n_dofs[i],
                    eig[i], eig_err[i], avg[i], mre[i]]))
