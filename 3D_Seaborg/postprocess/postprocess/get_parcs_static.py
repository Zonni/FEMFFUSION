#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 12:31:52 2019

@author: zonni
"""
from utils import parse_parcs_fluxes_3D, get_parcs_keff
import numpy as np

f = '../3D_CROCUS/crocus.parcsout'
#f = 'tests/test_3D_parcs.test'
flux_fs, flux_th = parse_parcs_fluxes_3D(f)

flux_fs = np.swapaxes(flux_fs, 2, 0)
flux_th = np.swapaxes(flux_th, 2, 0)

keff = get_parcs_keff(f)

import scipy.io as sio
a_dict = {'FLX1': flux_fs, 'FLX2':flux_th, 'keff': keff}

sio.savemat('3D_CROCUS_STATIC.mat', a_dict)
