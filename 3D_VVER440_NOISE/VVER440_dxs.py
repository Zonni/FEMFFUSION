#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Antoni Vidal
"""

dxs_file = 'VVER440_EFPD_007.dxs' 
n_materials = 20208
pert_mat = 4791
 
fp = open(dxs_file, "w")
print("dXS", file=fp)
print("#    dSigma_tr1       dSigma_a1        dSigma_f1      dSigma_12", file=fp)
print("#    real  imag       real  imag       real  imag     real  imag", file=fp)
print("#    dSigma_tr2       dSigma_a2        dSigma_f2", file=fp)
print("#    real  imag       real  imag       real  imag", file=fp)


for mat in range(n_materials):
    if (mat == pert_mat):
        print(str(mat+1) + " 0.0000 0.0000   0.1000 0.0000   0.0000 0.0000   0.0000 0.0000",  file=fp)
        print("  0.0000 0.0000  0.1000 0.0000   0.0000 0.0000 ",  file=fp)
    else:
        print(str(mat+1) + "   0.0000 0.0000   0.0000 0.0000   0.0000 0.0000   0.0000 0.0000",  file=fp)
        print("    0.0000 0.0000   0.0000 0.0000   0.0000 0.0000",  file=fp)


fp.close()
    

