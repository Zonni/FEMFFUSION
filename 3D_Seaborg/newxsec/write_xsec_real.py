#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 13:56:49 2023

@author: zonni
"""
import scipy.io as sio
import numpy

out_file = 'seaborg_p10_real.xsec'


# mat = sio.loadmat('AutoGenFullCoreReflected_res.mat')
mat = sio.loadmat('AutoGenFullCoreReflectedWithSmearedRodMovement_res.mat')
ng = len(mat['INF_TRANSPXS'][0])//2 # 
np = 6
print('Number of energy groups:', ng)
na = len(mat['INF_TRANSPXS']) 
print('Number of subdomains:', na)



n_planes = 10

def print_vector_xml(name, lis, file):
    """
    """
    print('<', name, '>', sep='', file=file)
    for i in lis:
        print(i, end=' ', file=f)
    print('', file=file)
    print('</', name, '>', sep='', file=file)


with open(out_file, 'w') as f:
    print('<?xml version="1.0" encoding="utf-8"?>', file=f)
    mats= list(range(na-1)) *n_planes
    # Replace 63 for 382
    for i in range(len(mats)):
        # replace hardik with shardul
        if mats[i] == 63:
            mats[i] = 382
 
    
    print_vector_xml('Composition', mats, file=f)

    print('<materials ngroups="', ng, '" nprecursors="', np,'">', file=f)
    
    # For each material 
    for mt in range(na):
        print('<mix id="', mt, '">', sep='', file=f)
        print('<name>', file=f)
        print('Material'+str(mt), file=f)
        print('</name>', file=f)


        print_vector_xml('Chi', mat['INF_CHIT'][mt][::2], file=f)
        #print_vector_xml('ChiP', mat['INF_CHIP'][mt][::2], file=f) Calculated in FEMFFUSION
        print_vector_xml('Nu', mat['INF_NUBAR'][mt][::2], file=f)
        print_vector_xml('NuSigF', mat['INF_NSF'][mt][::2], file=f)
        print_vector_xml('SigF', mat['INF_FISS'][mt][::2], file=f)
        print_vector_xml('SigmaA', mat['INF_ABS'][mt][::2], file=f)
        print_vector_xml('SigmaR', mat['INF_REMXS'][mt][::2], file=f)
        print_vector_xml('SigmaTR', mat['INF_TRANSPXS'][mt][::2], file=f)
        print_vector_xml('SigmaT', mat['INF_TOT'][mt][::2], file=f)
#        print_vector_xml('Chi', chi[:], file=f)
        # print_vector_xml('Velocities', velocities[:], file=f)
        print_vector_xml('Velocities', (1.0/mat['INF_INVV'][mt][::2]), file=f)
        print_vector_xml('Beta', (mat['BETA_EFF'][mt][::2][1:]), file=f)
        print_vector_xml('Lambda', (mat['LAMBDA'][mt][::2][1:]), file=f)
        
        
        scat = mat['INF_SP0'][mt][::2]
        chid = mat['INF_CHID'][mt][::2]
        
        print('<SigmaS>', file=f)
        for g1 in range(ng):
            for g2 in range(ng):
                 print(scat[g2*ng+g1], end=' ', file=f)
            print(';', file=f)
        print('</SigmaS>', file=f)
        
        print('<ChiD>', file=f)
        for p in range(np):
            for g in range(ng):
                  print(chid[g], end=' ', file=f)
            print(';', file=f)
        print('</ChiD>', file=f)
        
        print('</mix>', file=f)
    print('</materials>', file=f)
    
mt = 382
print(mat['INF_TOT'][mt][::2])
print(mat['INF_TRANSPXS'][mt][::2])