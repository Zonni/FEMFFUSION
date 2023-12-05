#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 13:56:49 2023

@author: zonni
"""
import scipy.io as sio
import numpy

out_file = '../periferal_fast/seaborg_periferal_p20.xsec'


# mat = sio.loadmat('AutoGenFullCoreReflected_res.mat')
mat = sio.loadmat('AutoGenFullCoreReflectedWithPeripheralSmearedRodMovement_res.mat')
ng = len(mat['INF_TRANSPXS'][0])//2 # 
np = 6
print('Number of energy groups:', ng)
na = len(mat['INF_TRANSPXS']) 
print('Number of subdomains:', na)

# Non Corrected
beta= [2.263340e-04, 1.223130e-03, 1.167820e-03, 2.670960e-03, 1.152390e-03, 4.785510e-04];
lambda_p = [1.334590e-02, 3.266710e-02, 1.209380e-01, 3.044400e-01, 8.563840e-01, 2.875980e+00];
betazero= [0.00, 0.00,0.00,0.00,0.00,0.00];
lambda_p_zero = [0.00, 0.00,0.00,0.00,0.00,0.00];
# velocities=[2.23517e+09, 4.98880e+08, 3.84974e+07, 5.12639e+06, 1.67542e+06, 7.26031e+05, 2.81629e+05, 8.81629e+04]; #To correct because they are wrong
#chi=[0.761563, 0.23826, 0.00017742, 0.0, 0.0, 0.0, 0.0, 0.0 ];
#chid=numpy.array([ 1.77957E-01,  8.18626E-01,   3.41793E-03,   0.0E+00,  0.0E+00, 0.0E+00, 0.0E+00, 0.0E+00]);

n_planes = 20
central_rod=66; #63 for central rodd


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
        if mats[i] == central_rod:
            mats[i] = 382
 
    
    print_vector_xml('Composition', mats, file=f)

    print('<materials ngroups="', ng, '" nprecursors="', np,'">', file=f)
    
    # For each material 
    for mt in range(na):
        print('<mix id="', mt, '">', sep='', file=f)
        print('<name>', file=f)
        print('Material'+str(mt), file=f)
        print('</name>', file=f)


        print_vector_xml('Chi', mat['INF_CHIP'][mt][::2], file=f)# CAMBIO A MARIO GENFOAM
        #print_vector_xml('ChiP', mat['INF_CHIP'][mt][::2], file=f) Calculated in FEMFFUSION 
        print_vector_xml('Nu', mat['INF_NUBAR'][mt][::2], file=f)
        print_vector_xml('NuSigF', mat['INF_NSF'][mt][::2], file=f)
        print_vector_xml('SigF', mat['INF_FISS'][mt][::2], file=f)
        print_vector_xml('SigmaA', mat['INF_ABS'][mt][::2], file=f)
        #print_vector_xml('SigmaR', mat['INF_TOT'][mt][::2], file=f)
        print_vector_xml('SigmaTR', mat['INF_TRANSPXS'][mt][::2], file=f)
        print_vector_xml('SigmaT', mat['INF_TOT'][mt][::2], file=f)
#        print_vector_xml('Chi', chi[:], file=f)
        # print_vector_xml('Velocities', velocities[:], file=f)
        print_vector_xml('Velocities', (1.0/mat['INF_INVV'][mt][::2]), file=f)
        if (mat['INF_CHIT'][mt][::2][0]<1e-12):
            print_vector_xml('Beta', betazero[:], file=f)
            print_vector_xml('Lambda', lambda_p_zero[:], file=f)
        else:
            print_vector_xml('Beta', beta[:], file=f)
            print_vector_xml('Lambda', lambda_p[:], file=f)
        
        scat = mat['INF_SP0'][mt][::2]
        chid = mat['INF_CHIP'][mt][::2] # CAMBIO A MARIO GENFOAM
        REM = mat['INF_TOT'][mt][::2].copy()
        for g in range(ng):
            REM[g] -=  scat[g*ng+g]
        #print(REM)
        print_vector_xml('SigmaR', REM, file=f)
        
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