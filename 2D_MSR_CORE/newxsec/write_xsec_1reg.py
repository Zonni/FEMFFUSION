#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 13:56:49 2023

@author: zonni
"""
import scipy.io as sio

out_file = 'msr_core_3reg.xsec'


# mat = sio.loadmat('AutoGenFullCoreReflected_res.mat')
mat = sio.loadmat('regions3.mat')
ng = len(mat['INF_TRANSPXS'][0])//2 # 
print('Number of energy groups:', ng)
na = len(mat['INF_TRANSPXS']) 
print('Number of subdomains:', na)

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
    mat_ids =  list(range(na-1)) + list(range(na-1)) + list(range(na-1))+ [127]
    print_vector_xml('Composition', mat_ids, file=f)

    print('<materials ngroups="', ng, '">', file=f)
    
    # For each material 
    for mt in range(na):
        print('<mix id="', mt, '">', sep='', file=f)
        print('<name>', file=f)
        print('Material'+str(mt), file=f)
        print('</name>', file=f)


        print_vector_xml('Chi', mat['INF_CHIT'][mt][::2], file=f)
        print_vector_xml('Nu', mat['INF_NUBAR'][mt][::2], file=f)
        print_vector_xml('NuSigF', mat['INF_NSF'][mt][::2], file=f)
        print_vector_xml('SigF', mat['INF_FISS'][mt][::2], file=f)
        print_vector_xml('SigmaA', mat['INF_ABS'][mt][::2], file=f)
        print_vector_xml('SigmaR', mat['INF_REMXS'][mt][::2], file=f)
        print_vector_xml('SigmaTR', mat['INF_TRANSPXS'][mt][::2], file=f)
        print_vector_xml('SigmaT', mat['INF_TOT'][mt][::2], file=f)
        
        scat = mat['INF_SP0'][mt][::2]
    
        
        print('<SigmaS>', file=f)
        for g1 in range(ng):
            for g2 in range(ng):
                 print(scat[g2*ng+g1], end=' ', file=f)
            print(';', file=f)
        print('</SigmaS>', file=f)
        
        print('</mix>', file=f)
    print('</materials>', file=f)
    
# mt = 381
# print(mat['INF_TOT'][mt][::2])
# print(mat['INF_TRANSPXS'][mt][::2])