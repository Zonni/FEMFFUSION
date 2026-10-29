#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 16:45:10 2026

@author: amanda
"""

import numpy as np
from collections import Counter

def parse_time_file_debug(filename):
    f = open(filename)
    line = ' '
    step, time, power = [], [], []
    while (line != ''):
        line = f.readline()
        if line[0:7] == 'Time in':
            line = f.readline()
            num = line.split()
            step.append(int(num[0]))
            time.append(float(num[1]))
        elif line[0:4] == 'Powe':
            power.append([])
            line = f.readline()
            for num in line.split():
                power[-1].append(float(num))
    f.close()
    return step, time, power

fname = "3D_Langenbuch_ds_fe3_01.out.time"   # tu FOM
step, time, power = parse_time_file_debug(fname)

lengths = [len(p) for p in power]
print("n bloques:", len(power))
print("longitudes unicas:", sorted(set(lengths)))
print("conteo por longitud:", Counter(lengths))

parse_time_file_debug('3D_Langenbuch_ds_fe3_01.out')