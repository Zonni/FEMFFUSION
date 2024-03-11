#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 17:39:31 2020

@author: zonni
"""

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from utils import parse_vtk_file, parse_vtk_grid
import matplotlib.pyplot as plt

import vtkplotlib as vpl




folder   = '../2D_test_vibration_hex/'
filename = folder + '2D_test.out.vtk'

[x, y, z] = parse_vtk_grid(filename)
stati_g1 = parse_vtk_file(filename, "Static_Flux_g1")
stati_g2 = parse_vtk_file(filename, "Static_Flux_g2")
noise_g1 = parse_vtk_file(filename, "Noise_g1_Magnitude")
noise_g2 = parse_vtk_file(filename, "Noise_g2_Magnitude")
phase_g1 = parse_vtk_file(filename, "Noise_g1_Phase")
phase_g2 = parse_vtk_file(filename, "Noise_g2_Phase")


reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(filename)
reader.Update()
data = reader.GetOutput()
points = data.GetPoints()
npts = points.GetNumberOfPoints()

triangles=  vtk_to_numpy(data.GetCells().GetData())
ntri = triangles.size//4  # number of cells
tri = np.take(triangles,[n for n in range(triangles.size) if n%4 != 0]).reshape(ntri,3)


n_levels = 10
# Print noise_g1
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
#ax1.tricontour(x, y, tri, noise_g1, levels=n_levels, linewidths=0.5, colors='k')
cntr2 = ax1.tricontourf(x, y, tri, noise_g1, levels=n_levels)
ax1.set_aspect('equal')
fig1.colorbar(cntr2, ax=ax1)
ax1.set_ylabel("y (cm)")
ax1.set_xlabel("x (cm)")
ax1.set_title("Relative Fast Noise Magnitude (\%)")