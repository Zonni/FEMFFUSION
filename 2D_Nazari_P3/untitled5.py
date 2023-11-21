#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 17:56:49 2023

@author: zonni
"""
import pyvista
import numpy as np

np.random.seed(12)
point_cloud = np.random.random((5, 3))
point_cloud[:, 2] = 0
point_cloud -= point_cloud.mean(0)
pdata = pyvista.PolyData(point_cloud)
pdata['values'] = np.random.random(5)
plane = pyvista.Plane()
plane.clear_data()
plane = plane.interpolate(pdata, sharpness=3.5)
sample = plane.sample_over_line((-0.5, -0.5, 0), (0.5, 0.5, 0))
pl = pyvista.Plotter()
_ = pl.add_mesh(
       pdata, render_points_as_spheres=True, point_size=50
   )
_ = pl.add_mesh(sample, scalars='values', line_width=10)
_ = pl.add_mesh(plane, scalars='values', style='wireframe')


# pyvista.set_default_active_scalars(mesh2)
# field, scalars = mesh2.active_scalars_info
line = sample.get_array('values')
print(line)