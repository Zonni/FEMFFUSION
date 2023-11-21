#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 18:01:41 2023

@author: zonni
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 17:56:49 2023

@author: zonni
"""
import pyvista
import numpy as np
import pyvista as pv

import sys
sys.path.append('../postprocess')
from utils import parse_vtk_file, parse_vtk_grid
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
from matplotlib import rcParams
plt.close('all')

scalar_name= 'phi_g1_eig_1'
file = 'P3_SDP1_Nazari.out.vtk'

point_a = (0, 4.5, 0)
point_b = (10.0, 4.5, 0)
resolution = 100





print(dots)

# Print phi_g1 line
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)

ax1.plot(dots/10)

ax1.set_xlabel("x (cm)")
ax1.set_ylabel("Scalar Flux (AU)")

ax1.grid(True)
ax1.legend()
