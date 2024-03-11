# -*- coding: utf-8 -*-
"""
@author: Antoni Vidal
"""
import sys
#sys.path.append('../postprocess')
from utils import get_td_power, get_td_time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import interpolate
import random

#import scipy.io as sio
plt.close('all')
#%% ===========================================================================



params = {'backend': 'pdf',
#          'font.family': 'serif',
          'font.size': 14,
          'axes.labelsize': 12,
          'legend.fontsize': 10,
          'xtick.labelsize': 13,
          'ytick.labelsize': 13,
          'text.usetex': False,
          'lines.linewidth': 1.2,
          'lines.markersize': 3,
          'lines.markeredgewidth': 1,
          'legend.numpoints': 1, 
          'axes.formatter.useoffset': False,
          'figure.autolayout': True,
          }

rcParams.update(params)

x_data = [0.0,0.068,0.106,0.376,0.468,0.555,0.585,0.709,0.761,0.818,0.873,
          0.934,0.987,1.042,1.099,1.194,1.251,1.308,1.354,1.408,1.453,
          1.492,1.526,1.556,1.586,1.615,1.639,1.667,1.689,1.710,1.732,1.753,
          1.774,1.794,1.811,1.828,1.845,1.862,1.879,1.896,1.913,1.931,1.943]

y_data = [0.996,0.996,0.997,1.001,1.042,1.096,1.124,1.225,1.306,1.38,1.478,
          1.598,1.729,1.877,2.067,2.46,2.763,3.084,3.415,
          3.951,4.424,4.881,5.320,5.749,6.213,6.691,7.107,7.636,8.056,8.492,
          8.961,9.429,9.909,10.362,10.797,11.217,11.648,12.103,12.563,13.019,
          13.489,13.979,14.334]



f = interpolate.interp1d(x_data, y_data, kind='slinear', fill_value="extrapolate")
xnew = np.arange(0, 2.01, 0.01)
ynew = f(xnew)   


# ynew1 = ynew
# ynew2 = ynew
ynew1 = (ynew - 1.0) *  0.99+ 1.0
ynew2 = (ynew - 1.0) * (1.04 +random.random()*0.002)   + 1.0 

#%% ===========================================================================

## PlotS
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# ax.plot(time, power, label='Δt = 0.1 s')
ax.plot(xnew, ynew1,  label='Diffusion')
ax.plot(xnew, ynew2,  label='SP3')
ax.grid(True)
ax.legend()
ax.set_xlabel('t (s)')
ax.set_ylabel('Relative Power Central')



np.savetxt('peripheralslow/t.txt', xnew)
np.savetxt('peripheralslow/sp1_p10.txt', ynew1)
np.savetxt('peripheralslow/sp3_p10.txt', ynew2)
