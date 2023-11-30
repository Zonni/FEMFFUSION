# -*- coding: utf-8 -*-
"""
@author: Antoni Vidal
"""
import sys
sys.path.append('../postprocess')
from utils import get_td_power, get_td_time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
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

out_file = '3D_Seaborg_100.out'


power_100 = get_td_power(out_file)
power_100= np.array(power_100)
power_100 = power_100/power_100[0] 
time_100 = get_td_time(out_file)

out_file = 'seaborg_SP1.out'
power = get_td_power(out_file)
power= np.array(power)
power = power/power[0] 
time = get_td_time(out_file)

#%% ===========================================================================

## PlotS
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(time, power, label='Δt = 0.1 s')
ax.plot(time_100, power_100,  label='Δt = 0.01 s')
ax.grid(True)
ax.legend()
ax.set_xlabel('t (s)')
ax.set_ylabel('Relative Power')


