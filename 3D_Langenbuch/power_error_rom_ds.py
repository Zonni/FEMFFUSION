# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal.
"""
import sys
#sys.path.append('../postprocess')
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

labels = ['ROM', 'Distributed']
out_files =[
    '3D_Langenbuch_rom_bars.out',
    '3D_Langenbuch_ds_bars_fe1.out'
    ]



powers = []
time = []

for i in range(len(out_files)):
    print(i)
    power = get_td_power(out_files[i])
    powers.append(np.array(power))
    time.append(get_td_time(out_files[i]))



#%% ===========================================================================

## PlotS
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# ax.plot(time, power, label='Δt = 0.1 s')
for i in range(len(out_files)):
    ax.plot(time[i], powers[i],  label=labels[i])
ax.grid(True)
ax.legend()
ax.set_xlabel('t (s)')
ax.set_ylabel('Relative Power')
fig.savefig('power_rom_ds.pdf', format='pdf')



