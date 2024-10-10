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
from scipy.interpolate import CubicSpline
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

# labels = ['FOM','POD-3','POD-5','POD-10']
# out_files =['3D_Langenbuch_ds_bars_fe3.out',
#     '3D_Langenbuch_rom3_bars_fe3.out',
#         '3D_Langenbuch_rom5_bars_fe3.out',
#         '3D_Langenbuch_rom10_bars_fe3.out'
#     ]

labels = ['FOM','POD-5','RPOD-5']
out_files =['3D_Langenbuch_ds_bars_fe3.out',
        '3D_Langenbuch_rom5_bars_fe3.out',
        '3D_Langenbuch_rrom5_bars_fe3.out'
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
error_power=[];
# ax.plot(time, power, label='Δt = 0.1 s')
for i in range(len(out_files)):
    ax.plot(time[i], powers[i],  label=labels[i])
    if i>0:
        x = time[i]
        y = powers[i]
        spl = CubicSpline(x, y)
        ynew=spl(time[0])
        error_power.append(np.mean(abs(powers[0]-ynew)/powers[0]))
    
ax.grid(True)
ax.legend()
ax.set_xlabel('t (s)')
ax.set_ylabel('Relative Power')
fig.savefig('langenbuch_rrom.pdf', format='pdf')



