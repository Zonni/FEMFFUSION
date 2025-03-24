# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal.
"""
import sys
sys.path.append('../postprocess')
from utils import get_td_power, get_td_time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.interpolate import CubicSpline
from utils import latex_row, parse_time_file, compare_distributions
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

markers = ['-', 'o', '*', '^', 'v',  'X', '+', "d", 'x']
labels = ['FOM','POD-STA-10','POD-STA-20', 'POD-MATS-10', 'LUPOD-MATS-10', 'LUPOD-MATS-60']
out_files =['1D_ROM_td.out', 
            '1D_ROM_pod_sta10.out', 
            '1D_ROM_pod_sta20.out',
            '1D_ROM_rampmat12_10.out',
            '1D_ROM_rampmat12_LUPOD10.out',
            '1D_ROM_rampmat12_LUPOD20.out',
            ]
# out_time_files =['1D_ROM_ds_time.out','1D_ROM_pod_mats_time.out'
#                  ]

#labels = ['FOM','POD-5','RPOD-5']
#out_files =['1D_ROM_ds.out'
#   ]
#out_time_files =['1D_ROM_ds_time.out'
#    ]


powers = []
time = []
# local_pows = []

for i in range(len(out_files)):
    print(i)
    power = get_td_power(out_files[i])
    powers.append(np.array(power))
    time.append(get_td_time(out_files[i]))
    


    # _, _, local_pow = parse_time_file(out_time_files[i])
    # local_pows.append(local_pow)


#%% ===========================================================================


## Plots
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
error_power=[];
# ax.plot(time, power, label='Δt = 0.1 s')
for i in range(len(out_files)):
    ax.plot(time[i], powers[i], markers[i],  label=labels[i], )
    

ax.grid(True)
ax.legend()
ax.set_xlabel('t (s)')
ax.set_ylabel('Relative Power')
fig.savefig('1D_ROM.pdf', format='pdf')



