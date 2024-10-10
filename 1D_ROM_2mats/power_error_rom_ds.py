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


labels = ['FOM','POD_3','POD_5','POD_10']
out_files =['1D_ROM_ds.out','1D_ROM_pod_mats_3.out', '1D_ROM_pod_mats_5.out', '1D_ROM_pod_mats_10.out'    ]


#out_time_files =['1D_ROM_ds_time.out','1D_ROM_pod_mats_time.out',
#     ]

# labels = ['FOM','POD-5','RPOD-5']
# out_files =['1D_ROM_ds.out'
#     ]
out_time_files =[ ]
    


powers = []
time = []
local_pows = []

for i in range(len(out_files)):
    print(i)
    power = get_td_power(out_files[i])
    powers.append(np.array(power))
    time.append(get_td_time(out_files[i]))
    

    if out_time_files:
        _, _, local_pow = parse_time_file(out_time_files[i])
        local_pows.append(local_pow)


#%% ===========================================================================


## Plots
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
error_power=[];
# ax.plot(time, power, label='Δt = 0.1 s')
for i in range(len(out_files)):
    print(i)
    ax.plot(time[i], powers[i],  label=labels[i])

    x = time[i]
    y = powers[i]
    spl = CubicSpline(x, y)
    ynew=spl(time[0])

    error_power = np.mean(abs(powers[0]-ynew)/powers[0]) *100
    err_local = []
    max_errs = []
    
    if local_pows:
        for t in range(len(local_pows[0])):
            mean_err, maxs, _, _ = compare_distributions(local_pows[i][t], local_pows[0][t] )
            err_local.append(mean_err)
            max_errs.append(maxs)
        mean_local_err = np.mean(err_local)
    
    
        data = [labels[i], error_power, mean_local_err, max(max_errs)]
        print(latex_row(data))
    

ax.grid(True)
ax.legend()
ax.set_xlabel('t (s)')
ax.set_ylabel('Relative Power')
fig.savefig('langenbuch_rrom.pdf', format='pdf')



