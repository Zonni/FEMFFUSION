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

labels = ['SP1']
out_files =[
    '../blockage/seaborg_SP1_p10.out'
    ]



powers = []
time = []

for i in range(len(out_files)):
    print(i)
    power = get_td_power(out_files[i])
    powers.append(np.array(power))
    powers[i] = powers[i]/powers[i][0] 
    time.append(get_td_time(out_files[i]))

# powers[0] = (powers[0] - 1.0) *  0.98 + 1.0
# powers[1] = (powers[1] - 1.0) *  1.01 + 1.0
# powers[2] = (powers[2] - 1.0) *  1.030 + 1.0
# powers[3] = (powers[3] - 1.0) *  1.040 + 1.0
# powers[4] = (powers[4] - 1.0) *  1.060 + 1.0

# out_file = 'seaborg_SP1.out'
# power = get_td_power(out_file)
# power= np.array(power)
# power = power/power[0] 
# time = get_td_time(out_file)

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
ax.set_ylabel('Relative Power Central')
fig.savefig('power_Central.pdf', format='pdf')


# np.savetxt('centralrod/t.txt', time[0])
# np.savetxt('centralrod/sp1_p5.txt', powers[0])
# np.savetxt('centralrod/sp1_p10.txt', powers[1])
# np.savetxt('centralrod/sp1_p20.txt', powers[2])
# np.savetxt('centralrod/sp1_p40.txt', powers[3])
# np.savetxt('centralrod/sp3_p10.txt', powers[4])


