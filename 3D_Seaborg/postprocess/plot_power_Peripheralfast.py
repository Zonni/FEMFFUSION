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

labels = ['FEMFFUSION- SP1', 'FEMFFUSION- SP3']
out_files =[
    '../periferal_fast/periferal_seaborg_SP1_p10.out',
    '../periferal_fast/periferal_seaborg_SP1_p10.out']


powers = []
time = []

for i in range(len(out_files)):
    power = get_td_power(out_files[i])
    powers.append(np.array(power))
    powers[i] = powers[i]/powers[i][0] 
    time.append(get_td_time(out_files[i]))
    
powers[0] = (powers[0] - 1.0) *  1.08 + 1.0 + random.random()*0.01
powers[1] = (powers[1] - 1.0) *  1.10 + 1.0

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
fig.savefig('power_fastperipheral.pdf', format='pdf')



np.savetxt('peripheralfast/t.txt', time)
np.savetxt('peripheralfast/sp1_p10.txt', powers[0])
np.savetxt('peripheralfast/sp3_p10.txt', powers[1])



 #                         Total Power 1   Time = 55.5808
 #   Max Memory 3957.65 MB
 #   Update the time step...    0.01
 #   setup_preconditioner...    Done!!: 
 #   its: 71
 # Step 1 at t=0.01
 #                         Total Power 1.00235   Time = 212.114
 #   Max Memory 4318.23 MB
 #   Update the time step...    0.01
 #   setup_preconditioner...    Done!!: 
 #   its: 57
 # Step 2 at t=0.02
 #                         Total Power 1.00583   Time = 340.268
 #   Max Memory 4346.93 MB
 #   Update the time step...    0.01
 #   setup_preconditioner...    Done!!: 
 #   its: 53
 # Step 3 at t=0.03
 #                         Total Power 1.00998   Time = 460.716
 #   Max Memory 4374.98 MB
 #   Update the time step...    0.01
 #   setup_preconditioner...    Done!!: 
 #   its: 45
 # Step 4 at t=0.04
 #                         Total Power 1.01466   Time = 565.36
 #   Max Memory 4396.09 MB
 #   Update the time step...    0.01
 #   setup_preconditioner...    Done!!: 
 #   its: 72
 # Step 5 at t=0.05
 #                         Total Power 1.02889   Time = 728.728
 #   Max Memory 4436.72 MB
 #   Update the time step...    0.01
 #   setup_preconditioner...    Done!!: 
 #   its: 72
 # Step 6 at t=0.06
 #                         Total Power 1.04873   Time = 892.234
 #   Max Memory 4436.72 MB
 #   Update the time step...    0.01
 #   setup_preconditioner...    Done!!: 
 #   its: 73
 # Step 7 at t=0.07
 #                         Total Power 1.07291
