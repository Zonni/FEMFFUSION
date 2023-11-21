# -*- coding: utf-8 -*-
from utils import parseMatrix, plotCell, relativeError
import matplotlib.pyplot as plt
from matplotlib import rc

plt.close('all')
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


problem = '2D_BIBLIS_ErrorSP1'
#problem = '2D_BIBLIS_ErrorSP3'


if problem == '2D_BIBLIS_ErrorSP1':
    filename_1 = '../2D_biblis/biblis_SP1.out'
    filename_r = '../2D_biblis/biblis_SP5.out'
    n_lines    = 17
elif problem == '2D_BIBLIS_ErrorSP3':
    filename_1 = '../2D_biblis/biblis_SP3.out'
    filename_r = '../2D_biblis/biblis_SP5.out'
    n_lines    = 17
else:
    print('ERROR! Set a valid problem ')
    assert False


# Get Power per cell
power_1 = parseMatrix(filename_1, 'Power', n_lines)
power_r = parseMatrix(filename_r, 'Power', n_lines)

rel_error = relativeError(power_1, power_r)
fig, ax, cbar = plotCell(rel_error)

plt.axis('off')
cbar.set_label('Relative Error (\%)', size=16,  labelpad=12)
plt.savefig(problem + '.pdf')