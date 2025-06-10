# -*- coding: utf-8 -*-
""" Some python utils to postprocess FEMFFUSION results."""


import numpy as np
from numpy import mean, square
from math import sqrt, pi
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN
import matplotlib as mpl
#import pyvista as pv

numbers = []
for i in range(10):
    numbers.append(str(i))
numbers.append('-')
numbers.append('+')

def remove_zeros(array):
    """ Remove zeros in an array.
    >>> remove_zeros([0.0, 3.0, 0.0, 5.0, 0.0])
    [3.0, 5.0]
    >>> remove_zeros([-1, 3.0, 0, 5])
    [-1, 3.0, 5]
    """
    return [i for i in array if abs(i) > 1e-5]

def rms(error):
    """ Root Mean Square of an error array.
    >>> rms([sqrt(2), 4])
    3.0
    >>> rms([sqrt(2), 0])
    1.0
    """
    return sqrt(mean(square(error)))


def parse_file_same_line(filename, begin=''):
    """ Parse a File and return what is in the same line.
    It must be defined a the begin title
    >>> filename = 'tests/test1.out.test'
    >>> parse_file_same_line(filename, begin='DoFs per Group:')
    [8.0]
    >>> parse_file_same_line(filename, 'Global Refinements:')
    [0.0]
    """
    f = open(filename)
    line = ' '
    out = []
    while (line != ''):
        line = f.readline()
        # Begin Found
        if line[0:len(begin)] == begin:
            for word in line.split():
                if len(word) > 0  and (word[0] in numbers ):
                    out.append(float(word))

    f.close()
    return out


def parse_file(filename, begin='', end='default', n_max_lines=-1):
    """ Parse a File and return what is below a begin title.
    It must be defined a the begin title and the end word or/and the maximum
    number of lines to read.
    >>> filename = 'tests/test1.out.test'
    >>> parse_file(filename, begin='Neutron Power', n_max_lines=1)
    [1.0, 1.0]
    """
    f = open(filename)
    line = ' '
    out = []
    n_lines_readed = 0
    while (line != ''):
        line = f.readline()
        # Begin Found
        if line[0:len(begin)] == begin:
            for line in f.readlines():
                if (line[0:5] != 'Plane'):
                    n_lines_readed += 1
                    # If n_max_lines is reached terminate
                    if n_lines_readed == (n_max_lines + 1):
                        f.close()
                        return out
                    for word in line.split():
                        # If end is found also terminate
                        if word[0:len(end)] == end and end != 'default':
                            f.close()
                            return out
                        # Append double numbers
                        elif len(word) > 3  and (word[0] in numbers ):
                            out.append(float(word))

    f.close()
    return out

def parse_time_file(filename):
    """ Parse a File 
    filename = 'tests/test1.out.test'
    parse_time_file(filename)
    [1.0, 1.0]
    """
    f = open(filename)
    line = ' '
    step = []
    time = []
    power = []
    while (line != ''):
        line = f.readline()
        # Begin Found
        if line[0:7] == 'Time in':
            line = f.readline()
            num =  line.split()
            step.append(int(num[0]))
            time.append(float(num[1]))
        elif line[0:4] == 'Powe':
            power.append([])
            line = f.readline()
            for num in line.split():
                # If end is found also terminate
                power[-1].append(float(num))



    f.close()
    return step, time, power

def parse_file_complex(filename, begin='', end='default', n_max_lines=-1):
    """ Parse a File and return what is below a begin title.
    It must be defined a the begin title and the end word or/and the maximum
    number of lines to read.
    >>> filename = 'tests/test1.out.test'
    >>> parse_file_complex(filename, begin='Flux Noise Group 1', n_max_lines=4)
    [(-1+2j), (-3+4j), (-5+6j), (-7+8j), (-1+2j), (-3+4j), (-5+6j), (-7+8j)]
    """
    f = open(filename)
    line = ' '
    out = []
    n_lines_readed = 0
    while (line != ''):
        line = f.readline()
        # Begin Found
        if line[0:len(begin)] == begin:
            for line in f.readlines():
                if (line[0:5] != 'Plane'):
                    n_lines_readed += 1
                    # If n_max_lines is reached terminate
                    if n_lines_readed == (n_max_lines + 1):
                        f.close()
                        return out
                    for word in line.split():
                        # If end is found also terminate
                        if word[0:len(end)] == end and end != 'default':
                            f.close()
                            return out
                        # Append double numbers
                        elif (word[0] in numbers ):
                            out.append(complex(word))

    f.close()
    return out


def parse_vtk_file(filename, title):
    """ Parse a File and return what is below a begin title.
    It must be defined a the begin title and the end word or/and the maximum
    number of lines to read.
    >>> filename = 'tests/vtk.test'
    >>> parse_vtk_file(filename, 'Neutronic_Power')
    array([1., 1., 1., 1.])
    """
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    array = VN.vtk_to_numpy(data.GetPointData().GetArray(title))
    return array

def parse_vtk_grid(filename):
    """ Return tu x, y, z np.arrays of the grid given in a vtk file
    >>> filename = 'tests/vtk.test'
    >>> parse_vtk_grid(filename)[0]
    array([0., 1., 1., 2.])
    >>> parse_vtk_grid(filename)[1]
    array([0., 0., 0., 0.])
    >>> parse_vtk_grid(filename)[2]
    array([0., 0., 0., 0.])
    """
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()

    x = np.zeros(data.GetNumberOfPoints())
    y = np.zeros(data.GetNumberOfPoints())
    z = np.zeros(data.GetNumberOfPoints())

    for i in range(data.GetNumberOfPoints()):
        x[i],y[i],z[i] = data.GetPoint(i)
    return x, y, z



def P2C(num):
    """Convert decimal point to decimal coma.
    >>> P2C(1.235)
    '1,235'
    >>> P2C('0.230')
    '0,230'
    """
    return str(num).replace('.',',')


def get_eigenvalues(filename):
    """ Return the eigenvalues in an .out file.
    >>> filename = 'tests/test1.out.test'
    >>> get_eigenvalues(filename)
    [1.588324]
    >>> get_eigenvalues('tests/test2.out.test')
    [1.332924]
    """
    eigs_title ='The Eigenvalues are:'
    eigs = parse_file(filename, eigs_title, end='Problem')

    return eigs


def compare_eig(eig, ref_eig):
    """ Return the Eigenvalue error in pcm.
    >>> compare_eig(0.50000, 1.50000)
    100000.0
    >>> compare_eig(1.0, 0.99998)
    2.0
    """
    # Round to avoid floating point arithmethic rounding error
    return round((abs(eig - ref_eig)*1e5), 3)


def latex_row(array):
    """Print an array as a latex table row.
    >>> latex_row(['Table', 2.0, 3.0])
    'Table & 2.0 & 3.0 \\\\\\\\'
    """
    out = ''
    for i in range(len(array)):
        if (i == len(array)-1):
            out += str(array[i]) + " \\\\"
        else:
            out += str(array[i]) + " & "
    return out


def compare_distributions(test, reference, weights='default', n_decimals=2,
                         ignore_zeros=False):
    """ compare_distributions compares 2 cell averageed distributions.
    It returns the mean relative error (%),
    the max absolute error,
    the root mean square relative error (%)
    and the peaking factor relative
    error (%). It admits mean averages
    >>> test = [1.2, 0.8, 1.2, 0.8]
    >>> ref = [1.0, 1.0, 1.0, 1.0]
    >>> compare_distributions(test, ref)
    (20.0, 0.2, 20.0, 20.0)
    >>> compare_distributions(ref, ref)
    (0.0, 0.0, 0.0, 0.0)
    >>> test = [0.0, 0.0, 0.8, 1.2]
    >>> ref = [0.0, 0.0, 1.0, 1.0]
    >>> compare_distributions(test, ref)
    (10.0, 0.2, 14.14, 20.0)

    It  can be set to ignore the reference zeros
    >>> test = [0.0, 0.0, 0.8, 1.2]
    >>> ref = [0.0, 0.0, 1.0, 1.0]
    >>> compare_distributions(test, ref, ignore_zeros=True)
    (20.0, 0.2, 20.0, 20.0)

    Or set a weighting to perform operations. Normaly the weightening
    to the volume of the cells.
    >>> test = [0.0, 0.0, 0.8, 1.1]
    >>> ref = [0.0, 0.0, 1.0, 1.0]
    >>> weights = [1.0, 2.0, 1.0, 2.0]
    >>> compare_distributions(test, ref, weights=weights)
    (6.67, 0.2, 10.0, 20.0)
    """
    # Hard copy the elements
    ref = reference[:]
    tes = test[:]

    assert(len(tes) == len(ref))
    if (weights == 'default'):
        weights=[1.0]*len(ref)
    assert(len(tes) == len(weights))

    # Remove zeros (in reverse order)
    if ignore_zeros is True:
        for i in range(len(ref)-1, -1, -1):
            if abs(ref[i]) < 1e-4:
                del tes[i]
                del weights[i]
                del ref[i]

    # Convert to python np.array
    tes = np.array(tes)
    ref = np.array(ref)
    weights = np.array(weights)
    weight_sum = sum(weights)
    #assert(abs(sum(weights*tes) - sum(weights * ref)) < 1.0)

    # Mean Relative Error
    mean_error = 0.0
    for i in range(len(ref)):
        if (abs(ref[i]) > 1e-4):
            mean_error += weights[i] * abs(ref[i]- tes[i])/(ref[i]*weight_sum)

    # Max error
    max_error_abs =  max(ref- tes)

    # RMS
    RMS = 0.0
    for i in range(len(ref)):
        if (abs(ref[i]) > 1e-4):
            RMS += weights[i]/weight_sum * square((ref[i]- tes[i])/ref[i])
    RMS = sqrt(RMS)

    # Power peaking Factor Error
    max_ref = max(ref)
    arg_max = np.argmax(ref)
    PPF = abs(tes[arg_max]-max_ref)/max_ref

    #   Round to avoid Floating point errors
    mean_error = round(mean_error*100, n_decimals)
    max_error_abs = round(max_error_abs, n_decimals)
    RMS = round(RMS*100, n_decimals)
    PPF = round(PPF * 100, n_decimals)

    return mean_error, max_error_abs, RMS, PPF


def comparePowerC5G7(test, reference):
    """ comparePowerC5G7 compares 2 cell averageed distributions following
    the methodology of the C5G7 benchmark.
    >>> test = [1.2, 0.8, 1.2, 0.8]
    >>> ref = [1.0, 1.0, 1.0, 1.0]
    >>> comparePowerC5G7(test, ref)
    (20.0, 20.0, 20.0, 20.0)
    """
    n_decimals = 2;
    assert(len(test) == len(reference))

    # Hard copy the elements
    ref = reference[:]
    tes = test[:]

    # Remove zeros (in reverse order)
    for i in range(len(ref)-1, -1, -1):
        if abs(ref[i]) < 1e-4:
            del tes[i]
            del ref[i]

    # Convert to python np.array
    tes = np.array(tes)
    ref = np.array(ref)
    weight_sum = len(ref)

    # Mean Relative Error
    mean_error = 0.0
    for i in range(len(ref)):
        mean_error += abs(ref[i]- tes[i])/(ref[i] *weight_sum)

    # Max Relative error
    max_error =  max(abs((ref- tes)/ref))

    # RMS
    RMS = 0.0
    for i in range(len(ref)):
        if (abs(ref[i]) > 1e-4):
            RMS += 1.0/weight_sum * square((ref[i]- tes[i])/ref[i])
    RMS = sqrt(RMS)

    # MRE
    MRE = 0.0
    sum_ref = sum(ref)
    for i in range(len(ref)):
        if (abs(ref[i]) > 1e-4):
            MRE += 1.0/(sum_ref) * abs((ref[i]- tes[i])/ref[i]) * ref[i]

    #   Round to avoid Floating point errors
    mean_error = round(mean_error * 100, n_decimals)
    max_error = round(max_error * 100, n_decimals)
    RMS = round(RMS * 100, n_decimals)
    MRE = round(MRE * 100, n_decimals)

    return mean_error, max_error, RMS, MRE

def homogeneous_reactor(tr1, tr2, a1, a2, s12, f1, f2, L1, L2, n):
    """ Return the fundamentalmental eigenvalue for an 2D homogeneous reactor
    with zero boundary conditions.
    >>> homogeneous_reactor(1.0, 1.0, 0.1, 0.1, 0.1, 0.25, 0.25, 2.0, 2.0, 1)
    0.14327
    >>> homogeneous_reactor(1.0, 1.0, 0.1, 0.1, 0.1, 0.25, 0.25, 2.0, 2.0, 2)
    0.05935
    """
    # Not always valid but for now...
    (range(1, 4),2)
    if (n==1):
        m=1
        l=1
    elif (n==2):
        m=1
        l=2
    elif (n==3):
        m=2
        l=1
    elif (n==4):
        m=2
        l=2
    else:
        print('ERROR n not valid')

    mu = (m*pi/L1)**2 + (l*pi/L2)**2;
    D1 = 1.0/(3*tr1)
    D2 = 1.0/(3*tr2)
    eig = (f1*(D2*mu+a2)+f2*s12)/((D2*mu+a2)*(a1+s12+D1*mu));
    return round(eig, 5)


def get_n_total_dofs(filename):
    """ Get the number of DoFs from a .out file
    >>> filename = 'tests/test1.out.test'
    >>> get_n_total_dofs(filename)
    16
    """
    out = parse_file_same_line(filename, begin='Total DoFs:')
    assert(len(out) == 1)
    return int(out[0])

def get_n_dofs(filename):
    """ Get the number of DoFs from a .out file
    >>> filename = 'tests/test1.out.test'
    >>> get_n_dofs(filename)
    8
    """
    out = parse_file_same_line(filename, begin='DoFs per Group:')
    assert(len(out) == 1)
    return int(out[0])

def get_n_refinements(filename):
    """ Get the number of refinements from a .out file
    >>> filename = 'tests/test1.out.test'
    >>> get_n_refinements(filename)
    0
    """
    out = parse_file_same_line(filename, 'Global Refinements:')
    assert(len(out) == 1)
    return int(out[0])

def get_fe_degree(filename):
    """ Get the fe degree from a .out file
    >>> filename = 'tests/test1.out.test'
    >>> get_fe_degree(filename)
    3
    """
    out = parse_file_same_line(filename, 'Degree of FE:')
    assert(len(out) == 1)
    return int(out[0])

def get_cpu_time(filename):
    """ Get the fe degree from a .out file
    >>> filename = 'tests/test1.out.test'
    >>> get_cpu_time(filename)
    0.000907
    """
    out = parse_file_same_line(filename, 'CPU Time:')
    assert(len(out) == 1)
    return out[0]

def get_n_cells(filename):
    """ Get the number of cells form a .out file
    >>> filename = 'tests/test1.out.test'
    >>> get_n_cells(filename)
    2
    """
    out = parse_file_same_line(filename, 'Number of active cells:')
    assert(len(out) == 1)
    return int(out[0])

def get_mesh_size(filename):
    """ Get the number of cells form a .out file
    >>> filename = 'tests/test1.out.test'
    >>> get_mesh_size(filename)
    [2, 2, 2]
    """
    size = parse_file_same_line(filename, begin='Mesh Size:')
    assert(len(size) == 3)
    size = [int(i) for i in size]
    return size
  
def get_power(filename):
    """ Get the number of cells form a .out file
    >>> filename = 'tests/test1.out.test'
    >>> get_power(filename)
    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    """
    size = get_mesh_size(filename)
    n_lines = size[1] * size[2]
    power = parse_file(filename, begin='Neutron Power', n_max_lines=n_lines)
    return power

def get_flux(filename, group):
    """ Get the flux vectors
    >>> filename = 'tests/test1.out.test'
    >>> get_flux(filename, 1)
    [10.45233, 10.45233, 10.45233, 10.45233, 10.45233, 10.45233, 10.45233, 10.45233]
    >>> get_flux(filename, 2)
    [1.881979, 1.881979, 1.881979, 1.881979, 1.881979, 1.881979, 1.881979, 2.881979]
    """
    size = get_mesh_size(filename)
    n_lines = size[1] * size[2]
    flux = parse_file(filename, begin='Group ' + str(group) +
                     ' flux', n_max_lines=n_lines)
    return flux

def get_delta_flux(filename, group):
    """ Get the delta flux vectors
    >>> filename = 'tests/test1.out.test'
    >>> get_delta_flux(filename, 1)
    [(-1+2j), (-3+4j), (-5+6j), (-7+8j), (-1+2j), (-3+4j), (-5+6j), (-7+8j)]
    >>> get_delta_flux(filename, 2)
    [(-1+2j), (-3+4j), (-5+6j), (-7+8j), (-1+2j), (-3+4j), (-5+6j), (-7+8j)]
    """
    size = get_mesh_size(filename)
    n_lines = size[1] * size[2]
    delta_flux = parse_file_complex(filename, 
                                    begin='Flux Noise Group ' + str(group),
                                    n_max_lines=n_lines)
    return delta_flux



def plot_hexagonal_assemblies(fig, ax, array, pitch, rows, 
                              norm='default'):
    """ Plot an array of values on an hexagonal reactor.
    Each value corresponds to an hexagonal assembly. """
    
    assert(sum(rows)==len(array))
    mx_rows = max(rows)
    n_rows = len(rows)

    ri = pitch / 2 # Radio circumferencia inscrita
    rc = 2/3*np.sqrt(3) * ri # Radio circumferencia circumscrita == lado

    ax.set_aspect('equal')
    xmax = pitch *(mx_rows)
    ymax = ((n_rows//2+1)+0.5*(n_rows//2))  * (2*rc) 
    
    margin = pitch/2
    ax.set_xlim([-margin, xmax+margin])
    ax.set_ylim([-margin, ymax+margin])
    
    # Set Colormap
    cmap = mpl.cm.get_cmap('RdBu_r')
    if (norm=='default'):
        norm=mpl.colors.Normalize(vmin=min(array), vmax=max(array))
    # Draw hexagons
    n = 0
    for j in range(len(rows)):
        for i in range(rows[j]):
            center_x = ((mx_rows - rows[j])/2 + i )  * pitch + pitch/2
            center_y = j * 3/2 * rc + rc
            hexagon = mpl.patches.RegularPolygon((center_x, center_y),
                                                 numVertices=6,
                                                 radius=rc, 
                                                 orientation=0, 
                                                 facecolor=cmap(norm(array[n])), 
                                                 alpha=0.9,
                                                 edgecolor='k')
            ax.add_patch(hexagon)
            n += 1
            
    # Add a colorbar
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)

def get_td_power(filename):
    """ Get the number of cells form a .out file

    """
    # size = get_mesh_size(filename)
    power = parse_file(filename, begin='Total Power vector', n_max_lines=1)
    return power

def get_td_time(filename):
    """ Get the number of cells form a .out file

    """
    # size = get_mesh_size(filename)
    time = parse_file(filename, begin='Time vector', n_max_lines=1)
    return time

    
def collapse(points, flux):
    """ Collapse to the x direction
    >>> x = np.array([0., 1.,  2., 0., 1., 2., 0., 1., 2.,  0.])
    >>> y = np.array([0., 1.,  2., 0., 1., 2., 1., 2., 3.,  4.])
    >>> z = np.array([0., 0.,  0., 0., 0., 0., 0., 0., 0.,  0.])
    >>> points = [x, y, z]
    >>> flux = np.array([1.0, 2.0, 1.0, 1.0, 2.0, 1.0, 1.0, 3.0, 6.0, 7.0])
    >>> collapse(points, flux)
    ([0.0, 1.0, 2.0], [3.0, 2.5, 3.5])
    """
    # Remove repeated points
    n_points = len(flux)
    points_set = set()
    flux2 = []
    x = []
    for i in range(n_points):
        p = (points[0][i], points[1][i], points[2][i])
        if p not in points_set:
            points_set.add(p)
            flux2.append(flux[i])
            x.append(points[0][i])
    
    # Collapse
    n_points = len(flux2)
    values = dict()
    num_times = dict()
    for i in range(n_points):
        values[x[i]] = values.get(x[i], 0.0)  + flux2[i]
        num_times[x[i]] = num_times.get(x[i], 0)  + 1
    for k,val in num_times.items():
        values[k] /=  val
    
    # Sort
    fl = [f for _,f in sorted(values.items())]
    x = sorted(values.keys())
    return x, fl

#def get_line(points, flux):
#    """ Collapse to the x direction
#    >>> x = np.array([0., 1.,  2., 0., 1., 2., 0., 1., 2.,  0.])
#    >>> y = np.array([0., 1.,  2., 0., 1., 2., 1., 2., 3.,  4.])
#    >>> z = np.array([0., 0.,  0., 0., 0., 0., 0., 0., 0.,  0.])
#    >>> points = [x, y, z]
#    >>> flux = np.array([1.0, 2.0, 1.0, 1.0, 2.0, 1.0, 1.0, 3.0, 6.0, 7.0])
#    >>> collapse(points, flux)
#    ([0.0, 1.0, 2.0], [3.0, 2.5, 3.5])
#    """
#
#    # Remove repeated points
#    n_points = len(flux)
#    points_set = set()
#    flux2 = []
#    x = []
#    for i in range(n_points):
#        p = (points[0][i], points[1][i], points[2][i])
#        if p not in points_set:
#            points_set.add(p)
#            flux2.append(flux[i])
#            x.append(points[0][i])
#    
#    # Collapse
#    n_points = len(flux2)
#    values = dict()
#    num_times = dict()
#    for i in range(n_points):
#        values[x[i]] = values.get(x[i], 0.0)  + flux2[i]
#        num_times[x[i]] = num_times.get(x[i], 0)  + 1
#    for i in range(len(num_times)):
#        values[x[i]] =  values[x[i]] / num_times[x[i]]
#    
#    print num_times
#    # Sort
#    x = sorted(values.keys())
#    fl = [f for _,f in sorted(values.items())]
#    return x, fl



def parse_fluxes_parcs_out_1D(fileout):
    """ """
    ##########################################################################
    # Get n_assemblies  through  geo_dim sizes in Parcs
    #n_assemblies = 0
    head =  'geo_dim'

    with open(fileout) as fp:  
       for line in fp:
           line =  line.strip()
           if (line[0:len(head)] == head):  
               line = line.split()
               n_assemblies = int(line[1])
               break
                

    
    ##########################################################################
    # Get x through  grid sizes in Parcs
    head =  'grid_x'
    delta_x = []
    with open(fileout) as fp:  
       for line in fp:
           line =  line.strip()
           if (line[0:len(head)] == head):  
               line = line.split()[1:]

               for l in line:
                   l = l.split('*', 2)
                   if len(l)== 1:
                       delta_x.append(float(l[0]))
                   elif len(l)== 2:
                       for i in range(int(l[0])):
                           delta_x.append(float(l[1]))
                   else:
                      print('Invalid grid' + l )
                      assert(False)
               break
    
                
    assert(len(delta_x) == n_assemblies)            
    
    delta_x = np.array(delta_x)
    x = []
    for i in range(len(delta_x)):
        x.append(sum(delta_x[0:i]) + delta_x[i]/2)
    
    
    ##########################################################################
    # Get Flux shapes
    head =  ' At Simulation Time'
    head2 = ' Planar Flux at Plane'
    flux_fs = []
    flux_th = []
    t = []
    with open(fileout) as fp:  
       for line in fp:
    #        print(line)
            if (line[0:len(head)] == head):
                line = line.split()
                for i, l in enumerate(line):
                    if l[-1] == '=':
                        t.append(float(line[i+1]))
                
            if (line[0:len(head2)] == head2):
                fp.readline()
                fp.readline()
                # Fast Flux
                line = fp.readline()
                line = line.split()
                line = [float(i) for i in line]
                flux_fs.append(np.array(line[1:]))
                # Thermal Flux
                line = fp.readline()
                line = line.split()
                line = [float(i) for i in line]
                flux_th.append(np.array(line[0:]))
    return x, flux_fs, flux_th



def parse_line_fluxes_parcs_2D(fileout, n_line):
    """Parse an .out from parcs getting the x, y and the time dependent fluxes.
    >>> file = 'tests/test_2D_parcs.test' 
    >>> x, flux_fs, flux_th = parse_line_fluxes_parcs_2D(file, 2)
    >>> x
    array([ 5., 15., 25.])
    >>> flux_fs
    [array([1., 1., 1.]), array([1., 1., 1.]), array([1., 1., 1.])]
    >>> flux_th
    [array([3., 3., 3.]), array([3., 3., 3.]), array([3., 3., 3.])]
    """
    ##########################################################################
    # Get n_assemblies  through  geo_dim sizes in Parcs
    #n_assemblies = 0
    head =  'geo_dim'

    with open(fileout) as fp:  
       for line in fp:
           line =  line.strip()
           if (line[0:len(head)] == head):  
               line = line.split()
               n_assemblies_x = int(line[1])
               break
                    
    ##########################################################################
    # Get x through  grid sizes in Parcs
    head =  'grid_x'
    delta_x = []
    with open(fileout) as fp:  
       for line in fp:
           line =  line.strip()
           if (line[0:len(head)] == head):  
               line = line.split()[1:]

               for l in line:
                   l = l.split('*', 2)
                   if len(l)== 1:
                       delta_x.append(float(l[0]))
                   elif len(l)== 2:
                       for i in range(int(l[0])):
                           delta_x.append(float(l[1]))
                   else:
                      print('Invalid grid' + l )
                      assert(False)
               break
    
                
    assert(len(delta_x) == n_assemblies_x)            
    
    delta_x = np.array(delta_x)
    x = []
    for i in range(len(delta_x)):
        x.append(sum(delta_x[0:i]) + delta_x[i]/2)
    x = np.array(x)
    
    ##########################################################################
    # Get Flux shapes
    head = ' Planar Flux at Plane'
    flux_fs = []
    flux_th = []
    with open(fileout) as fp:  
       for line in fp:
            if (line[0:len(head)] == head):
                fp.readline()
                fp.readline()
                # Fast Flux
                for line in range(3*(n_line-1)):
                    line = fp.readline()

                line = fp.readline()

                line = line.split()
                line = [float(i) for i in line]
                flux_fs.append(np.array(line[1:]))
                # Thermal Flux
                line = fp.readline()
                line = line.split()
                line = [float(i) for i in line]
                flux_th.append(np.array(line[0:]))
                
    return x, flux_fs, flux_th

def parse_parcs_geometry(fileout):
    """Parse the geometry from a parcs .out file.
    >>> x, y = parse_parcs_geometry('tests/test_2D_parcs.test')
    >>> x
    array([ 5., 15., 25.])
    >>> y
    array([ 5., 15., 25.])
    """
    ##########################################################################
    # Get n_assemblies  through  geo_dim sizes in Parcs
    #n_assemblies = 0
    head =  'geo_dim'

    with open(fileout) as fp:  
       for line in fp:
           line =  line.strip()
           if (line[0:len(head)] == head):  
               line = line.split()
               n_assemblies_x = int(line[1])
               n_assemblies_y = int(line[2])
               #n_assemblies_z = int(line[3])
               break
                    
    ##########################################################################
    # Get x through  grid sizes in Parcs
    head =  'grid_x'
    delta_x = []
    with open(fileout) as fp:  
       for line in fp:
           line =  line.strip()
           if (line[0:len(head)] == head):  
               line = line.split()[1:]

               for l in line:
                   l = l.split('*', 2)
                   if len(l)== 1:
                       delta_x.append(float(l[0]))
                   elif len(l)== 2:
                       for i in range(int(l[0])):
                           delta_x.append(float(l[1]))
                   else:
                      print('Invalid grid x' + l )
                      assert(False)
               break
    
                
    assert(len(delta_x) == n_assemblies_x)            
    
    delta_x = np.array(delta_x)
    x = []
    for i in range(len(delta_x)):
        x.append(sum(delta_x[0:i]) + delta_x[i]/2)
    x = np.array(x)
        
    ##########################################################################
    # Get x through  grid sizes in Parcs
    head =  'grid_y'
    delta_y = []
    with open(fileout) as fp:  
       for line in fp:
           line =  line.strip()
           if (line[0:len(head)] == head):  
               line = line.split()[1:]

               for l in line:
                   l = l.split('*', 2)
                   if len(l)== 1:
                       delta_y.append(float(l[0]))
                   elif len(l)== 2:
                       for i in range(int(l[0])):
                           delta_y.append(float(l[1]))
                   else:
                      print('Invalid grid y' + l )
                      assert(False)
               break
               
    assert(len(delta_y) == n_assemblies_y)            
    
    delta_y = np.array(delta_y)
    y = []
    for i in range(len(delta_y)):
        y.append(sum(delta_y[0:i]) + delta_y[i]/2)
    y = np.array(y)
            
    return x, y

def parse_parcs_time(fileout):
    """Parse the geometry from a parcs .out file.
    >>> t = parse_parcs_time('tests/test_2D_parcs.test')
    >>> t
    array([0.   , 0.001, 0.002])
    """
    ##########################################################################
    # Get n_assemblies  through  geo_dim sizes in Parcs
    head= ' At Simulation Time = '

    t = [0.0]

    with open(fileout) as fp:  
       for line in fp:         
            if (line[0:len(head)] == head):
                line = line.split()
                t.append(float(line[4]))
    t = np.array(t)
    return t

    
def parse_fluxes_parcs_2D(fileout):
    """Parse an .out from parcs getting the x, y and the time dependent fluxes.
    >>> flux_fs, flux_th = parse_fluxes_parcs_2D('tests/test_2D_parcs.test')
    >>> flux_fs
    array([[[1., 1., 1.],
            [1., 1., 1.],
            [1., 1., 8.]],
    <BLANKLINE>
           [[1., 1., 1.],
            [1., 1., 1.],
            [1., 1., 6.]],
    <BLANKLINE>
           [[1., 1., 1.],
            [1., 1., 1.],
            [1., 1., 1.]]])
    >>> flux_th
    array([[[3., 1., 3.],
            [3., 3., 3.],
            [3., 7., 3.]],
    <BLANKLINE>
           [[3., 1., 3.],
            [3., 3., 3.],
            [3., 3., 3.]],
    <BLANKLINE>
           [[3., 1., 3.],
            [3., 3., 3.],
            [3., 3., 6.]]])
    """
    ##########################################################################
    # Get n_assemblies  through  geo_dim sizes in Parcs
    #n_assemblies = 0
    head =  'geo_dim'

    with open(fileout) as fp:  
       for line in fp:
           line =  line.strip()
           if (line[0:len(head)] == head):  
               line = line.split()
               n_assemblies_x = int(line[1])
               n_assemblies_y = int(line[2])
#               n_assemblies_z = int(line[3])
               break

    ##########################################################################
    # Get n_assemblies  through  geo_dim sizes in Parcs
    #n_assemblies = 0
    head =  'rad_conf'
    rad_conf = []
    with open(fileout) as fp:  
       for line in fp:
           line =  line.strip()
           if (line[0:len(head)] == head):
               for line in range(n_assemblies_y):
                   line = fp.readline()
                   line = line.split()
                   line = [int(i) for i in line]
                   rad_conf.append(line[0:])
               break
           
    #print(rad_conf)
    ##########################################################################
    # Get Flux shapes
    head_flux = ' Planar Flux at Plane'
    flux_fs = []
    flux_th = []

    with open(fileout) as fp:  
       for line in fp:         
            if (line[0:len(head_flux)] == head_flux):
                flux_fs_line = []
                flux_th_line = []
                fp.readline()
                fp.readline()
                for line in range(n_assemblies_y):
                    # Fast Flux              
                    line = fp.readline()
                    line = line.split()
                    line = [float(i) for i in line]
                    flux_fs_line.append(line[1:])
                    
                    # Thermal Flux
                    line = fp.readline()
                    line = line.split()
                    line = [float(i) for i in line]
                    flux_th_line.append(line[0:])
                    
                    # Blank line
                    line = fp.readline()
                

                flux_fs.append(np.array(flux_fs_line))
                flux_th.append(np.array(flux_th_line)) 
                
    ##########################################################################
    n_steps = len(flux_fs)
    
    flux_fs2 = np.zeros([n_steps, n_assemblies_y, n_assemblies_x])
    flux_th2 = np.zeros([n_steps, n_assemblies_y, n_assemblies_x])
    # Put 
    for st in range(n_steps):
        for ny in range(n_assemblies_y):
            nx1 = 0
            for nx in range(n_assemblies_x):
                if rad_conf[ny][nx] == 0:
                    
                    flux_fs2[st][ny][nx] = 0.0
                    flux_th2[st][ny][nx] = 0.0
                else:
                    flux_fs2[st][ny][nx] = flux_fs[st][ny][nx1]
                    flux_th2[st][ny][nx] = flux_th[st][ny][nx1]
                    nx1 += 1 
    
    return flux_fs2, flux_th2



def parse_fluxes_parcs_static(fileout):
    """Get the static fluxes from a PARCS .out.
    >>> flux_fs, flux_th = parse_fluxes_parcs_static('tests/test_2D_parcs.test')
    >>> flux_fs
    array([1., 1., 1., 1., 1., 1., 1., 1., 8.])
    >>> flux_th
    array([3., 1., 3., 3., 3., 3., 3., 7., 3.])
    """
    ##########################################################################
    # Get n_assemblies  through  geo_dim sizes in Parcs
    #n_assemblies = 0
    head =  'geo_dim'

    with open(fileout) as fp:  
       for line in fp:
           line =  line.strip()
           if (line[0:len(head)] == head):  
               line = line.split()
               n_assemblies_y = int(line[2])
               n_assemblies_z = int(line[3])
               break
           
    ##########################################################################
    # Get Flux shapes
    head2 = ' Planar Flux at Plane'
    flux_fs = []
    flux_th = []
    with open(fileout) as fp:  
       for line in fp:          
            if (line[0:len(head2)] == head2):
                fp.readline()
                fp.readline()
                for line in range(n_assemblies_z):
                    for line in range(n_assemblies_y):
                        # Fast Flux
                        line = fp.readline()
                        line = line.split()
                        line = [float(i) for i in line]
                        for flx in range(1, len(line)):
                            flux_fs.append(line[flx])
                        
                        # Thermal Flux
                        line = fp.readline()
                        line = line.split()
                        line = [float(i) for i in line]
                        for flx in line:
                            flux_th.append(flx)
                            
                        # Blank line
                        line = fp.readline()
                break
            
    flux_fs = np.array(flux_fs)
    flux_th = np.array(flux_th)
    
    return flux_fs, flux_th


def get_parcs_keff(fileout):
    """Get keff form a parcs .out file.
    >>> get_parcs_keff('tests/test_2D_parcs.test')
    1.0
    >>> get_parcs_keff('tests/test_3D_parcs.test')
    1.053634
    """
    return parse_file_same_line(fileout, '   K-Effective: ')[0]

def parse_parcs_fluxes_3D(fileout):
    """Parse an .out from parcs getting the x, y and the time dependent fluxes.
    >>> flux_fs, flux_th = parse_parcs_fluxes_3D('tests/test_3D_parcs.test')
    >>> flux_fs
    array([[[ 1.  ,  2.  ,  3.  ],
            [ 4.  ,  5.27,  6.27],
            [ 7.27,  8.27,  9.27]],
    <BLANKLINE>
           [[10.27, 11.27, 12.27],
            [13.27, 14.27, 15.27],
            [16.27, 17.27, 18.27]],
    <BLANKLINE>
           [[19.27, 20.27, 21.27],
            [22.27, 23.27, 24.27],
            [25.27, 26.27, 27.28]]])
    >>> flux_th
    array([[[ 5.  ,  6.  , 77.  ],
            [ 7.  ,  0.56, 10.56],
            [ 8.56, 10.56, 10.56]],
    <BLANKLINE>
           [[10.56, 10.56, 10.56],
            [10.56, 10.56, 10.56],
            [10.56, 10.56, 10.56]],
    <BLANKLINE>
           [[10.56, 10.56, 10.56],
            [10.56, 10.56, 10.56],
            [10.56, 10.56, 10.56]]])
    >>> flux_fs[0][0][2]
    3.0
    >>> flux_fs[0][1][2]
    6.27
    >>> flux_fs[1][0][0]
    10.27
    """
    ##########################################################################
    # Get n_assemblies  through  geo_dim sizes in Parcs
    #n_assemblies = 0
    head =  'geo_dim'

    with open(fileout) as fp:  
       for line in fp:
           line =  line.strip()
           if (line[0:len(head)] == head):  
               line = line.split()
               n_assemblies_x = int(line[1])
               n_assemblies_y = int(line[2])
#               n_assemblies_z = int(line[3])
               break

    ##########################################################################
    # Get n_assemblies  through  geo_dim sizes in Parcs
    #n_assemblies = 0
    head =  'rad_conf'
    rad_conf = []
    with open(fileout) as fp:  
       for line in fp:
           line =  line.strip()
           if (line[0:len(head)] == head):
               for line in range(n_assemblies_y):
                   line = fp.readline()
                   line = line.split()
                   line = [int(i) for i in line]
                   rad_conf.append(line[0:])
               break
           
    #print(rad_conf)
    ##########################################################################
    # Get Flux shapes
    head_flux = ' Planar Flux at Plane'
    flux_fs = []
    flux_th = []

    with open(fileout) as fp:  
       for line in fp:         
            if (line[0:len(head_flux)] == head_flux):
                flux_fs_line = []
                flux_th_line = []
                fp.readline()
                fp.readline()
                for line in range(n_assemblies_y):
                    # Fast Flux              
                    line = fp.readline()
                    line = line.split()
                    line = [float(i) for i in line]
                    flux_fs_line.append(line[1:])
                    
                    # Thermal Flux
                    line = fp.readline()
                    line = line.split()
                    line = [float(i) for i in line]
                    flux_th_line.append(line[0:])
                    
                    # Blank line
                    line = fp.readline()
                

                flux_fs.append(np.array(flux_fs_line))
                flux_th.append(np.array(flux_th_line)) 
                
    ##########################################################################
    n_steps = len(flux_fs)
    
    flux_fs2 = np.zeros([n_steps, n_assemblies_y, n_assemblies_x])
    flux_th2 = np.zeros([n_steps, n_assemblies_y, n_assemblies_x])
    # Put 
    for st in range(n_steps):
        for ny in range(n_assemblies_y):
            nx1 = 0
            for nx in range(n_assemblies_x):
                if rad_conf[ny][nx] == 0:
                    
                    flux_fs2[st][ny][nx] = 0.0
                    flux_th2[st][ny][nx] = 0.0
                else:
                    flux_fs2[st][ny][nx] = flux_fs[st][ny][nx1]
                    flux_th2[st][ny][nx] = flux_th[st][ny][nx1]
                    nx1 += 1 
    
    return flux_fs2, flux_th2

def remove_repeated_data_point(x,y,z, data):
    """ Remove repeated data
    """
    dict_data = {}
    for p in range(len(data)):
        dict_data[(x[p], y[p], z[p])] = data[p]
    return np.array(list(dict_data.values()))

def remove_repeated_point(x,y,z):
    """ Remove repeated point
    """
    dict_x = {}
    dict_y = {}
    dict_z = {}    
    for p in range(len(x)):
        dict_x[(x[p], y[p], z[p])] = x[p]
        dict_y[(x[p], y[p], z[p])] = y[p]
        dict_z[(x[p], y[p], z[p])] = z[p]
        
    return (np.array(list(dict_x.values())),
            np.array(list(dict_y.values())), 
            np.array(list(dict_z.values())))


def  parse_vtk_over_line(file, scalar_name, point_a, point_b, resolution=1000) :
    """Sample a dataset onto a line.

    Parameters
    ----------
    pointa : sequence[float]
        Location in ``[x, y, z]``.

    pointb : sequence[float]
        Location in ``[x, y, z]``.
        
    scalar_name : 
        Name of scalar value, example ``phi_g1_eig_1``.

    resolution : int, optional
        Number of pieces to divide line into. Defaults to number of cells
        in the input mesh. Must be a positive integer.
        
    >>> point_a = (0, 4.5, 0)
    >>> point_b = (10.0, 4.5, 0)
    >>> x, y = parse_vtk_over_line("tests/2Dtest.vtk",'phi_g1_eig_1', point_a, point_b, resolution=2)
    >>> print(x)
    [ 0.  5. 10.]
    >>> print(y)
    [14.1851   12.3342    0.970584]
   """
    print('Not working as pyvista is not installed TODO')
    return
    # mesh = pv.read(file)
    # line = pv.Line(point_a, point_b, resolution=resolution)
    # line = line.sample(mesh)
    
    # # Get x
    # xa, ya, za = point_a
    # xb, yb, zb = point_b
    # distance = np.sqrt((xb - xa)**2 +  (yb - ya)**2 + (zb - zb)**2)
    # x = np.linspace(0, distance, resolution+1)
    
    # return x, np.array(line.get_array(scalar_name))



if __name__ == "__main__":
    import doctest
    doctest.testmod()
