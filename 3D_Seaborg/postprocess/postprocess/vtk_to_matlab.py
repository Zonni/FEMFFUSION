# -*- coding: utf-8 -*-
"""
We use this script for ploting 1d vtk files with matplotlib instead
of paraview. We first use the `vtk` library to read the data,
and then we convert the data to `numpy` vector, to be ploted with
`matplotlib`.

You have to install python-vtk:
  sudo apt-get install python-vtk

The package is documented in
  http://www.vtk.org/doc/nightly/html/index.html

The function that we are using is documented in
  http://www.vtk.org/doc/nightly/html/classvtkUnstructuredGridReader.html

And the data object returned is documented here
  http://www.vtk.org/doc/nightly/html/classvtkUnstructuredGrid.html
"""

from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN

def parseVtkFile(filename, title):
    """ Parse a File and return what is below a begin title.
    It must be defined a the begin title and the end word or/and the maximum
    number of lines to read.
    >>> filename = 'tests/vtk.test'
    >>> parseVtkFile(filename, 'Neutron_Power')
    array([ 1.,  1.,  1.,  1.])
    """
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    array = VN.vtk_to_numpy(data.GetPointData().GetArray(title))
    return array
    
def parseVtkGrid(filename):
    """ Return the x, y, z np.arrays of the grid given in a vtk file
    >>> filename = 'tests/vtk.test'
    >>> parseVtkGrid(filename 'Neutron_Power')
    array([ 1.,  1.,  1.,  1.])
    """
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    grid = []

    for i in range(data.GetNumberOfPoints()):
        grid.append(data.GetPoint(i))
    return grid
    
def vectoToMatlab(file_handler, vector_name, vector):
    file_handler.write(vector_name + ' = ' + str(vector) + ';\n')
    
    
def removeDuplicates(grid, values):
    output_grid = []
    output_values = []
    for j in range(len(values)):
        output_values.append([])
        
    for i in range(len(grid)):
        p = grid[i]
        if p not in output_grid:
            output_grid.append(p)
            for j in range(len(values)):
                output_values[j].append(values[j][i])
            
    return output_grid, output_values
    


## Parameters
#vtk_file = "reactores/1D12.out.vtk"
#file_matlab = '1D12.m'
#
#f1=open(file_matlab, 'w+')
## parse from vtk file
#grid = parseVtkGrid(vtk_file)
#phi1 = parseVtkFile(vtk_file, "Phi1_Eig_1")
#phi2 = parseVtkFile(vtk_file, "Phi2_Eig_1")
## Remove duplicates in colliding cells
#grid, [phi1, phi2] = removeDuplicates(grid, [phi1, phi2])
#
## unpack grid
#x = []
#y = []
#z = []
#for xi, yi, zi in grid:
#    x.append(xi)
#    y.append(yi)
#    z.append(zi)
#
## Print to matlab file
#vectoToMatlab(f1, 'phi1', phi1)
#vectoToMatlab(f1, 'phi2', phi2)
#vectoToMatlab(f1, 'x', x)
#vectoToMatlab(f1, 'y', y)
#vectoToMatlab(f1, 'z', z)
#f1.close()
#

