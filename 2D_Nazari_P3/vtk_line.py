"""
.. _plot_over_line_example:

Plot Over Line
~~~~~~~~~~~~~~

Plot the values of a dataset over a line through that dataset
"""

# sphinx_gallery_thumbnail_number = 2
import pyvista as pv
from pyvista import examples


file = 'P3_SDP1_Nazari.out.vtk'
mesh = pv.read(file)

file2 = 'P3_SDP1.out.vtk'
mesh2 = pv.read(file2)


###############################################################################
# Flat Surface
# ++++++++++++
import collections.abc
from typing import Optional, Sequence, Union
import warnings

import matplotlib.pyplot as plt
import numpy as np

import pyvista
from pyvista import FieldAssociation, _vtk
from pyvista.core.errors import VTKVersionError
from pyvista.core.filters import _get_output, _update_alg
from pyvista.errors import AmbiguousDataError, MissingDataError
from pyvista.utilities import (
    NORMALS,
    abstract_class,
    assert_empty_kwargs,
    generate_plane,
    get_array,
    get_array_association,
    transformations,
    wrap,
)
from pyvista.utilities.cells import numpy_to_idarr

# Make two points to construct the line between
pointa = [0, 8.967, 0.0]
pointb = [10.0, 8.967, 0.0]


###############################################################################
# # Run the filter and produce a line plot
# mesh.plot_over_line(
#     a,
#     b,
#     resolution=10000,
#     title="Elevation Profile",
#     ylabel="Height above sea level",
#     figsize=(10, 5),
#     figure=False
# )



def sample_over_line(self, pointa, pointb, resolution=None, tolerance=None, progress_bar=False):
    """Sample a dataset onto a line.

    Parameters
    ----------
    pointa : sequence[float]
        Location in ``[x, y, z]``.

    pointb : sequence[float]
        Location in ``[x, y, z]``.

    resolution : int, optional
        Number of pieces to divide line into. Defaults to number of cells
        in the input mesh. Must be a positive integer.

    tolerance : float, optional
        Tolerance used to compute whether a point in the source is in a
        cell of the input.  If not given, tolerance is automatically generated.

    progress_bar : bool, default: False
        Display a progress bar to indicate progress.

    Returns
    -------
    pyvista.PolyData
        Line object with sampled data from dataset.

    Examples
    --------
    Sample over a plane that is interpolating a point cloud.

    >>> import pyvista
    >>> import numpy as np
    >>> np.random.seed(12)
    >>> point_cloud = np.random.random((5, 3))
    >>> point_cloud[:, 2] = 0
    >>> point_cloud -= point_cloud.mean(0)
    >>> pdata = pyvista.PolyData(point_cloud)
    >>> pdata['values'] = np.random.random(5)
    >>> plane = pyvista.Plane()
    >>> plane.clear_data()
    >>> plane = plane.interpolate(pdata, sharpness=3.5)
    >>> sample = plane.sample_over_line((-0.5, -0.5, 0), (0.5, 0.5, 0))
    >>> pl = pyvista.Plotter()
    >>> _ = pl.add_mesh(
    ...     pdata, render_points_as_spheres=True, point_size=50
    ... )
    >>> _ = pl.add_mesh(sample, scalars='values', line_width=10)
    >>> _ = pl.add_mesh(plane, scalars='values', style='wireframe')
    >>> pl.show()

    """
    if resolution is None:
        resolution = int(self.n_cells)
    # Make a line and sample the dataset
    line = pyvista.Line(pointa, pointb, resolution=resolution)
    sampled_line = line.sample(self, tolerance=tolerance, progress_bar=progress_bar)
    return sampled_line



mesh2.plot_over_line(
#     a,
#     b,
#     resolution=10000,
#     title="y=" + str(y_line) + " cm",
#     ylabel="Scalar Flux (AU)"
#     xlabel=" x (cm)",
#     figsize=(10, 5),
#     figure=False,
# )
tolerance = None
resolution = 100
progress_bar = False

# Sample on line
sampled = DataSetFilters.sample_over_line(
    mesh2, pointa, pointb, resolution, tolerance, progress_bar=progress_bar
)

# Get variable of interest
if scalars is None:
    pyvista.set_default_active_scalars(mesh2)
    field, scalars = mesh2.active_scalars_info
values = sampled.get_array(scalars)
distance = sampled['Distance']

for i in range(values.shape[1]):
    print(values[:, i])
# Plot it in 2D
# for i in range(values.shape[1]):
#     plt.plot(distance, values[:, i], label=f'Component {i}')






