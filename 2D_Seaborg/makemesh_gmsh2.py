# ------------------------------------------------------------------------------
# Make mesh
# ------------------------------------------------------------------------------

# The Python API is entirely defined in the `gmsh.py' module (which contains the
# full documentation of all the functions in the API):
import gmsh
import sys
import math

# Before using any functions Gmsh API must be initialized:
gmsh.initialize()
gmsh.option.setNumber("Mesh.Algorithm", 8)
gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", 2)
gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
gmsh.option.setNumber("Mesh.Smoothing", 10)

msh = gmsh.model.geo  # rename namespace

# Next we add a new model named "t1" (
gmsh.model.add("seaborg")
output_msh = "seaborg.msh"

pins = [1,2,1,2,1]

radius = 25  # cm
pitch = 66   # cm
centers = [[0, 0],     [-57.2,  33], [0,  66],
           [57.2, 33],  [57.2, -33], [0, -66],
           [-57.2,-33]]
lc = pitch/50.0  # cm

points =[
        [24,  19.05255888325765, 99],
        [26, -19.05255888325765, 99], 
        [13, -38.14744111674236, 66],
        [15, -76.25255888325765, 66],
        [16, -95.30511776651531, 33],
        [17, -76.25255888325765, 0],
        [71, -95.30511776651531, -33], 
        [72, -76.25255888325765, -66], 
        [60, -38.1051177665153, -66],
        [61, -19.05255888325765, -99], 
        [58,  19.05255888325765, -99],
        [50,  38.14744111674236, -66],
        [47,  76.25255888325765, -66],
        [45,  95.30511776651531, -33],
        [36,  76.25255888325765, 0],
        [34,  95.30511776651531, 33],
        [35,  76.25255888325765, 66],
        [23,  38.1051177665153, 66]]



def make_circle(cx, cy, radius, id):
    """Make a circle inside an  hexagon
    cx - Center of the circle and the hexagon
    cy - Center of the circle and the hexagon
    pitch - Pitch of the hexagonal
    radius - Radius of the cercle
    id - Id of the hexagon, the circle will have and if
    """
    # Circle
    # Inner Fuel Circle
    p0 = msh.addPoint(cx, cy, 0, lc)         # Center
    p7 = msh.addPoint(cx, cy - radius, 0, lc)  # Bottom
    p8 = msh.addPoint(cx + radius, cy, 0, lc)  # Right
    p9 = msh.addPoint(cx, cy + radius, 0, lc)  # Top
    p10 = msh.addPoint(cx - radius, cy, 0, lc)  # Left

    c1 = msh.addCircleArc(p7, p0, p8)
    c2 = msh.addCircleArc(p8, p0, p9)
    c3 = msh.addCircleArc(p9, p0, p10)
    c4 = msh.addCircleArc(p10, p0, p7)
    ll_circ = msh.addCurveLoop([c1, c2, c3, c4])

    circ = msh.addPlaneSurface([ll_circ])
    gmsh.model.addPhysicalGroup(2, [circ], id)
    
    return circ;


id = 1
circles = []

for cent in centers:
    c = make_circle(cent[0], cent[1], radius, id)
    circles.append(c)
    id += 1 
    
######################################################################
# HEXAGONAL BOUNDARY

points_ids = []
for p in points:
    points_ids.append(msh.addPoint(p[1], p[2], 0.0, lc))

bc_lines = []  
for i in range(len(points_ids)-1):
    bc_lines.append(msh.addLine(points_ids[i], points_ids[i+1]))
bc_lines.append(msh.addLine(points_ids[-1], points_ids[0]))

ll_hexagon = msh.addCurveLoop(bc_lines)
circles.insert(0 , ll_hexagon)

hexagon = msh.addPlaneSurface(circles)  
gmsh.model.addPhysicalGroup(2, [hexagon], 8)

gmsh.model.addPhysicalGroup(1, bc_lines, 1)
######################################################################


msh.synchronize()

# We can then generate a 2D mesh...
gmsh.model.mesh.generate()
gmsh.write(output_msh)


# To visualize the model we can run the graphical user interface with
# `gmsh.fltk.run()'. Here we run it only if the "-nopopup" is not provided in
# the command line arguments:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()


gmsh.finalize()
