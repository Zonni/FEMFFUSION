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
gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", 7.5)
gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
gmsh.option.setNumber("Mesh.Smoothing", 10)

msh = gmsh.model.geo  # rename namespace

# Next we add a new model named "t1" (
gmsh.model.add("seaborg")
output_msh = "seaborg.msh"

pins = [1,2,1,2,1]

radius = 25  # cm
pitch = 66   # cm
centers = [[0, 0],     [-57.1576,  33], [0,  66],
           [57.1576, 33],  [57.1576, -33], [0, -66],
           [-57.1576,-33]]
lc = pitch/100.0  # cm
points = []
lines = []

def new_point(px, py):
    """Create a new point and return its point id in case it does not exists
    Otherwise return its point id.
    """
    for p in points:
        dist = (p[1] - px)**2 + (p[2] - py)**2
        if (dist < 0.1):
            return p[0]
    
    point_id = msh.addPoint(px, py, 0, lc)  # Rigth
    points.append((point_id, px, py))
    return point_id

def new_line(p1, p2):
    """Create a new line and return its id in case it does not exists
    Otherwise return its point id.
    """
    for l in lines:
        if (l[1]== p1 and p2==l[2]): # Change by a dictionary
            return l[0]
    
    line_id = msh.addLine(p1, p2)  # Rigth
    lines.append((line_id, p1, p2))
    return line_id

def make_hexagon_and_circle(cx, cy, pitch, radius, id):
    """Make a circle inside an  hexagon
    cx - Center of the circle and the hexagon
    cy - Center of the circle and the hexagon
    pitch - Pitch of the hexagonal
    radius - Radius of the cercle
    id - Id of the hexagon, the circle will have and if
    """

    ri = pitch / 2        # Radius Incircle
    rc = 2.0/3.0*math.sqrt(3) * ri  # Radio Excircles == Side

    p1 = new_point(cx+rc,   cy+0.0) # Rigth
    p2 = new_point(cx+rc/2, cy+ri)  # Top Right
    p3 = new_point(cx-rc/2, cy+ri)  # Top Left
    p4 = new_point(cx-rc,  cy+0.0)  # Left
    p5 = new_point(cx-rc/2, cy-ri)  # Bottom Left
    p6 = new_point(cx+rc/2, cy-ri)  # Bottom Right

    l1 = new_line(p1, p2)  # Bottom Right
    l2 = new_line(p2, p3)  # Right
    l3 = new_line(p3, p4)  # Top Right
    l4 = new_line(p4, p5)  # Top Left
    l5 = new_line(p5, p6)  # Left
    l6 = new_line(p6, p1)  # Left Right

    ll_hexagon = msh.addCurveLoop([l1, l2, l3, l4, l5, l6])

    ################################################################
    # Circle
    # Inner Fuel Circle
    p0 = new_point(cx, cy)           # Center
    p7 = new_point(cx, cy - radius)  # Bottom
    p8 = new_point(cx + radius, cy)  # Right
    p9 = new_point(cx, cy + radius)  # Top
    p10 = new_point(cx - radius, cy) # Left

    c1 = msh.addCircleArc(p7,  p0, p8)
    c2 = msh.addCircleArc(p8,  p0, p9)
    c3 = msh.addCircleArc(p9,  p0, p10)
    c4 = msh.addCircleArc(p10, p0, p7)
    ll_circ = msh.addCurveLoop([c1, c2, c3, c4])

    circ = msh.addPlaneSurface([ll_circ])
    gmsh.model.addPhysicalGroup(2, [circ], id)

    ################################################################
    hexagon = msh.addPlaneSurface([ll_hexagon, ll_circ])  # hexagon with a hole
    gmsh.model.addPhysicalGroup(2, [hexagon], id+1)
    
    print(f"hexagon at center ({cx}, {cy}) with id {id} with lines")
    print([p1, p2, p3, p4, p5, p6])




id = 1
for cent in centers:
    make_hexagon_and_circle(cent[0], cent[1], pitch, radius, id)
    id += 2

# Boundary
#bc_lines = [12, 13, 14, 63, 64, 65, 54, 55,
#            56, 45, 46, 41, 36, 31, 32, 21, 22, 23]
#gmsh.model.addPhysicalGroup(1, bc_lines, 1)

msh.synchronize()

# We can then generate a 2D mesh...
gmsh.model.mesh.generate()
gmsh.write(output_msh)

# To visualize the model:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()


gmsh.finalize()