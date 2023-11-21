import sys
import os
import gmsh  #header file for the python API
import math

#Use the below as a main reference (and the gmsh API docs)
#https://gitlab.onelab.info/gmsh/gmsh/-/blob/gmsh_4_10_5/tutorials/python/t16.py

#-----initialisation-----------------------------------------------------------#
gmsh.initialize(sys.argv)
#clear any pre-existing file of the same name
filename = "FullCoreReflected2D"
if os.path.isfile(filename + ".msh"):
    print("deleting old .msh")
    os.remove(filename + ".msh")

gmsh.open(filename + ".msh") #file must be opened prior to defining anything
gmsh.model.add(filename)

gmsh.logger.start()

#-----geometry-----------------------------------------------------------------#

#remember gmshToFoam assumes dimensions are in metres

#general hexagon dimensions
m = 1
cm2m = 1.0
#---target mesh size in vicinity of points
lc = 10*cm2m #2.5*cm2m #1.25*cm2m

#General geometry parameters
hex_halfpitch = 7.5*cm2m #half pitch defined as shortest distance between opposite sides
hex_fullpitch = 2.0*hex_halfpitch
hex_sidelength = 2.0/math.sqrt(3)*hex_halfpitch
row_dist = 1.5*hex_sidelength
tube_radius = 5.0*cm2m
tube_thickness = 0.5*cm2m
core_height = 100.0*cm2m #Arbitrary height; just meant for 2D case in openfoam
rod_radius = 1.0*cm2m #All absorber rods have the same size


#Copy the hexcentres dimensions from the Serpent Deck Generator script:
#%%%%%%%%%%%%%%%%
N_rows = 13 #Just define the number of rows in the full lattice
core_width = N_rows*2*hex_halfpitch
reflector_min_thickness = 30.0*cm2m #30 cm thick at its thinnest point
#reflector_radius = round((core_width/2.0 + reflector_min_thickness)*100.0)*cm2m #should give a radius of 127 cm 
reflector_radius = 127.0
print("TOTAL RADIUS OF THE MODELLED CORE:")
print(reflector_radius)
#%%%%%%%%%%%%%%%%
N_hexes_per_row = []
for i in range((N_rows+1)//2):
    N_hexes_per_row.append( (N_rows+1)//2 + i )
for i in range((N_rows-1)//2):
    N_hexes_per_row.append( (N_rows-1) - i )
#The above makes this list but for any max N_rows ---- 
#N_hexes_per_row = [7,8,9,10,11,12,13,12,11,10,9,8,7]

print("number of hex cells per row: " )
print(N_hexes_per_row)
print("--------------")
#define lists for hexes and fuel tubes. The fuel tube elements should be lists of the tube layers
#List of lists of coordinates of each hex, each row also in its own list
hexcentres = []
for n in range(N_rows):
     hexcentres.append([]) #prep the x and y coords of each hex

for i,N in enumerate(N_hexes_per_row):
    for n in range(N):
        hexcentres[i].append([0,0])
    #print(hexcentres[i])
#Now have the set-up for coords of all hex cells. Just need to populate.
#Treating the centre of the full core as the origin, loop through all the hexes row by row from bottom to top
for i,N in enumerate(N_hexes_per_row):
    for n in range(N):
        if N%2 > 0: #if the row has an odd number of hexes
            #works for both positive and negative cases
            hexcentres[i][n][0] = -((N-1)//2 - n)*hex_fullpitch #x-coord
            hexcentres[i][n][1] = -((N_rows-1)//2 - i)*row_dist #y-coord
        else: #else if the number of hexes is even
            hexcentres[i][n][0] = -((N-1)//2 - n + 0.5)*hex_fullpitch #x-coord
            hexcentres[i][n][1] = -((N_rows-1)//2 - i)*row_dist #y-coord

    #print("row number " + str(i+1))
    #print(hexcentres[i])


#Define the tubes first
#Each tube is defined as an extruded circle with 4 points around its perimeter.
tubes = []#list of tubes containing all points that define each tube
tubelines = []
tubecurveloops = []
tubesurfs = []
#cyl_radii = [tube_radius]
#cyl_radii = [tube_radius*3.0/4.0, tube_radius-tube_thickness, tube_radius] #can add more annuli if desired
cyl_radii = [tube_radius-tube_thickness, tube_radius] #can add more annuli if desired
#Obtain all the coordinates for the relevant points first
tubecount = 0
for row in hexcentres:
    for hex in row:
        tubes.append([])
        for r in cyl_radii:
            points = [ [hex[0], hex[1]], [ hex[0], hex[1] - r ], [ hex[0] + r, hex[1] ], [ hex[0],hex[1] + r ], [ hex[0] - r, hex[1] ]  ]
            tubes[tubecount].append(points)
        tubecount = tubecount + 1
#Go ahead and make the circles that correspond to each tube
#print(len(tubes[0]))
#print(len(tubes[-1]))
tubesurfs = [None]*len(cyl_radii)*sum(N_hexes_per_row) #prepare the list to hold tubesurfs in the right order
for i in range(len(tubes)):
    tubecurveloops.append([])
    tubelines.append([])
    #tubesurfs.append([])
    #tubes[i] contains all the points for all annuli at a certain tube
    for j,points in enumerate(tubes[i]): #5 pts at each given radius at index j
        #print(points)
        for k,point in enumerate(points):
            #print(point)
            #print(j)
            #print(k)
            gmsh.model.geo.addPoint(point[0],point[1],0, lc, 5*len(cyl_radii)*i + 5*j + k+1 ) #len(cyl_radii) = len(tubes[i])
            #add the index of each point as an added "coordinate"
            tubes[i][j][k].append(5*len(cyl_radii)*i + 5*j + k+1)
        #connect the points in arcs
        gmsh.model.geo.addCircleArc(tubes[i][j][1][2],tubes[i][j][0][2],tubes[i][j][2][2], 4*len(cyl_radii)*i +4*j+1) #start point, centrepoint, endpoint
        gmsh.model.geo.addCircleArc(tubes[i][j][2][2],tubes[i][j][0][2],tubes[i][j][3][2], 4*len(cyl_radii)*i +4*j+2)
        gmsh.model.geo.addCircleArc(tubes[i][j][3][2],tubes[i][j][0][2],tubes[i][j][4][2], 4*len(cyl_radii)*i +4*j+3)
        gmsh.model.geo.addCircleArc(tubes[i][j][4][2],tubes[i][j][0][2],tubes[i][j][1][2], 4*len(cyl_radii)*i +4*j+4)
        #Save the arcs (lines) forming each tube in a separate list
        tubelines[i].append([4*len(cyl_radii)*i +4*j+1, 4*len(cyl_radii)*i +4*j+2, 4*len(cyl_radii)*i +4*j+3, 4*len(cyl_radii)*i +4*j+4])
        #add the curve loop and plane surface

        temp_tubecurveloop = gmsh.model.geo.addCurveLoop(tubelines[i][-1], i*len(cyl_radii)+j+1) #return the object
        tubecurveloops[i].append(temp_tubecurveloop)
        if j==0: #if innermost fuel region
            temp_tubesurf = gmsh.model.geo.addPlaneSurface([i*len(cyl_radii)+j+1], i*len(cyl_radii)+j+1) #remember all gmsh objects are independently numbered
            #tubesurfs[i] = temp_tubesurf #inner fuel regions first
            #tubesurfs.append(temp_tubesurf)
        else:
            temp_tubesurf = gmsh.model.geo.addPlaneSurface([i*len(cyl_radii)+j, -(i*len(cyl_radii)+j+1)], i*len(cyl_radii)+j+1) #make the surf between this curveloop and the last

        tubesurfs[j*sum(N_hexes_per_row) + i] = temp_tubesurf #when j = 0, the first 127 indices of tubesurfs get populated. Then the next 127, and so on.
            #tubesurfs.append(temp_tubesurf)


print("first point at first tube, first radius, ")
print(tubes[0][0][0])
print("second point at first tube, first radius, ")
print(tubes[0][0][1])
print(len(tubes[0][0]))
print(tubes[1][0][0])

#Hexes now. Define the vertices as points without duplicating.
#additions to coords of hexcentres to get all points of a hex
# p0 = [0.0, 0.0]
p1 = [0.0, -hex_sidelength]
p2 = [+hex_halfpitch, -hex_sidelength/2.]
p3 = [+hex_halfpitch, +hex_sidelength/2.]
p4 = [0.0,  +hex_sidelength]
p5 = [-hex_halfpitch, +hex_sidelength/2.]
p6 = [-hex_halfpitch, -hex_sidelength/2.]
hex_part = [p1, p2, p3, p4, p5, p6]
#Again obtain the coordinates first
hexes = [] #build a basic list of lists of all the coordinates of all hex vertices.
for row in hexcentres:
    for hex in row:
        hexverts = []
        for p in hex_part:
            vert_coords = [sum(x) for x in zip(p,hex)] #add elementwise to get local node coords
            #if node_coords not in points:
            hexverts.append(vert_coords)
        hexes.append(hexverts)
#Now iterate through the "points" list and build the hexnodes list/dict
nodes = {} #dictionary linking unique point element labels to coordinates.
lines = {} #similarly need a dict for the lines forming all hexagons. Save as tuples of node element
hexnodes = [] #here save the nodes that form each hexagon in the same ordering as usual; duplicates allowed
hexlines = [] #save the lines that form each hexagon in the same ordering
hexcurveloops = []
hexsurfs = []
nodecount = 1 + 5*len(cyl_radii)*len(tubes) #start the node count after the last tube node
linecount = 1 + 4*len(cyl_radii)*len(tubes)
for i in range(len(hexes)):
    temp_hexnodes = []
    for j,vert in enumerate(hexes[i]):
        vert = [round(vert[0],4),round(vert[1],4)] #round off the values to make the dict search work
        #Save the vertex to the list of unique nodes
        #Create the point element based on the unique node
        if vert not in nodes.values():
            #----
            gmsh.model.geo.addPoint(vert[0],vert[1],0, lc, nodecount) #Since this is a new unique node, add the point element
            #----
            nodes[nodecount] = vert
            temp_hexnodes.append(nodecount)
            nodecount = nodecount + 1
        else:
            #print("Duplicate node found")
            temp_hexnodes.append( list(nodes.keys())[list(nodes.values()).index(vert)] ) #Get the index of the vert in the unique hexnodes and append it to this hex's nodes
    hexnodes.append(temp_hexnodes)

    #Now make the line elements that constitute each hexagon
    temp_hexlines = []
    for k in range(len(hexnodes[i])): #step through each list of nodes in each hexagon. The item i will have just been added by the previous loop
        if k < (len(hexnodes[i]) - 1):
            line = [ hexnodes[i][k],hexnodes[i][k+1] ] #By construction, the lines are formed with the conventional orientation
        else:
            line = [ hexnodes[i][k],hexnodes[i][0] ]#edge case for the last line
        if line not in lines.values():
            if line[::-1] not in lines.values(): #need to also make sure it's not unique just because of opposite orientation
                #-----
                gmsh.model.geo.addLine(line[0],line[1], linecount)
                #-----
                lines[linecount] = line
                temp_hexlines.append(linecount) #save the line numbers in the usual ordering
                linecount = linecount + 1
            else: #if the reversed line already exists
                #print("Duplicate line with opposite orientation found")
                temp_hexlines.append( -1*(list(lines.keys())[list(lines.values()).index(line[::-1])]) ) #Add the negative sign to denote the opposite orientation
        else:
            #print("Duplicate line found")
            temp_hexlines.append( list(lines.keys())[list(lines.values()).index(line)] )
    hexlines.append(temp_hexlines)
    #Go ahead and make the curveloop and planesurface
    #print("hex lines and tube lines that make the hex channel surf")
    #print(hexlines[i]+tubelines[i])
    temp_hexcurveloop = gmsh.model.geo.addCurveLoop(hexlines[i] + tubelines[i][-1], 4*len(cyl_radii)*len(tubes) + i+1) #note this approach of combining all the lines in one loop generates a warning, but it hasn't caused issues with the hex assembly model
    #print("index of the temp hex curve loop:")
    #print(len(tubelines)+i+1)
    hexcurveloops.append(temp_hexcurveloop)
    temp_hexsurf = gmsh.model.geo.addPlaneSurface([4*len(cyl_radii)*len(tubes) + i+1], 4*len(cyl_radii)*len(tubes) + i+1)
    hexsurfs.append(temp_hexsurf)


#print("Unique nodes forming hexes:")
#print(nodes)
#print("List of nodes making up each hex: ")
#print(hexnodes)
print("list of lines making up each hex:")
print(hexlines)

#---Reflector-related: pick out the hex lines that form the exterior of the core lattice. 

#for now just pick out the lower part
hexlines_lower = []
for i in range(N_hexes_per_row[0]):
    hexlines_lower.append(hexlines[i][5]) #save the last and first lines of each hex which make up the lower part. 
    hexlines_lower.append(hexlines[i][0]) 


hexlines_right = []
for i in range(len(N_hexes_per_row)):
    if i < (len(N_hexes_per_row)-1)/2.0: #if lower rows
        if i == 0:
            hexlines_right.append(hexlines[sum(N_hexes_per_row[:i]) + N_hexes_per_row[i] - 1][1] )
        else:
            hexlines_right.append(hexlines[sum(N_hexes_per_row[:i]) + N_hexes_per_row[i] - 1][0] )
            hexlines_right.append(hexlines[sum(N_hexes_per_row[:i]) + N_hexes_per_row[i] - 1][1] )
    elif i == (len(N_hexes_per_row)-1)/2.0: #if middle row
        hexlines_right.append(hexlines[sum(N_hexes_per_row[:i]) + N_hexes_per_row[i] - 1][0] ) 
        hexlines_right.append(hexlines[sum(N_hexes_per_row[:i]) + N_hexes_per_row[i] - 1][1] ) 
        hexlines_right.append(hexlines[sum(N_hexes_per_row[:i]) + N_hexes_per_row[i] - 1][2] ) 
    elif i > (len(N_hexes_per_row)-1)/2.0: #if upper rows
        if i == (len(N_hexes_per_row) - 1):
            hexlines_right.append(hexlines[sum(N_hexes_per_row[:i]) + N_hexes_per_row[i] - 1][1] )
        else:            
            hexlines_right.append(hexlines[sum(N_hexes_per_row[:i]) + N_hexes_per_row[i] - 1][1] )
            hexlines_right.append(hexlines[sum(N_hexes_per_row[:i]) + N_hexes_per_row[i] - 1][2] )

hexlines_upper = []
#for i in range (len(hexlines) - 1 - N_hexes_per_row[-1], len(hexlines)):
for i in range(sum(N_hexes_per_row[:-1]), sum(N_hexes_per_row)):
    hexlines_upper.append(hexlines[i][3])
    hexlines_upper.append(hexlines[i][2])
hexlines_upper.reverse() #occurs in-place

hexlines_left = []
for i in range(len(N_hexes_per_row)): 
    if i > (len(N_hexes_per_row)-1)/2.0: #if upper rows
        if i == (len(N_hexes_per_row) - 1):
            hexlines_left.append(hexlines[sum(N_hexes_per_row[:i])][4] )
        else:
            hexlines_left.append(hexlines[sum(N_hexes_per_row[:i])][4] )
            hexlines_left.append(hexlines[sum(N_hexes_per_row[:i])][3] )
    elif i == (len(N_hexes_per_row)-1)/2.0: # if middle row
        hexlines_left.append(hexlines[sum(N_hexes_per_row[:i])][5] )
        hexlines_left.append(hexlines[sum(N_hexes_per_row[:i])][4] )
        hexlines_left.append(hexlines[sum(N_hexes_per_row[:i])][3] )
        
    elif i < (len(N_hexes_per_row)-1)/2.0: #if lower rows
        if i == 0:
            hexlines_left.append(hexlines[sum(N_hexes_per_row[:i])][4] )
        else:
            hexlines_left.append(hexlines[sum(N_hexes_per_row[:i])][5] )
            hexlines_left.append(hexlines[sum(N_hexes_per_row[:i])][4] )
hexlines_left.reverse()

lattice_border_lines = hexlines_lower + hexlines_right + hexlines_upper + hexlines_left
#lattice_border_lines.remove(1922) 
print("----")
print("lines that make up the lattice border:")
print(lattice_border_lines)
print("----")
#Now define the reflector outer boundary
parc1 = gmsh.model.geo.addPoint(0, -1*reflector_radius, 0, lc*4) 
parc2 = gmsh.model.geo.addPoint(reflector_radius, 0, 0, lc*4) 
parc3 = gmsh.model.geo.addPoint(0, reflector_radius, 0, lc*4) 
parc4 = gmsh.model.geo.addPoint(-1*reflector_radius, 0, 0, lc*4) 

#Remember each element in "tubes" corresponds to one of the channels. 
#Each element contains a list of lists of 5 points, one for each radial layer. 
#Then within each sublist, the first element is the centrepoint of that circle/annulus. The index of that point was saved as a 3rd 'coordinate'.
centrepoint = tubes[len(tubes)//2][0][0][-1] 
arc1 = gmsh.model.geo.addCircleArc(parc1, centrepoint, parc2) #start, centre, end
arc2 = gmsh.model.geo.addCircleArc(parc2, centrepoint, parc3)
arc3 = gmsh.model.geo.addCircleArc(parc3, centrepoint, parc4)
arc4 = gmsh.model.geo.addCircleArc(parc4, centrepoint, parc1)

reflector_outerboundary = [arc1, arc2, arc3, arc4]


#---Tried the 'proper way' above but it led to annoying duplication of surfaces. This way is simpler. 
reflector_curveloop = gmsh.model.geo.addCurveLoop(reflector_outerboundary + lattice_border_lines)
reflector_surface = gmsh.model.geo.addPlaneSurface([reflector_curveloop]) 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#----Now the physical surfs
ptags_tubes = []
ptags_hexes = []
for i in range(len(tubesurfs)):
    temp_ptag = gmsh.model.addPhysicalGroup(2,[tubesurfs[i]],i+1)
    ptags_tubes.append(temp_ptag)
    gmsh.model.setPhysicalName(2,ptags_tubes[-1], name = str(i+1))
for i in range(len(hexsurfs)):
    temp_ptag = gmsh.model.addPhysicalGroup(2,[hexsurfs[i]], len(tubesurfs) + i+1)
    ptags_hexes.append(temp_ptag)
    gmsh.model.setPhysicalName(2, ptags_hexes[-1], name = str(len(tubesurfs) + i+1))
    
    
ptag_reflector = gmsh.model.addPhysicalGroup(1, reflector_outerboundary, 1) 
gmsh.model.setPhysicalName(1, ptag_reflector, name="outer")

# Reflector
ptag_reflectorbase = gmsh.model.addPhysicalGroup(2, [reflector_surface], len(tubesurfs) + len(hexsurfs) + 1 )
gmsh.model.setPhysicalName(2, ptag_reflectorbase, name="reflector_base")

# #@@@@@@@@@@@@@@@@@@@@@@@@@@@
gmsh.model.geo.synchronize() 
# #@@@@@@@@@@@@@@@@@@@@@@@@@@@


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#gmsh.option.setNumber('Mesh.CharacteristicLengthMax', lc*4)
#gmsh.option.setNumber('Mesh.CharacteristicLengthMin', 0.008)

#print(gmsh.option.getString("General.BuildOptions"))
#gmsh.option.setNumber("Mesh.Algorithm", 2)
gmsh.option.setNumber("Mesh.Algorithm", 2) #option 3 no longer works after upgrading to Ubuntu 22.04... Option 2 is automatic algorithm selection
gmsh.option.setNumber("Mesh.Format", 2)
# gmsh.option.setNumber("Mesh.RecombineAll", 1)
#optional additional paramters for FEMFFUSION
# gmsh.option.setNumber("Mesh.Algorithm3D", 1) #1=Delauney
# gmsh.option.setNumber("Mesh.Recombine3DAll", 1) #1=Delauney
#  #recombine triangle elements to form quadrangles. This is important because Deal.II does not support triangles or tetrahedra.


# gmsh.option.setNumber("Mesh.Algorithm", 8) #option 3 no longer works after upgrading to Ubuntu 22.04... Option 2 is automatic algorithm selection
# gmsh.option.setNumber("Mesh.Format", 2)
# gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2);
# # gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
# gmsh.option.setNumber("Mesh.RecombineAll", 1)

print(gmsh.option.getNumber("Mesh.Algorithm"))
#%%%%%%%%%%%%%%%%%
print("Generating mesh...")
gmsh.model.mesh.generate(2) #produce 3D mesh automatically
print("Done!")
gmsh.option.setNumber('Mesh.MshFileVersion', 2.0)
gmsh.write(filename + ".msh") #might need some special formatting tricks, ASCII etc to get gmshToFoam to work properly


#%%%%%%%%%%%%%%%%%
# finish up
log = gmsh.logger.get()
print("Logger has recorded " + str(len(log)) + " lines")
print(log)

gmsh.logger.stop()

if '-nopopup' not in sys.argv: #open the gmsh GUI unless specified otherwise
    gmsh.fltk.run()


gmsh.finalize()
