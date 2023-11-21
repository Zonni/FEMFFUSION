cm = 1;
pin_size = 0.4 * cm;
radius = 0.15 * cm;
lc = 3*cm;
n_pinsx = 1;
n_pinsy = 1;


pitch = 15 * cm;
ri = pitch / 2;        // Radius Incircle  
rc = 2.0/3.0*Sqrt(3) * ri; // Radio Excircles == Side

// Make an hexagon
Macro Hexagon

	cx =  0.0;  // Center position x
	cy =  0.0;  // Center position y

	p0 = newp; Point(p0) = {cx+0.0,  cy+0.0,  0, lc}; // Center
	p1 = newp; Point(p1) = {cx+0.0,   cy-rc,  0, lc}; // Bottom
	p2 = newp; Point(p2) = {cx+ri,  cy-rc/2,  0, lc}; // Bottom Right
	p3 = newp; Point(p3) = {cx+ri,  cy+rc/2,  0, lc}; // Top Right
	p4 = newp; Point(p4) = {cx+0.0,   cy+rc,  0, lc}; // Top 
	p5 = newp; Point(p5) = {cx-ri,  cy+rc/2,  0, lc}; // Top Left
	p6 = newp; Point(p6) = {cx-ri,  cy-rc/2,  0, lc}; // Bottom Left


	l1 = newc; Line(l1) = {p1, p2}; // Bottom Right
	l2 = newc; Line(l2) = {p2, p3}; // Right
	l3 = newc; Line(l3) = {p3, p4}; // Top Right
	l4 = newc; Line(l4) = {p4, p5}; // Top Left
	l5 = newc; Line(l5) = {p5, p6}; // Left
	l6 = newc; Line(l6) = {p6, p1}; // Left Right

	ll = newreg; Line Loop(ll) = {l1, l2, l3, l4, l5, l6};


	Call CircleHole;

	hexagon = newreg;
	Plane Surface(hexagon) = ll;
	Physical Surface(id) = {hexagon};
	Printf("Hexagon %g has number %g, Line loop=%g", id, hexagon, ll);
	id = id + 1;

Return

// Make a cicle hole 
Macro CircleHole

	cx =  0.0;  // Center position x
	cy =  0.0;  // Center position y
	//radius  

	//Inner Fuel Circle
	p5 = newp; Point(p5) = {cx, cy, 0, lc};          // Center
	p6 = newp; Point(p6) = {cx, cy - radius, 0, lc}; // Bottom
	p7 = newp; Point(p7) = {cx + radius, cy, 0, lc}; // Right
	p8 = newp; Point(p8) = {cx, cy + radius, 0, lc}; // Top	
	p9 = newp; Point(p9) = {cx - radius, cy, 0, lc}; // Left

	c1 = newreg; Circle(c1) = {p6, p5, p7};
	c2 = newreg; Circle(c2) = {p7, p5, p8};
	c3 = newreg; Circle(c3) = {p8, p5, p9};
	c4 = newreg; Circle(c4) = {p9, p5, p6};

	ll = newreg; Line Loop(ll) = {c1, c2, c3, c4};
	//id = pin_x + (n_pinsx * pin_y) + 1;
	//theloops[id] = ll;
	 
	circle = newreg;
	Plane Surface(circle) = ll;

	Physical Surface(id) = {circle};
	Printf("Circle %g has number %g, Line loop=%g", id, circle, ll) ;
	id = id + 1;


Return

// 
id = 1;

// Fuel Circles
y = 0;
For pin_y In {0:n_pinsy-1}
	x = 0;
	For pin_x In {0:n_pinsx-1}
		Call Hexagon;

		//Call CircleHole;

		x += pin_size;
		EndFor
		y += pin_size;
EndFor


//Physical Line(0) = {4};
//Physical Line(1) = {2};
//Physical Line(2) = {1};
//Physical Line(3) = {3};


Mesh.Algorithm = 5;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 3;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 10;
