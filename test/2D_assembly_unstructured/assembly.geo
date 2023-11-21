cm = 1;
pin_size = 0.4 * cm;
fuel_radius = 0.15 * cm;
lc = 0.2*cm;
n_pinsx = 5;
n_pinsy = 5;

assembly_sizex = n_pinsx * pin_size;
assembly_sizey = n_pinsy * pin_size;


// Assembly Contour
// Square 
Point(1) = {0, 0, 0, lc};
Point(2) = {assembly_sizex, 0, 0, lc};
Point(3) = {0, assembly_sizey, 0, lc};
Point(4) = {assembly_sizex, assembly_sizey, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

// 
Macro CircleHole
	//Inner Fuel Circle
	p5 = newp; Point(p5) = {x+pin_size/2, y+pin_size/2, 0, lc};
	p6 = newp; Point(p6) = {x+pin_size/2, y+pin_size/2 - fuel_radius, 0, lc};
	p7 = newp; Point(p7) = {x+pin_size/2 + fuel_radius, y+pin_size/2, 0, lc};
	p8 = newp; Point(p8) = {x+pin_size/2, y+pin_size/2 + fuel_radius, 0, lc};
	p9 = newp; Point(p9) = {x+pin_size/2 - fuel_radius, y+pin_size/2, 0, lc};

	c1 = newreg; Circle(c1) = {p6, p5, p7};
	c2 = newreg; Circle(c2) = {p7, p5, p8};
	c3 = newreg; Circle(c3) = {p8, p5, p9};
	c4 = newreg; Circle(c4) = {p9, p5, p6};

	ll = newreg; Line Loop(ll) = {c1, c2, c3, c4};
	pin = pin_x + (n_pinsx * pin_y) + 1;
	theloops[pin] = ll;

	
	// We then store the line loops identification numbers in a list for later
    	// reference (we will need these to define the final volume):
	thehole = newreg ;
	Plane Surface(thehole) = ll;
	Physical Surface(pin) = {thehole};
	Printf("Pin %g has number %g, Line loop=%g", pin, thehole, ll) ;

Return



// Fuel Circles
y = 0;
For pin_y In {0:n_pinsy-1}
	x = 0;
	For pin_x In {0:n_pinsx-1}
		Call CircleHole;
		
		x += pin_size;
		EndFor
		y += pin_size;
EndFor

// Moderator
ll = newreg;
theloops[0] = ll;
Line Loop(ll) = {1, 2, 3, 4};
ps = newreg;
Plane Surface(ps) = theloops[];
Physical Surface(pin+1) = {ps};

Physical Line(0) = {4};
Physical Line(1) = {2};
Physical Line(2) = {1};
Physical Line(3) = {3};

