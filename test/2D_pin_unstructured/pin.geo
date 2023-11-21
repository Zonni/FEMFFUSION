cm = 1;
lc = cm;
pin_size = 2* cm;
fuel_radius = 0.7* cm;

// Square 
Point(1) = {0, 0, 0, lc};
Point(2) = {pin_size, 0, 0, lc};
Point(3) = {0, pin_size, 0, lc};
Point(4) = {pin_size, pin_size, 0, lc};
Point(5) = {pin_size/2, pin_size/2, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};


//Inner Fuel Circle
Point(6) = {pin_size/2, pin_size/2 - fuel_radius, 0, lc};
Point(7) = {pin_size/2 + fuel_radius, pin_size/2, 0, lc};
Point(8) = {pin_size/2, pin_size/2 + fuel_radius, 0, lc};
Point(9) = {pin_size/2 - fuel_radius, pin_size/2, 0, lc};

Circle(6) = {6, 5, 7};
Circle(7) = {7, 5, 8};
Circle(8) = {8, 5, 9};
Circle(9) = {9, 5, 6};

// Moderator
Line Loop(11) = {1, 2, 3, 4, -6, -7, -8, -9};
Plane Surface(12) = {11};
Physical Surface("Moderator", 2) = {12};
Physical Line(0) = {4};
Physical Line(1) = {2};
Physical Line(2) = {1};
Physical Line(3) = {3};


// Fuel
Line Loop(12) = {6, 7, 8, 9};
Plane Surface(13) = {12};
Physical Surface("Fuel", 1) = {13};

