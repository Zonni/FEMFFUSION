// Gmsh project created on Thu Jun 27 10:28:22 2013
lc = 100 ;
cm =1 ;    // 0.01 
//   cell  1 
l =13.62546635* cm ;
ll = 40.87639905 * cm ; // ll = l*3
l5 =20.438199525 ;
r =11.8 * cm ;
rr =23.6   * cm ;  // rr = r*2
rr1 =23.6   * cm ;  // rr = r*2
rr2 =23.6   * cm ;  // rr = r*2
rr3 =23.6   * cm ;  // rr = r*2
rr4 =23.6   * cm ;  // rr = r*2

s = 6.812733175* cm ;  // R- 0.5 l = 0.5 *l
ls = 20.438199525 * cm ;   // ll = l+s
Point(1) = {0, 0, 0, lc} ;
Point(2) = {r, -s, 0, lc} ;
Point(3) = {rr, 0, 0, lc} ;
Point(4) = {rr, l, 0, lc} ;
Point(5) = {r, ls, 0, lc} ;
Point(6) = {0, l, 0, lc} ;

/////

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};


Line Loop(7) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {7};

For t In {1:4}
  rr += 23.6  ;
my_new_surfs[] = Translate {rr-23.6 , 0, 0} { Duplicata{ Surface{1}; }  };
EndFor
 For t In {1:4} //  13/11/9 column
  rr1 +=23.6  ;
my_new_surfs[] = Translate {rr1-r*3,l5, 0} { Duplicata{ Surface{1}; }  };
EndFor
For t In {1:4} //  13/11/9 column
  rr2 +=23.6  ;
my_new_surfs[] = Translate {rr2-r*3,-l5, 0} { Duplicata{ Surface{1}; }  };
EndFor
//

For t In {1:3} //  13/11/9 column
  rr3 +=23.6  ;
my_new_surfs[] = Translate {rr3-r*2,l5*2, 0} { Duplicata{ Surface{1}; }  };
EndFor
For t In {1:3} //  13/11/9 column
  rr4 +=23.6  ;
my_new_surfs[] = Translate {rr4-r*2,-l5*2, 0} { Duplicata{ Surface{1}; }  };
EndFor

Physical Surface(1) = {76};
Physical Surface(2) = {83};
Physical Surface(3) = {89};
Physical Surface(4) = {32};
Physical Surface(5) = {39};
Physical Surface(6) = {45};
Physical Surface(7) = {51};
Physical Surface(8) = {1};
Physical Surface(9) = {8};
Physical Surface(10) = {14};
Physical Surface(11) = {20};
Physical Surface(12) = {26};
Physical Surface(13) = {57};
Physical Surface(14) = {64};
Physical Surface(15) = {68};
Physical Surface(16) = {72};
Physical Surface(17) = {95};
Physical Surface(18) = {102};
Physical Surface(19) = {106};
Physical Line(1) = {6, 1, 63, 58, 101, 96, 97, 103, 104, 107, 108, 109, 74, 75, 28, 29, 30, 54, 55, 92, 93, 94, 87, 88, 80, 81, 82, 37, 38, 5};
