run mat;


M = [A_11 -A_12  A_13 -A_14 ;
     A_21  A_22  A_23  A_24 ;
     A_31 -A_32  A_33 -A_34 ;
     A_41  A_42  A_43  A_44 ;
];

norm(rhs,2)
x = M\rhs;

norm(x, 2)

W = A_44;
T = A_34;


mat = -W-T; 
issymmetric(mat)
eig(mat)
max(eig(mat))


mat = T-W; 
issymmetric(mat)
eig(mat)
max(eig(mat))

