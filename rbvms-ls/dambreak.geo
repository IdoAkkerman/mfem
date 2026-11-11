// Parameter
x = 1.0;
y = 1.0;


// Pointss
p1=newp; Point(p1) = {0.0, 0.0, 0.0, 1.0};
p2=newp; Point(p2) = {  x, 0.0, 0.0, 1.0};
p3=newp; Point(p3) = {  x,   y, 0.0, 1.0};
p4=newp; Point(p4) = {0.0,   y, 0.0, 1.0};

// Lines
l11=newl; Line(l11) = {p1, p2};
l12=newl; Line(l12) = {p2, p3};
l13=newl; Line(l13) = {p3, p4};
l14=newl; Line(l14) = {p4, p1};

cl1=newcl; Curve Loop(cl1) = {l11, l12, l13, l14};
s1=news; Plane Surface(s1)={cl1};

// Physical
Physical Line(1) = {l11, l12, l13, l14};
Physical Surface(100) = {s1};
