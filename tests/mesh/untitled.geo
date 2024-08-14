//+
Point(1) = {0., 0, 0., 1.0};
Point(2) = {1., 0, 0., 1.0};
//+
Point(3) = {1, 0.2, 0., 1.0};
//+
Point(4) = {0., 0.2, 0., 1.0};
//+
Point(5) = {0., 0., 0.2, 1.0};
//+
Point(6) = {1, 0., 0.2, 1.0};
//+
Point(7) = {1, 0.2, 0.2, 1.0};
//+
Point(8) = {0., 0.2, 0.2, 1.0};
//+
Line(1) = {5, 6};
//+
Line(2) = {6, 7};
//+
Line(3) = {7, 8};
//+
Line(4) = {8, 5};
//+
Line(5) = {5, 1};
//+
Line(6) = {1, 4};
//+
Line(7) = {4, 3};
//+
Line(8) = {3, 2};
//+
Line(9) = {2, 1};
//+
Line(10) = {3, 7};
//+
Line(11) = {2, 6};
//+
Line(12) = {4, 8};
//+
Curve Loop(1) = {12, 4, 5, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, 4, 1, 2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, -2, -11, -8};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {7, 8, 9, 6};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {1, -11, 9, -5};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {3, -12, 7, 10};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {6, 2, 1, 5, 3, 4};
//+
Volume(1) = {1};
//+
Physical Volume("air", 13) = {1};
//+
Physical Surface("inlet", 14) = {3};
//+
Transfinite Curve {3, 1, 7, 9} = 5 Using Progression 1;
