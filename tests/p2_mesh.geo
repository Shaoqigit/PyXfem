// Gmsh project created on Tue Sep 24 22:47:00 2024
//+
Point(1) = {-0., 0., 0, 1.0};
//+
Point(2) = {1.0, 0., 0, 1.0};
//+
Point(3) = {0., 1.0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 1};
//+
Curve Loop(1) = {3, 1, 2};
//+
Plane Surface(1) = {1};
