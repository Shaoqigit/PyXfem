//+
Point(1) = {-0.5, 0, 0, 1.0};
//+
Point(2) = {0, 0, 0, 1.0};
//+
Point(3) = {0.5, 0, 0, 1.0};
//+
Point(4) = {0.5, 1, 0, 1.0};
//+
Point(5) = {0., 1, 0, 1.0};
//+
Point(6) = {-0.5, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Line(7) = {5, 2};
//+
Curve Loop(1) = {6, 1, -7, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, 4, 7, 2};
//+
Plane Surface(2) = {2};
//+
Physical Surface("mat1", 8) = {1};
//+
Physical Surface("mat2", 9) = {2};
//+
Physical Curve("left", 10) = {6};
//+
Physical Curve("left_bot", 11) = {1};
//+
Physical Curve("left_top", 12) = {5};
//+
Physical Curve("right_top", 13) = {4};
//+
Physical Curve("right_bot", 14) = {2};
//+
Physical Curve("right", 15) = {3};
