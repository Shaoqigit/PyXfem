//+
Point(1) = {-0.5, 0, 0., 1.0};
Point(2) = {0., 0, 0., 1.0};
Point(3) = {0.5, 0, 0., 1.0};
//+
Point(4) = {0.5, 0.1, 0., 1.0};
Point(5) = {0., 0.1, 0., 1.0};
//+
Point(6) = {-0.5, 0.1, 0., 1.0};
//+
Point(7) = {-0.5, 0., 0.1, 1.0};
Point(8) = {-0., 0., 0.1, 1.0};
//+
Point(9) = {0.5, 0., 0.1, 1.0};
//+
Point(10) = {0.5, 0.1, 0.1, 1.0};
Point(11) = {0., 0.1, 0.1, 1.0};
//+
Point(12) = {-0.5, 0.1, 0.1, 1.0};
//+


//+
Line(1) = {7, 8};
//+
Line(2) = {8, 9};
//+
Line(3) = {9, 10};
//+
Line(4) = {10, 11};
//+
Line(5) = {11, 12};
//+
Line(6) = {12, 7};
//+
Line(7) = {11, 8};
//+
Line(8) = {12, 6};
//+
Line(9) = {6, 1};
//+
Line(10) = {1, 7};
//+
Line(11) = {1, 2};
//+
Line(12) = {2, 3};
//+
Line(13) = {3, 4};
//+
Line(14) = {4, 10};
//+
Line(15) = {3, 9};
//+
Line(16) = {4, 5};
//+
Line(17) = {5, 2};
//+
Line(18) = {5, 6};
//+
Line(19) = {2, 8};
//+
Line(20) = {5, 11};
//+
Curve Loop(1) = {9, 10, -6, 8};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 1, -7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, -18, 20, 5};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {18, 9, 11, -17};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {11, 19, -1, -10};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {16, 17, 12, 13};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {16, 20, -4, -14};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {4, 7, 2, 3};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {2, -15, -12, 19};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {15, 3, -14, -13};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {19, -7, -20, 17};
//+
Plane Surface(11) = {11};
//+
Surface Loop(1) = {5, 4, 3, 1, 2, 11};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {9, 8, 7, 6, 10, 11};
//+
Volume(2) = {2};
//+
Physical Surface("inlet", 21) = {1};
//+
Physical Volume("air", 22) = {1};
//+
Physical Volume("foam", 23) = {2};
//+
//+
Transfinite Curve {11, 1, 18, 5, 12, 2, 16, 4} = 10 Using Progression 1;
