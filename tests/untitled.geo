//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 0.5, 0, 2*Pi};
//+
Dilate {{0, 0, 0}, {0.5, 0.5, 0.5}} {
  Duplicata { Curve{1}; }
}
//+
Curve Loop(1) = {2};
//+
Curve Loop(2) = {1};

//+
Curve Loop(3) = {2};
//+
Plane Surface(1) = {3};
//+
Curve Loop(4) = {1};
//+
Curve Loop(5) = {2};
//+
Plane Surface(2) = {4, 5};

//+
Recombine Surface {2};
