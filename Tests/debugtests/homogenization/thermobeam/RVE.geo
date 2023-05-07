// Gmsh project created on Mon Feb 27 21:30:30 2023
SetFactory("OpenCASCADE");
//+
L = DefineNumber[ 2, Name "Parameters/L" ];
//+
h = DefineNumber[ 1, Name "Parameters/h" ];
//+
b = DefineNumber[ 1, Name "Parameters/b" ];
//+
nx = DefineNumber[ 6, Name "Parameters/nx" ];
//+
ny = DefineNumber[ 6, Name "Parameters/ny" ];
//+
nz = DefineNumber[ 6, Name "Parameters/nz" ];
//+
Point(1) = {-L/2, -b/2, -h/2, 1.0};
//+
Point(2) = {-L/2, b/2, -h/2, 1.0};
//+
Point(3) = {-L/2, b/2, h/2, 1.0};
//+
Point(4) = {-L/2, -b/2, h/2, 1.0};
//+
Point(5) = {0, -b/2, -h/2, 1.0};//+
Point(6) = {0, b/2, -h/2, 1.0};
//+
Point(7) = {0, b/2, h/2, 1.0};
//+
Point(8) = {0, -b/2, h/2, 1.0};
//+
Point(9) = {L/2, -b/2, -h/2, 1.0};
//+
Point(10) = {L/2, b/2, -h/2, 1.0};
//+
Point(11) = {L/2, b/2, h/2, 1.0};
//+
Point(12) = {L/2, -b/2, h/2, 1.0};
//+
Point(13) = {0, 0, 0, 1.0};
//+
Point(14) = {0, 0, 0, 1.0};
//+
Line(1) = {4, 8};
//+
Line(2) = {8, 5};
//+
Line(3) = {5, 1};
//+
Line(4) = {1, 4};
//+
Line(5) = {4, 3};
//+
Line(6) = {3, 2};
//+
Line(7) = {2, 1};
//+
Line(8) = {2, 6};
//+
Line(9) = {6, 10};
//+
Line(10) = {10, 11};
//+
Line(11) = {11, 7};
//+
Line(12) = {7, 6};
//+
Line(13) = {7, 3};
//+
Line(14) = {7, 8};
//+
Line(15) = {5, 6};
//+
Line(16) = {10, 9};
//+
Line(17) = {9, 5};
//+
Line(18) = {9, 12};
//+
Line(19) = {12, 11};
//+
Line(20) = {12, 8};
//+
Curve Loop(1) = {5, 6, 7, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {8, -15, 3, -7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {12, -8, -6, -13};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {3, 4, 1, 2};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {1, -14, 13, -5};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {14, 2, 15, -12};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {20, 2, -17, 18};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {17, 15, 9, 16};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {12, 9, 10, 11};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {11, 14, -20, 19};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {18, 19, -10, 16};
//+
Plane Surface(11) = {11};
//+
Transfinite Curve {3, 1, 13, 8, 11, 20, 9, 17} = nx Using Progression 1;
//+
Transfinite Curve {4, 6, 2, 12, 18, 10} = nz Using Progression 1;
//+
Transfinite Curve {5, 7, 14, 15, 19, 16} = ny Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {4};
//+
Transfinite Surface {5};
//+
Transfinite Surface {3};
//+
Transfinite Surface {2};
//+
Transfinite Surface {6};
//+
Transfinite Surface {10};
//+
Transfinite Surface {9};
//+
Transfinite Surface {8};
//+
Transfinite Surface {7};
//+
Transfinite Surface {11};
//+
Surface Loop(1) = {5, 4, 2, 3, 1, 6};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {10, 9, 8, 7, 11, 6};
//+
Volume(2) = {2};
//+
Transfinite Volume{1};
//+
Transfinite Volume{2};
//+
Physical Volume("volelements", 21) = {1, 2};
//+
Physical Surface("interface", 22) = {6};
//+
Physical Point("interfacePoint", 23) = {13};
//+
Physical Point("constraintPoint", 24) = {14};
