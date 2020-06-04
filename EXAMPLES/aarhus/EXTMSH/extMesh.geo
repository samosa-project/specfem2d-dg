// Display. //---------------//
General.Axes = 1; // Show axes.
Mesh.Lines = 1; // Display lines.
Geometry.LineNumbers = 1; // Display line numbers.
Geometry.PointNumbers = 1; // Display point numbers
Geometry.SurfaceNumbers = 1; // Display surface numbers.
Mesh.Points = 1; // Display points.
Mesh.PointNumbers = 0; // Display point numbers
Mesh.LineNumbers = 0; // Display line numbers.
Mesh.SurfaceEdges = 1; // Display edges.
Mesh.SurfaceFaces = 1; // Display surfaces.
Mesh.SurfaceNumbers = 0; // Display surface numbers.

// Elements. //--------------//
Mesh.SubdivisionAlgorithm = 1; // Meshing algorithm: all quads.
Mesh.RecombineAll = 1; // Recombine all triangular elements.

// Meshing algorithm. //-----//
Mesh.Algorithm = 9; // 1: MeshAdapt. 5: Delaunay for quads. 9: structured (experimental).
Mesh.ElementOrder = 1; // Element order.
Mesh.RandomFactor=1e-3; // Perturbate every points' positions in order to avoid 3 aligned points. Default is at 1e-9, maximum is 1e-3.

// Geometry. //--------------//
// Define elements' length.
dxair = 35.e-3; // lambda_min = c_min/f_max = 250/4500 = 55e-3 => dx = lambda_min/1.5 = 35e-3

tunnel_length = 10.;
tunnel_diamet = 2.1;
tunnel_separa = 1.;
speaker_x = 0.;
mic1 = 0.3;
mic2 = 0.5;
mic3 = 1.;
mic4 = 2.;
mic5 = 3.;
thicc = 2.*dxair;
xmin = -1.;
//xmax = tunnel_length
xmax = mic5+1.;
zmin = -0.5*tunnel_separa;
zmax = +0.5*tunnel_separa;

dxtop=20;

// Define geometry keypoints.
Point(1) = {xmin,zmin,0};
Point(2) = {xmax,zmin,0};
Point(3) = {xmax,zmax,0};
Point(4) = {xmin,zmax,0};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
//Line Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};

Point(1011) = {xmin-thicc,zmin,0};
Point(1012) = {xmin-thicc,zmin-thicc,0};
Point(1013) = {xmin,zmin-thicc,0};

Point(2011) = {xmax,zmin-thicc,0};
Point(2012) = {xmax+thicc,zmin-thicc,0};
Point(2013) = {xmax+thicc,zmin,0};

Point(3011) = {xmax+thicc,zmax,0};
Point(3012) = {xmax+thicc,zmax+thicc,0};
Point(3013) = {xmax,zmax+thicc,0};

Point(4011) = {xmin,zmax+thicc,0};
Point(4012) = {xmin-thicc,zmax+thicc,0};
Point(4013) = {xmin-thicc,zmax,0};

Line(72) = {1013,2011};
Line(73) = {2013,3011};
Line(74) = {3013,4011};
Line(75) = {4013,1011};

Line(77) = {2011, 2};
Line(78) = {1, 1013};
Line(79) = {3, 3013};
Line(80) = {4011, 4};

Line(82) = {1, 1011};
Line(83) = {1011, 1012};
Line(84) = {1012, 1013};
Line(85) = {4, 4013};
Line(86) = {4011, 4012};
Line(87) = {4012, 4013};
Line(88) = {2011, 2012};
Line(89) = {2012, 2013};
Line(90) = {2013, 2};
Line(91) = {3, 3011};
Line(98) = {3011, 3012};
Line(99) = {3012, 3013};

nelts = (xmax-xmin)/dxair;
Transfinite Line{1, 3, 72, 74} = nelts/2+1;
nelts = (zmax-zmin)/dxair;
Transfinite Line{2, 4, 73, 75} = nelts/2+1;
//nelts = (thicc)/dxair;
Transfinite Line{78, 77, 79, 80, 82, 83, 84, 90, 88, 89, 91, 98, 99, 85, 86, 87} = 1;

i = 1;
Line Loop(i) = {1,2,3,4}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {72, 77, -1, 78}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-2, -90, 73, -91}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-3, 79, 74, 80}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-4, 85, 75, -82}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-78, 82, 83, 84}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {88, 89, 90, -77}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {91, 98, 99, -79}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-80, 86, 87, -85}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;

// Set element size in the neighbourhood of some points.
//Mesh.CharacteristicLengthFactor = 1; // Scaling of all elements' sizes.
//Characteristic Length {1, 2} = 2*dxbot;
//Characteristic Length {4, 3} = 2*dxtop;

// ******************************** //
// Finally, "name" lines and surfaces in order to be able to transfer this to SPECFEM standards.
// Use these whatever was done before.
// ******************************** //
Physical Line("Top") = {99, 74, 86}; // Associate line 4 (RHS) to "Top" (LHS) SPECFEM boundary. Warning: if many lines, list them in counter-clockwise order.
Physical Line("Left") = {87, 75, 83}; // Associate line 4 (RHS) to "Left" (LHS) SPECFEM boundary. Warning: if many lines, list them in counter-clockwise order.
Physical Line("Bottom") = {84, 72, 88}; // Associate line 4 (RHS) to "Bottom" (LHS) SPECFEM boundary. Warning: if many lines, list them in counter-clockwise order.
Physical Line("Right") = {89, 73, 98}; // Associate line 4 (RHS) to "Right" (LHS) SPECFEM boundary. Warning: if many lines, list them in counter-clockwise order.
Physical Surface("M1") = {1}; // Associate surfaces (RHS) to SPECFEM material (LHS).
Physical Surface("M2") = {2, 3, 4, 5, 6, 7, 8, 9};
Recombine Surface {1, 2, 3, 4, 5, 6, 7, 8, 9};
