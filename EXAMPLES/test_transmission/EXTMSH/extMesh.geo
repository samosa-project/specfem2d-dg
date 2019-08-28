// Display. //---------------//
General.Axes = 0; // Show axes.
Mesh.Lines = 1; // Display lines.
Geometry.LineNumbers = 1; // Display line numbers.
Mesh.LineNumbers = 0; // Display line numbers.
Mesh.Points = 1; // Display points.
Geometry.PointNumbers = 1; // Display point numbers
Mesh.PointNumbers = 0; // Display point numbers
Mesh.SurfaceEdges = 1; // Display edges.
Mesh.SurfaceFaces = 1; // Display surfaces.
Geometry.SurfaceNumbers = 1; // Display surface numbers.
Mesh.SurfaceNumbers = 0; // Display surface numbers.

// Elements. //--------------//
Mesh.SubdivisionAlgorithm = 1; // Meshing algorithm: all quads.
Mesh.RecombineAll = 1; // Recombine all triangular elements.

// Meshing algorithm. //-----//
Mesh.Algorithm = 9; // 1: MeshAdapt. 5: Delaunay for quads. 9: structured (experimental).
Mesh.ElementOrder = 2; // Element order.
Mesh.RandomFactor=1e-3; // Perturbate every points' positions in order to avoid 3 aligned points. Default is at 1e-9, maximum is 1e-3.

// Geometry. //--------------//
xmin = -150.;
xmax = 75.;
//xmax = -xmin;

triangle_halfbase = 50.;

zmin = -75.;
zmax = 75.;

zinterface = 0.;
angle_rad = 30. *3.14159265358979323846/180.;
triangle_top = triangle_halfbase*Tan(angle_rad);
interface_triangle_top = zinterface + triangle_top;
shift_triangle_horizontal = +0.;

triangle_peak_xpos = 0.+shift_triangle_horizontal;

// Define elements' length.
dxGlob = 5.;

// Define geometry keypoints.
i = 1;
Point(i) = {xmin,zmin,0}; i=i+1;
Point(i) = {xmin+2.*dxGlob,zmin,0}; i=i+1;
Point(i) = {triangle_peak_xpos,zmin,0}; i=i+1;
Point(i) = {xmax-2.*dxGlob,zmin,0}; i=i+1;
Point(i) = {xmax,zmin,0}; i=i+1;

Point(i) = {xmax,0.5*(zmin+zinterface),0}; i=i+1; // cut solid right
Point(i) = {xmax-2.*dxGlob,0.5*(zmin+zinterface),0}; i=i+1; // cut solid right
Point(i) = {0.,0.5*(zmin+zinterface),0}; i=i+1; // under triangle
//Point(i) = {triangle_peak_xpos,zinterface-triangle_top,0}; i=i+1; // under triangle
Point(i) = {xmin+2.*dxGlob,0.5*(zmin+zinterface),0}; i=i+1; // cut solid left
Point(i) = {xmin,0.5*(zmin+zinterface),0}; i=i+1; // cut solid left

Point(i) = {xmin,zinterface,0}; i=i+1;
Point(i) = {xmin+2.*dxGlob,zinterface,0}; i=i+1;
Point(i) = {triangle_peak_xpos-triangle_halfbase,zinterface,0}; i=i+1; // triangle
Point(i) = {triangle_peak_xpos,interface_triangle_top,0}; i=i+1; // triangle
Point(i) = {triangle_peak_xpos+triangle_halfbase,zinterface,0}; i=i+1; // triangle
Point(i) = {xmax-2.*dxGlob,zinterface,0}; i=i+1;
Point(i) = {xmax,zinterface,0}; i=i+1;

//Point(i) = {xmax,0.5*(zmax+zinterface),0}; i=i+1; // cut air right
//Point(i) = {xmin,0.5*(zmax+zinterface),0}; i=i+1; // cut air left
Point(i) = {xmax,interface_triangle_top,0}; i=i+1; // cut air right
Point(i) = {xmax-2.*dxGlob,interface_triangle_top,0}; i=i+1; // cut air right
Point(i) = {xmin+2.*dxGlob,interface_triangle_top,0}; i=i+1; // cut air left
Point(i) = {xmin,interface_triangle_top,0}; i=i+1; // cut air left

Point(i) = {xmin,zmax,0}; i=i+1;
Point(i) = {xmin+2.*dxGlob,zmax,0}; i=i+1;
Point(i) = {triangle_peak_xpos,zmax,0}; i=i+1;
Point(i) = {xmax-2.*dxGlob,zmax,0}; i=i+1;
Point(i) = {xmax,zmax,0}; i=i+1;


// Define lines joining points. Warning: prefer counter-clockwise direction to prevent inverted elements. Warning: do not duplicate lines if two or more materials are needed.
i = 1;
Line(i) = {1,2}; i=i+1; // lower left left
Line(i) = {2,9}; i=i+1;
Line(i) = {9,10}; i=i+1;
Line(i) = {10,1}; i=i+1;

Line(i) = {2,3}; i=i+1; // lower left
Line(i) = {3,8}; i=i+1;
Line(i) = {8,9}; i=i+1;
//-2

Line(i) = {3,4}; i=i+1; // lower right
Line(i) = {4,7}; i=i+1;
Line(i) = {7,8}; i=i+1;
//-6

Line(i) = {4,5}; i=i+1; // lower right right
Line(i) = {5,6}; i=i+1;
Line(i) = {6,7}; i=i+1;
//-9

//-3 // low left left
Line(i) = {9,12}; i=i+1;
Line(i) = {12,11}; i=i+1;
Line(i) = {11,10}; i=i+1;
//-7 // low left
Line(i) = {8,13}; i=i+1;
Line(i) = {13,12}; i=i+1;
//-14

Line(i) = {8,15}; i=i+1; // triangle
Line(i) = {15,14}; i=i+1;
Line(i) = {14,13}; i=i+1;
//-17

//-10
Line(i) = {7,16}; i=i+1; // low right
Line(i) = {16,15}; i=i+1;
//-19

//-13
Line(i) = {6,17}; i=i+1; // low right right
Line(i) = {17,16}; i=i+1;
//-22

//-15
Line(i) = {12,20}; i=i+1;
Line(i) = {20,21}; i=i+1;
Line(i) = {21,11}; i=i+1;

//-18
//-21
Line(i) = {14,20}; i=i+1;
//-26

//-20
//-23
Line(i) = {16,19}; i=i+1;
Line(i) = {19,14}; i=i+1;

//-25
Line(i) = {17,18}; i=i+1;
Line(i) = {18,19}; i=i+1;
//-30

//-27
Line(i) = {20,23}; i=i+1;
Line(i) = {23,22}; i=i+1;
Line(i) = {22,21}; i=i+1;

//-29
Line(i) = {14,24}; i=i+1;
Line(i) = {24,23}; i=i+1;
//-34

//-31
Line(i) = {19,25}; i=i+1;
Line(i) = {25,24}; i=i+1;
//-37

//-33
Line(i) = {18,26}; i=i+1;
Line(i) = {26,25}; i=i+1;
//-39

//Line(i) = {}; i=i+1;

// ******************************** //
// If varying the elements' shapes is acceptable, use the following lines.
// ******************************** //
/*
// Define a line loop with the previously defined lines. Warning: always user counter-clockwise direction to prevent inverted elements. Note: use "-" before line number to reverse its direction.
Line Loop(1) = {1,2,3,4,5,6,7};
Line Loop(2) = {8,9,10,-6,-5,-4,-3};
// Define a surface with the previously defined line loop.
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// Set element size in the neighbourhood of some points.
Mesh.CharacteristicLengthFactor = 1; // Scaling of all elements' sizes.
Characteristic Length {1, 2, 3, 4, 5, 6, 7, 8, 9} = dxGlob;
*/

// ******************************** //
// If a strictly structured (constant dX) mesh is wanted, uncomment and configure the following lines instead.
// ******************************** //

// horizontals
c1[]=Point{1}; c2[]=Point{2}; nelts=Abs(c1[0]-c2[0])/dxGlob;
Transfinite Line{1,3,15,27,35} = nelts/2+1;
c1[]=Point{4}; c2[]=Point{5}; nelts=Abs(c1[0]-c2[0])/dxGlob;
Transfinite Line{11,13,25,33,42} = nelts/2+1;

c1[]=Point{2}; c2[]=Point{3}; nelts=Abs(c1[0]-c2[0])/dxGlob;
Transfinite Line{5,7,18,29,38} = nelts/2+1;
c1[]=Point{3}; c2[]=Point{4}; nelts=Abs(c1[0]-c2[0])/dxGlob;
Transfinite Line{8,10,23,31,40} = nelts/2+1;

//verticals
//c1[]=Point{7}; c2[]=Point{8}; nelts=Abs(c1[1]-c2[1])/dxGlob;
c1[]=Point{13}; c2[]=Point{14}; nelts=Sqrt((c1[0]-c2[0])^2 + (c1[1]-c2[1])^2)/(0.75*dxGlob); // triangle
//c1[]=Point{11}; c2[]=Point{21}; nelts=Sqrt((c1[0]-c2[0])^2 + (c1[1]-c2[1])^2)/(0.75*dxGlob); // triangle
Transfinite Line{16,14,17,19,22,24,28,26,21,20,30,32} = nelts/2+1;
c1[]=Point{1}; c2[]=Point{10}; nelts=Abs(c1[1]-c2[1])/dxGlob;
Transfinite Line{4,2,6,9,12} = nelts/2+1;
c1[]=Point{21}; c2[]=Point{22}; nelts=Abs(c1[1]-c2[1])/dxGlob;
Transfinite Line{36,34,37,39,41} = nelts/2+1;



// Define a line loop with the previously defined lines. Warning: always user counter-clockwise direction to prevent inverted elements. Note: use "-" before line number to reverse its direction.
i = 1;
Line Loop(i) = {1,2,3,4}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {5,6,7,-2}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {8,9,10,-6}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {11,12,13,-9}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1; // triangle
Line Loop(i) = {-3,14,15,16}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-7,17,18,-14}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-17,19,20,21}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-10,22,23,-19}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-13,24,25,-22}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;

Line Loop(i) = {-15,26,27,28}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-18,-21,29,-26}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-20,-23,30,31}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-25,32,33,-30}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-27,34,35,36}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-29,37,38,-34}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-31,39,40,-37}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
Line Loop(i) = {-33,41,42,-39}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;
//Line Loop(i) = {}; Plane Surface(i) = {i}; Transfinite Surface{i}; i=i+1;


// ******************************** //
// Finally, "name" lines and surfaces in order to be able to transfer this to SPECFEM standards.
// Use these whatever was done before.
// ******************************** //
Physical Line("Top") = {35,38,40,42}; // Associate line 4 (RHS) to "Top" (LHS) SPECFEM boundary. Warning: if many lines, list them in counter-clockwise order.
Physical Line("Left") = {4,16,28,36}; // Associate line 4 (RHS) to "Left" (LHS) SPECFEM boundary. Warning: if many lines, list them in counter-clockwise order.
Physical Line("Bottom") = {1,5,8,11}; // Associate line 4 (RHS) to "Bottom" (LHS) SPECFEM boundary. Warning: if many lines, list them in counter-clockwise order.
Physical Line("Right") = {12,24,32,41}; // Associate line 4 (RHS) to "Right" (LHS) SPECFEM boundary. Warning: if many lines, list them in counter-clockwise order.
Physical Surface("M1") = {10,11,12,13,14,15,16,17}; // Associate surface 1 (RHS) to SPECFEM material 1 (LHS).
Physical Surface("M2") = {1,2,3,4,5,6,7,8,9}; // Associate surface 2 (RHS) to SPECFEM material 2 (LHS).
Recombine Surface {1,2,3,4,5,6,7,8,9};
Recombine Surface {10,11,12,13,14,15,16,17};
