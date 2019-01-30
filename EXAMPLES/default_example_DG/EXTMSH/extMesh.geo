// Display. //---------------//
General.Axes = 0; // Show axes.
Mesh.Lines = 1; // Display lines.
Geometry.LineNumbers = 1; // Display line numbers.
Mesh.LineNumbers = 1; // Display line numbers.
Mesh.Points = 1; // Display points.
Geometry.PointNumbers = 1; // Display point numbers
Mesh.PointNumbers = 1; // Display point numbers
Mesh.SurfaceEdges = 1; // Display edges.
Mesh.SurfaceFaces = 1; // Display surfaces.
Geometry.SurfaceNumbers = 1; // Display surface numbers.
Mesh.SurfaceNumbers = 1; // Display surface numbers.

// Elements. //--------------//
Mesh.SubdivisionAlgorithm = 1; // Meshing algorithm: all quads.
Mesh.RecombineAll = 1; // Recombine all triangular elements.

// Meshing algorithm. //-----//
Mesh.Algorithm = 9; // 1: MeshAdapt. 5: Delaunay for quads. 9: structured (experimental).
Mesh.ElementOrder = 2; // Element order.
Mesh.RandomFactor=1e-3; // Perturbate every points' positions in order to avoid 3 aligned points. Default is at 1e-9, maximum is 1e-3.

// Geometry. //--------------//
xmin=-50;
xmax=-xmin;
zmin=-50;
zmax=-zmin;

// Define elements' length.
dxbot=10;
dxtop=20;

// Define geometry keypoints.
Point(1) = {xmin,zmin,0};
Point(2) = {xmax,zmin,0};
Point(3) = {xmax,zmax,0};
Point(4) = {xmin,zmax,0};
// Define lines joining points. Warning: prefer counter-clockwise direction to prevent inverted elements. Warning: do not duplicate lines if two or more materials are needed.
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// ******************************** //
// If varying the elements' shapes is acceptable, use the following lines.
// ******************************** //
// Define a line loop with the previously defined lines. Warning: always user counter-clockwise direction to prevent inverted elements. Note: use "-" before line number to reverse its direction.
Line Loop(1) = {1,2,3,4};
// Define a surface with the previously defined line loop.
Plane Surface(1) = {1};
// Set element size in the neighbourhood of some points.
Mesh.CharacteristicLengthFactor = 1; // Scaling of all elements' sizes.
Characteristic Length {1, 2} = 2*dxbot;
Characteristic Length {4, 3} = 2*dxtop;

// ******************************** //
// If a strictly structured (constant dX) mesh is wanted, uncomment and configure the following lines instead.
// ******************************** //
/*
nelts=5; // Choose number of elements to create on each of the following lines.
Transfinite Line{1,3} = nelts/2+1;
nelts=8; // Choose number of elements to create on each of the following lines.
Transfinite Line{2,4} = nelts/2+1;
// Define a line loop with the previously defined lines. Warning: always user counter-clockwise direction to prevent inverted elements. Note: use "-" before line number to reverse its direction.
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Transfinite Surface{1};
*/

// ******************************** //
// Finally, "name" lines and surfaces in order to be able to transfer this to SPECFEM standards.
// Use these whatever was done before.
// ******************************** //
Physical Line("Top") = {4}; // Associate line 4 (RHS) to "Top" (LHS) SPECFEM boundary. Warning: if many lines, list them in counter-clockwise order.
Physical Line("Left") = {3}; // Associate line 4 (RHS) to "Left" (LHS) SPECFEM boundary. Warning: if many lines, list them in counter-clockwise order.
Physical Line("Bottom") = {1}; // Associate line 4 (RHS) to "Bottom" (LHS) SPECFEM boundary. Warning: if many lines, list them in counter-clockwise order.
Physical Line("Right") = {2}; // Associate line 4 (RHS) to "Right" (LHS) SPECFEM boundary. Warning: if many lines, list them in counter-clockwise order.
Physical Surface("M1") = {1}; // Associate surface 1 (RHS) to SPECFEM material 1 (LHS).
Recombine Surface {1};
