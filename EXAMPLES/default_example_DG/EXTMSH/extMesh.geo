// Display. //---------------//
General.Axes = 0; // Show axes.
Mesh.Lines = 1; // Display lines.
Mesh.LineNumbers = 1; // Display line numbers.
Mesh.Points = 1; // Display points.
Mesh.PointNumbers = 1; // Display point numbers
Mesh.SurfaceEdges = 1; // Display edges.
Mesh.SurfaceFaces = 1; // Display surfaces.
Mesh.SurfaceNumbers = 1; // Display surface numbers.

// Elements. //--------------//
Mesh.SubdivisionAlgorithm = 1; // Meshing algorithm: all quads.
Mesh.RecombineAll = 1; // Recombine all triangular elements.

// Meshing algorithm. //-----//
Mesh.Algorithm = 9; // 1: MeshAdapt. 5: Delaunay for quads. 9: structured (experimental).
Mesh.ElementOrder = 1; // Element order.

// Geometry. //--------------//
xmin=-10;
xmax=-xmin;
zmin=0;
zmax=20;

// Define elements' length.
Mesh.CharacteristicLengthFactor = 1; // Scaling of all elements' sizes.
dxbot=1;
dxtop=3;

// Define geometry keypoints.
Point(1) = {xmin,zmin,0};
Point(2) = {xmax,zmin,0};
Point(3) = {xmax,zmax,0};
Point(4) = {xmin,zmax,0};
// Define lines joining points. Warning: prefer counter-clockwise direction to prevent inverted elements.
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
// Define a line loop with the previously defined lines. Warning: prefer counter-clockwise direction to prevent inverted elements.
Line Loop(5) = {1,2,3,4};
// Define a surface with the previously defined line loop.
Plane Surface(6) = {5};

Physical Line("Top") = {4};
Physical Line("Left") = {3};
Physical Line("Bottom") = {1};
Physical Line("Right") = {2};

Physical Surface("M1") = {6};

Characteristic Length {1, 2} = 2*dxbot; // Set mesh characteristic size at some points. *2 to prevent it from being too small due to first subdivision.
Characteristic Length {4, 3} = 2*dxtop;

Recombine Surface {6};
