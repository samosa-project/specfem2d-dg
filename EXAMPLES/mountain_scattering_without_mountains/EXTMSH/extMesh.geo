// This geofile was created by '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_mountains_3.00L0/prepare_geofile'.
// No topography.


// Display. //---------------//
General.Axes = 1; // Show axes.
Mesh.Lines = 1; // Display lines.
Geometry.LineNumbers = 1; // Display line numbers.
Mesh.LineNumbers = 0; // Display line numbers.
Mesh.Points = 0; // Display points.
Geometry.PointNumbers = 0; // Display point numbers.
Mesh.PointNumbers = 0; // Display point numbers
Mesh.SurfaceEdges = 1; // Display edges.
Mesh.SurfaceFaces = 1; // Display surfaces.
Geometry.SurfaceNumbers = 1; // Display surface numbers.
Mesh.SurfaceNumbers = 0; // Display surface numbers.
// Elements. //--------------//
Mesh.SubdivisionAlgorithm = 1; // Meshing algorithm: all quads.
Mesh.RecombineAll = 1; // Recombine all triangular elements.
// Meshing algorithm. //-----//
Mesh.Algorithm = 5; // 1: MeshAdapt. 5: Delaunay for quads. 8: Delaunay (experimental). 9: structured (packing of parallelograms, experimental).
Mesh.ElementOrder = 1; // Element order.
Mesh.RandomFactor = 0.000100; // Perturbate every points positions in order to avoid 3 aligned points. Default is at 1e-9, maximum is 1e-3.
// Some Variables. //--------//
dx_int_1 = 1519.515000;
dx_int_2 = 110.000000;
dx_int_3 = 132.000000;
// Geometry. //--------------//
Point(1) = {-23000.000000, -15000.000000, 0};
Point(2) = {23000.000000, -15000.000000, 0};
Point(3) = {-23000.000000, 0.000000, 0};
Point(4) = {23000.000000, 0.000000, 0};
Point(5) = {-23000.000000, 15000.000000, 0};
Point(6) = {23000.000000, 15000.000000, 0};
Line(1) = {1, 2};
Line(2) = {5, 6};
Line(3) = {1, 3};
Line(4) = {3, 5};
Line(5) = {2, 4};
Line(6) = {4, 6};
Point(7) = {0.000000, 0.000000, 0};
Line(7) = {3, 7};
Line(8) = {7, 4};
Line Loop(1) = {1,5,-8,-7,-3};
Plane Surface(1) = {1};
Physical Surface("M2") = {1};
Line Loop(2) = {7,8,6,-2,-4};
Plane Surface(2) = {2};
Physical Surface("M1") = {2};
Characteristic Length {1,2} = 2*dx_int_1;
Characteristic Length {5,6} = 2*dx_int_3;
Characteristic Length {3,4,7} = 2*dx_int_2;
Physical Line("Top") = {2};
Physical Line("Bottom") = {1};
Physical Line("Left") = {3,4};
Physical Line("Right") = {5,6};
Recombine Surface {1,2};
