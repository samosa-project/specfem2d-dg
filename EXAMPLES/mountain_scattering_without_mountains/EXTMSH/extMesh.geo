// Display. //---------------//
General.Axes = 0; // Show axes.
Mesh.Lines = 1; // Display lines.
Geometry.LineNumbers = 1; // Display line numbers.
Mesh.LineNumbers = 0; // Display line numbers.
Mesh.Points = 1; // Display points.
Geometry.PointNumbers = 1; // Display point numbers.
Mesh.PointNumbers = 0; // Display point numbers
Mesh.SurfaceEdges = 1; // Display edges.
Mesh.SurfaceFaces = 1; // Display surfaces.
Geometry.SurfaceNumbers = 1; // Display surface numbers.
Mesh.SurfaceNumbers = 0; // Display surface numbers.
// Elements. //--------------//
Mesh.SubdivisionAlgorithm = 1; // Meshing algorithm: all quads.
Mesh.RecombineAll = 1; // Recombine all triangular elements.
// Meshing algorithm. //-----//
Mesh.Algorithm = 5; // 1: MeshAdapt. 5: Delaunay for quads. 9: structured (experimental).
Mesh.ElementOrder = 1; // Element order.
Mesh.RandomFactor = 0.000100; // Perturbate every points positions in order to avoid 3 aligned points. Default is at 1e-9, maximum is 1e-3.
// Geometry. //--------------//
// Physical parameters.
dxAirBot = 110.000000;
dxAirTop = 132.000000;
dxGrdInts1 = 2592.000000;
dxGrdInts2 = 1620.000000;
dxGrdInts3 = 762.000000;
// Define geometry keypoints.
Point(1) = {-28000.000000, -41140.000000, 0}; // corner
Point(2) = {28000.000000, -41140.000000, 0}; // corner
Point(3) = {-28000.000000, -31140.000000, 0}; // corner
Point(4) = {28000.000000, -31140.000000, 0}; // corner
Point(5) = {-28000.000000, -16560.000000, 0}; // corner
Point(6) = {28000.000000, -16560.000000, 0}; // corner
Point(7) = {-28000.000000, 0.000000, 0}; // corner
Point(8) = {28000.000000, 0.000000, 0}; // corner
Point(9) = {0.000000, 0.000000, 0}; // mountain point 
Point(10) = {-28000.000000, 30000.000000, 0};
Point(11) = {28000.000000, 30000.000000, 0};
// Define lines joining points. Warning: prefer counter-clockwise direction to prevent inverted elements. Warning: do not duplicate lines if two or more materials are needed.
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Line(5) = {4, 6};
Line(6) = {6, 5};
Line(7) = {5, 3};
Line(8) = {6, 8};
Line(9) = {8, 9};
Line(10) = {9, 7};
Line(11) = {8, 11};
Line(12) = {11, 10};
Line(13) = {10, 7};
Line(14) = {7, 5};
Line Loop(1) = {11,12,13,-10,-9};
Line Loop(4) = {1, 2, 3, 4};
Line Loop(3) = {-3, 5, 6, 7};
Line Loop(2) = {-6, 8, 9, 10, 14}; 
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Mesh.CharacteristicLengthFactor = 1; // Scaling of all elements' sizes. 
Characteristic Length {7, 8, 9} = 2*dxAirBot; 
Characteristic Length {10, 11} = 2*dxAirTop; 
Characteristic Length {1, 2} = 2*dxGrdInts1; 
Characteristic Length {3, 4} = 2*dxGrdInts2; 
Characteristic Length {5, 6} = 2*dxGrdInts3; 
Physical Line("Top") = {12};
Physical Line("Left") = {13, 14,7,4};
Physical Line("Bottom") = {1};
Physical Line("Right") = {2,5,8, 11}; 
Physical Surface("M1") = {1};
Physical Surface("M2") = {2};
Physical Surface("M3") = {3};
Physical Surface("M4") = {4};
Recombine Surface {1,1,2,3,4};
