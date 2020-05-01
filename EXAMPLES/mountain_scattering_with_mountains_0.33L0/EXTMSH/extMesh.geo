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
Mesh.Algorithm = 5; // 1: MeshAdapt. 5: Delaunay for quads. 8: Delaunay (experimental). 9: structured (packing of parallelograms, experimental).
Mesh.ElementOrder = 1; // Element order.
Mesh.RandomFactor = 0.000100; // Perturbate every points positions in order to avoid 3 aligned points. Default is at 1e-9, maximum is 1e-3.
// Geometry. //--------------//
// Physical parameters.
dxAirBot = 93.500000;
dxAirTop = 132.000000;
dxGrdInts1 = 2592.000000;
dxGrdInts2 = 1620.000000;
dxGrdInts3 = 762.000000;
// Define geometry keypoints.
Point(1) = {-23000.000000, -41140.000000, 0}; // corner
Point(2) = {23000.000000, -41140.000000, 0}; // corner
Point(3) = {-23000.000000, -31140.000000, 0}; // corner
Point(4) = {23000.000000, -31140.000000, 0}; // corner
Point(5) = {-23000.000000, -16560.000000, 0}; // corner
Point(6) = {23000.000000, -16560.000000, 0}; // corner
Point(7) = {-23000.000000, 0.000000, 0}; // corner
Point(8) = {23000.000000, 0.000000, 0}; // corner
Point(9) = {18000.000000, 0.000000, 0}; // mountain point 
Point(10) = {17500.000000, 1500.000000, 0}; // mountain point 
Point(11) = {17000.000000, 0.000000, 0}; // mountain point 
Point(12) = {16500.000000, 1500.000000, 0}; // mountain point 
Point(13) = {16000.000000, 0.000000, 0}; // mountain point 
Point(14) = {15500.000000, 1500.000000, 0}; // mountain point 
Point(15) = {15000.000000, 0.000000, 0}; // mountain point 
Point(16) = {14500.000000, 1500.000000, 0}; // mountain point 
Point(17) = {14000.000000, 0.000000, 0}; // mountain point 
Point(18) = {13500.000000, 1500.000000, 0}; // mountain point 
Point(19) = {13000.000000, 0.000000, 0}; // mountain point 
Point(20) = {12500.000000, 1500.000000, 0}; // mountain point 
Point(21) = {12000.000000, 0.000000, 0}; // mountain point 
Point(22) = {11500.000000, 1500.000000, 0}; // mountain point 
Point(23) = {11000.000000, 0.000000, 0}; // mountain point 
Point(24) = {10500.000000, 1500.000000, 0}; // mountain point 
Point(25) = {10000.000000, 0.000000, 0}; // mountain point 
Point(26) = {9500.000000, 1500.000000, 0}; // mountain point 
Point(27) = {9000.000000, 0.000000, 0}; // mountain point 
Point(28) = {8500.000000, 1500.000000, 0}; // mountain point 
Point(29) = {8000.000000, 0.000000, 0}; // mountain point 
Point(30) = {7500.000000, 1500.000000, 0}; // mountain point 
Point(31) = {7000.000000, 0.000000, 0}; // mountain point 
Point(32) = {6500.000000, 1500.000000, 0}; // mountain point 
Point(33) = {6000.000000, 0.000000, 0}; // mountain point 
Point(34) = {5500.000000, 1500.000000, 0}; // mountain point 
Point(35) = {5000.000000, 0.000000, 0}; // mountain point 
Point(36) = {4500.000000, 1500.000000, 0}; // mountain point 
Point(37) = {4000.000000, 0.000000, 0}; // mountain point 
Point(38) = {3500.000000, 1500.000000, 0}; // mountain point 
Point(39) = {3000.000000, 0.000000, 0}; // mountain point 
Point(40) = {2500.000000, 1500.000000, 0}; // mountain point 
Point(41) = {2000.000000, 0.000000, 0}; // mountain point 
Point(42) = {1500.000000, 1500.000000, 0}; // mountain point 
Point(43) = {1000.000000, 0.000000, 0}; // mountain point 
Point(44) = {500.000000, 1500.000000, 0}; // mountain point 
Point(45) = {0.000000, 0.000000, 0}; // mountain point 
Point(46) = {-500.000000, 1500.000000, 0}; // mountain point 
Point(47) = {-1000.000000, 0.000000, 0}; // mountain point 
Point(48) = {-1500.000000, 1500.000000, 0}; // mountain point 
Point(49) = {-2000.000000, 0.000000, 0}; // mountain point 
Point(50) = {-2500.000000, 1500.000000, 0}; // mountain point 
Point(51) = {-3000.000000, 0.000000, 0}; // mountain point 
Point(52) = {-3500.000000, 1500.000000, 0}; // mountain point 
Point(53) = {-4000.000000, 0.000000, 0}; // mountain point 
Point(54) = {-4500.000000, 1500.000000, 0}; // mountain point 
Point(55) = {-5000.000000, 0.000000, 0}; // mountain point 
Point(56) = {-5500.000000, 1500.000000, 0}; // mountain point 
Point(57) = {-6000.000000, 0.000000, 0}; // mountain point 
Point(58) = {-6500.000000, 1500.000000, 0}; // mountain point 
Point(59) = {-7000.000000, 0.000000, 0}; // mountain point 
Point(60) = {-7500.000000, 1500.000000, 0}; // mountain point 
Point(61) = {-8000.000000, 0.000000, 0}; // mountain point 
Point(62) = {-8500.000000, 1500.000000, 0}; // mountain point 
Point(63) = {-9000.000000, 0.000000, 0}; // mountain point 
Point(64) = {-9500.000000, 1500.000000, 0}; // mountain point 
Point(65) = {-10000.000000, 0.000000, 0}; // mountain point 
Point(66) = {-10500.000000, 1500.000000, 0}; // mountain point 
Point(67) = {-11000.000000, 0.000000, 0}; // mountain point 
Point(68) = {-11500.000000, 1500.000000, 0}; // mountain point 
Point(69) = {-12000.000000, 0.000000, 0}; // mountain point 
Point(70) = {-12500.000000, 1500.000000, 0}; // mountain point 
Point(71) = {-13000.000000, 0.000000, 0}; // mountain point 
Point(72) = {-13500.000000, 1500.000000, 0}; // mountain point 
Point(73) = {-14000.000000, 0.000000, 0}; // mountain point 
Point(74) = {-14500.000000, 1500.000000, 0}; // mountain point 
Point(75) = {-15000.000000, 0.000000, 0}; // mountain point 
Point(76) = {-15500.000000, 1500.000000, 0}; // mountain point 
Point(77) = {-16000.000000, 0.000000, 0}; // mountain point 
Point(78) = {-16500.000000, 1500.000000, 0}; // mountain point 
Point(79) = {-17000.000000, 0.000000, 0}; // mountain point 
Point(80) = {-17500.000000, 1500.000000, 0}; // mountain point 
Point(81) = {-18000.000000, 0.000000, 0}; // mountain point 
Point(82) = {-23000.000000, 15000.000000, 0};
Point(83) = {23000.000000, 15000.000000, 0};
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
Line(10) = {9, 10};
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 13};
Line(14) = {13, 14};
Line(15) = {14, 15};
Line(16) = {15, 16};
Line(17) = {16, 17};
Line(18) = {17, 18};
Line(19) = {18, 19};
Line(20) = {19, 20};
Line(21) = {20, 21};
Line(22) = {21, 22};
Line(23) = {22, 23};
Line(24) = {23, 24};
Line(25) = {24, 25};
Line(26) = {25, 26};
Line(27) = {26, 27};
Line(28) = {27, 28};
Line(29) = {28, 29};
Line(30) = {29, 30};
Line(31) = {30, 31};
Line(32) = {31, 32};
Line(33) = {32, 33};
Line(34) = {33, 34};
Line(35) = {34, 35};
Line(36) = {35, 36};
Line(37) = {36, 37};
Line(38) = {37, 38};
Line(39) = {38, 39};
Line(40) = {39, 40};
Line(41) = {40, 41};
Line(42) = {41, 42};
Line(43) = {42, 43};
Line(44) = {43, 44};
Line(45) = {44, 45};
Line(46) = {45, 46};
Line(47) = {46, 47};
Line(48) = {47, 48};
Line(49) = {48, 49};
Line(50) = {49, 50};
Line(51) = {50, 51};
Line(52) = {51, 52};
Line(53) = {52, 53};
Line(54) = {53, 54};
Line(55) = {54, 55};
Line(56) = {55, 56};
Line(57) = {56, 57};
Line(58) = {57, 58};
Line(59) = {58, 59};
Line(60) = {59, 60};
Line(61) = {60, 61};
Line(62) = {61, 62};
Line(63) = {62, 63};
Line(64) = {63, 64};
Line(65) = {64, 65};
Line(66) = {65, 66};
Line(67) = {66, 67};
Line(68) = {67, 68};
Line(69) = {68, 69};
Line(70) = {69, 70};
Line(71) = {70, 71};
Line(72) = {71, 72};
Line(73) = {72, 73};
Line(74) = {73, 74};
Line(75) = {74, 75};
Line(76) = {75, 76};
Line(77) = {76, 77};
Line(78) = {77, 78};
Line(79) = {78, 79};
Line(80) = {79, 80};
Line(81) = {80, 81};
Line(82) = {81, 7};
Line(83) = {8, 83};
Line(84) = {83, 82};
Line(85) = {82, 7};
Line(86) = {7, 5};
Line Loop(1) = {83,84,85,-82,-81,-80,-79,-78,-77,-76,-75,-74,-73,-72,-71,-70,-69,-68,-67,-66,-65,-64,-63,-62,-61,-60,-59,-58,-57,-56,-55,-54,-53,-52,-51,-50,-49,-48,-47,-46,-45,-44,-43,-42,-41,-40,-39,-38,-37,-36,-35,-34,-33,-32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9};
Line Loop(4) = {1, 2, 3, 4};
Line Loop(3) = {-3, 5, 6, 7};
Line Loop(2) = {-6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 86}; 
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Mesh.CharacteristicLengthFactor = 1; // Scaling of all elements' sizes. 
Characteristic Length {7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81} = 2*dxAirBot; 
Characteristic Length {82, 83} = 2*dxAirTop; 
Characteristic Length {1, 2} = 2*dxGrdInts1; 
Characteristic Length {3, 4} = 2*dxGrdInts2; 
Characteristic Length {5, 6} = 2*dxGrdInts3; 
Physical Line("Top") = {84};
Physical Line("Left") = {85, 86,7,4};
Physical Line("Bottom") = {1};
Physical Line("Right") = {2,5,8, 83}; 
Physical Surface("M1") = {1};
Physical Surface("M2") = {2};
Physical Surface("M3") = {3};
Physical Surface("M4") = {4};
Recombine Surface {1,1,2,3,4};
