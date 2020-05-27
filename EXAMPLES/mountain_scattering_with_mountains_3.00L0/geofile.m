clear all;
close all;
clc;


interfaces = [-41140, -31140, -16560, 0, 15e3]; igrd = find(interfaces==0);
dx = [2700, 1800, 800, 110, 132];
xminmax = [-1,1]*23e3;

% wo_0__w033L_1__w3L_2__w033Lheight_3 = 0;
% wo_0__w033L_1__w3L_2__w033Lheight_3 = 1;
% wo_0__w033L_1__w3L_2__w033Lheight_3 = 2;
% wo_0__w033L_1__w3L_2__w033Lheight_3 = 3;
% wo_0__w033L_1__w3L_2__w033Lheight_3 = 4; % realistic, pyrenees

for wo_0__w033L_1__w3L_2__w033Lheight_3 = 0:4
  switch(wo_0__w033L_1__w3L_2__w033Lheight_3)
    case 1
      % output to "with 0.33L" folder
      outputFile = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_mountains_0.33L0/EXTMSH/extMesh.geo';
      LTopo = 1000;
      nPerio = 36;
      peakHeight = 1500;
      nptperperio = 30;
    case 2
      outputFile = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_mountains_3.00L0/EXTMSH/extMesh.geo';
      LTopo = 9000;
      nPerio = 4;
      peakHeight = 1500;
      nptperperio = 30;
    case 3
      % output to "with 0.33L height adjusted" folder
      outputFile = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_mountains_0.33L0_lower/EXTMSH/extMesh.geo';
      LTopo = 1000;
      nPerio = 36;
      peakHeight = (1500/9000)*LTopo; % keep the same angle as in the 3L0 simulation
      nptperperio = 30;
    case 0
      % output to "without" folder
      outputFile = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_without_mountains/EXTMSH/extMesh.geo';
      LTopo = -1;
      peakHeight = -1;
      nPerio = -1;
      nptperperio = -1;
    case 4
      outputFile = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_pyrenees/EXTMSH/extMesh.geo';
      LTopo = -1;
      nPerio = -1;
      peakHeight = -1;
      nptperperio = -1;
    otherwise
      error('eeorkrkerlk');
  end

  % build surface
  switch(wo_0__w033L_1__w3L_2__w033Lheight_3)
    case 0
      x = 0;
      z = 0;
    case 4
      % load
      load('/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_pyrenees/cross_section_pyrenees.mat');
      x = rq-0.5*range(rq);
      xstop = max(x)+1e3;
      xminmax = [-1,1]*xstop;
      z = elvq;
      dx = [2700, 1800, 800, 200, 250];
  %     pause
    otherwise
      xstop = nPerio*0.5*LTopo;
      n = nPerio*nptperperio;
      x = linspace(-xstop, xstop, n);
      z = peakHeight * 0.5 * (1-cos(x*2*pi/LTopo)) .* (x>=-xstop) .* (x<=xstop);
  end
  surface_xz = [x; z];

  % Mesh.
  meshAlgorithm = 5;

  % Creation. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fid = fopen(outputFile, 'w');

  printDisplay(fid);
  printMesh(fid, meshAlgorithm);

  % Print some variables.
  fprintf(fid, '// Some Variables. //--------//\n');
  dx_int = [];
  for i=1:numel(interfaces)
    dx_int(i).var = ['dx_int_', num2str(i)];
    dx_int(i).val = dx(i);
    defineVariable(fid, dx_int(i).var, dx_int(i).val);
  end

  fprintf(fid, '// Geometry. //--------------//\n');

  % Points.
  gip = 1;
  for i=1:numel(interfaces)
    for j=1:2
      printPoint(fid, gip, xminmax(j), interfaces(i));
      gip = gip + 1;
    end
    if(i==igrd)
      % save last two ids as being borders of interface
      surface_ends = [gip-1, gip]-1;
    end
  end

  % Lines.
  gil = 1;
  lids = ((1:numel(interfaces))*2-1);
  lines_list_horiz = []; % save horizontal lines
  for i=lids
    if(i==surface_ends(1))
      % don't do this one, save relative position in list of lines
      lines_list_horiz = [lines_list_horiz; [-1, -1, -1]];
    else
      printLine(fid, gil, i, i+1);
      lines_list_horiz = [lines_list_horiz; [gil, i, i+1]];
      gil = gil + 1;
    end
  end
  vids = [lids(1:end-1), lids(1:end-1)+1];
  lines_list_verti = [];
  for i=vids
    lines_list_verti = [lines_list_verti; [gil, i, i+2]];
    printLine(fid, gil, i, i+2);
    gil = gil+1;
  end

  % custom surface interface, points
  pts_list_csi = [];
  for i=1:size(surface_xz, 2)
    pts_list_csi = [pts_list_csi, gip];
    printPoint(fid, gip, surface_xz(1, i), surface_xz(2, i));
    gip = gip + 1;
  end

  % custom surface interface, lines
  lines_list_csi = [];
  printLine(fid, gil, surface_ends(1), pts_list_csi(1));
  lines_list_csi = [lines_list_csi; [gil, surface_ends(1), pts_list_csi(1)]];
  gil = gil+1;
  for i=1:(numel(pts_list_csi)-1)
    printLine(fid, gil, pts_list_csi(i), pts_list_csi(i+1));
    lines_list_csi = [lines_list_csi; [gil, pts_list_csi(i), pts_list_csi(i+1)]];
    gil = gil + 1;
  end
  printLine(fid, gil, pts_list_csi(end), surface_ends(2));
  lines_list_csi = [lines_list_csi; [gil, pts_list_csi(end), surface_ends(2)]];
  gil = gil+1;

  % Build line loops and physical surfaces.
  gill=1;
  for l=1:(size(lines_list_horiz, 1)-1)
    hl = lines_list_horiz(l, 1);
    potbot = lines_list_horiz(l, :);
    pottop = lines_list_horiz(l+1, :);
    if(potbot(1)==-1)
      botlist = lines_list_csi(:, 1)';
    else
      botlist = potbot(1);
    end
    if(pottop(1)==-1)
      toplist = -fliplr(lines_list_csi(:, 1)');
    else
      toplist = -pottop(1);
    end
    % find left/right
    if(pottop(1)==-1)
      vl_id_l = find(lines_list_verti(:,3)==surface_ends(1)); % find the line finishing at the leftmost surface interface point
      vl_id_r = find(lines_list_verti(:,3)==surface_ends(2)); % find the line finishing at the rightmost surface interface point
    else
      vl_id_l = find(lines_list_verti(:, 3)==pottop(2)); % find a vertical line ending where the top line starts
      vl_id_r = find(lines_list_verti(:, 3)==pottop(3)); % find a vertical line ending where the top line end
    end
    left = -lines_list_verti(vl_id_l, 1);
    right = lines_list_verti(vl_id_r, 1);
    printLineLoopAndPhySurf(fid, gill, [botlist, right, toplist, left]);
    gill = gill + 1;
  end

  % Specify mesh size.
  for l=1:size(lines_list_horiz, 1)
    hl = lines_list_horiz(l, :);
    if(hl(1)~=-1)
      printCharacteristicLength(fid, hl(2:3), dx_int(l).var);
    end
  end
  printCharacteristicLength(fid, sort(unique([pts_list_csi, surface_ends])), dx_int(igrd).var);

  fclose(fid);
end

function printCharacteristicLength(fid, ptslist, var)
  list = sprintf('%d,', ptslist);
  list(end)=[];
  fprintf(fid, ['Characteristic Length {',list,'} = 2*',var,';\n']);
end
function printLineLoopAndPhySurf(fid, i, lines)
  list = sprintf('%d,', lines);
  list(end)=[];
  fprintf(fid, ['Line Loop(',num2str(i),') = {',list,'};\n']);
  fprintf(fid, 'Plane Surface(%d) = {%d};\n', i, i);
end
function printLine(fid, i, p1, p2)
  fprintf(fid, 'Line(%d) = {%d, %d};\n', i, p1 , p2);
end
function printPoint(fid, i, x, z)
  fprintf(fid, 'Point(%d) = {%f, %f, 0};\n', i, x, z);
end

function defineVariable(fid, var, value)
  fprintf(fid, [var, ' = %f;\n'], value);
end

function printMesh(fid, meshAlgorithm)
  meshRandomFactor = 1e-4;

  fprintf(fid, '// Elements. //--------------//\n');
  fprintf(fid, 'Mesh.SubdivisionAlgorithm = 1; // Meshing algorithm: all quads.\n');
  fprintf(fid, 'Mesh.RecombineAll = 1; // Recombine all triangular elements.\n');
  fprintf(fid, '// Meshing algorithm. //-----//\n');
  fprintf(fid, 'Mesh.Algorithm = %d; // 1: MeshAdapt. 5: Delaunay for quads. 8: Delaunay (experimental). 9: structured (packing of parallelograms, experimental).\n', meshAlgorithm);
  fprintf(fid, 'Mesh.ElementOrder = 1; // Element order.\n');
  fprintf(fid, 'Mesh.RandomFactor = %f; // Perturbate every points positions in order to avoid 3 aligned points. Default is at 1e-9, maximum is 1e-3.\n', meshRandomFactor);
end

function printDisplay(fid)
  showAxes = 1;
  displayLine = 1;
  displayLineNumber = 1;
  displayMeshLineNumber = 0;
  displayMeshPoints = 0;
  displayPointsNumber = 1;
  displayMeshPointsNumber = 0;
  displaySurfaceEdge = 1;
  displaySurfaceFaces = 1;
  displaySurfaceNumber = 1;
  displayMeshSurfaceNumber = 0;

  fprintf(fid, '// Display. //---------------//\n');
  fprintf(fid, 'General.Axes = %d; // Show axes.\n', showAxes);
  fprintf(fid, 'Mesh.Lines = %d; // Display lines.\n', displayLine);
  fprintf(fid, 'Geometry.LineNumbers = %d; // Display line numbers.\n', displayLineNumber);
  fprintf(fid, 'Mesh.LineNumbers = %d; // Display line numbers.\n', displayMeshLineNumber);
  fprintf(fid, 'Mesh.Points = %d; // Display points.\n', displayMeshPoints);
  fprintf(fid, 'Geometry.PointNumbers = %d; // Display point numbers.\n', displayPointsNumber);
  fprintf(fid, 'Mesh.PointNumbers = %d; // Display point numbers\n', displayMeshPointsNumber);
  fprintf(fid, 'Mesh.SurfaceEdges = %d; // Display edges.\n', displaySurfaceEdge);
  fprintf(fid, 'Mesh.SurfaceFaces = %d; // Display surfaces.\n', displaySurfaceFaces);
  fprintf(fid, 'Geometry.SurfaceNumbers = %d; // Display surface numbers.\n', displaySurfaceNumber);
  fprintf(fid, 'Mesh.SurfaceNumbers = %d; // Display surface numbers.\n', displayMeshSurfaceNumber);
end