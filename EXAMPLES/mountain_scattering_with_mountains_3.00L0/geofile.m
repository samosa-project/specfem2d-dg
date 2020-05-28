clear all;
close all;
clc;

% Prepare ground model.
depthmax = -40e3;
lat = 0.5*(45.05+44.78);
lon = 0.5*(5.08+7.42);
[s, res] = system(['/home/l.martire/Documents/software/LITHO1.0/access_litho -p ',sprintf('%.4f',lat),' ',sprintf('%.4f',lat),'']);
res = split(res, newline);
res(end) = [];
model = zeros(1, 6);
c = 1;
for i = numel(res):-1:1
  cur_z_r_vp_vs_qk_qm = str2num(regexp(res{i}, ['( +-?[0-9]+\.[0-9]* ){6}'], 'match', 'once'));
  if(not(isempty(cur_z_r_vp_vs_qk_qm)))
    model(c, :) = cur_z_r_vp_vs_qk_qm;
    c = c + 1;
  end
end
model(:,1) = - (model(:,1) - min(model(:,1)));
model(logical([0; diff(model(:,1))~=0]), :) = []; % make sure no duplicate interface
model(model(:,1)<depthmax, :) = []; % remove layers too deep
model = flipud(model);
itooclose = find(diff(model(:,1))<1000); model(itooclose+[1], :) = []; % find layers too close to each other, remove the top one

% Deduce interfaces and dx.
f0 = 2; % [hz]
nptsperwavelength = 2;
% interfaces = [-41140, -31140, -16560, 0, 15e3];
interfaces = [depthmax, model(:,1)', 15e3];
igrd = find(interfaces==0);
% dx = [2700, 1800, 800, 110, 132];
dx = [(model(:,3)'/f0)/nptsperwavelength, 110, 132];
xminmax = [-1,1]*23e3;

% Print it.
NMATERIALS = size(model, 1)+1;
disp('# Number of models.');
disp(['nbmodels                        = ',num2str(NMATERIALS),'']);
disp('#   acoustic:    model_number 1  rho  Vp   0   0   0   QKappa Qmu 0   0   0   0    0    0');
disp('#   elastic:     model_number 1  rho  Vp   Vs  0   0   QKappa Qmu 0   0   0   0    0    0');
% disp('# Remark: for elastic media, Vp and Vs must be the unrelaxed velocities (following the viscoelastic behaviour at infinite frequency). Attenuation can be taken into account by setting the parameters from the ''Attenuation'' section above as wanted.');
% disp('1 1 1.164  349.    0. 0 0   10.   10. 0 0 0 0 0 0 # air');
fmt = '%7.2f';
disp([num2str(1),' 1 ',sprintf(fmt, 1.4),'  ',sprintf(fmt, 340),' ',sprintf(fmt, 0),' 0 0 ',sprintf(fmt, 9999),' ',sprintf(fmt, 9999),' 0 0 0 0 0 0 # Air.']);
disp(['# LITHO1.0 (doi 10.1002/2013JB010626) model computed at (',sprintf('%.4f',lat),'°N, ',sprintf('%.4f',lon),'°E):']);
model = flipud(model);
for i=1:(NMATERIALS-1)
  rho = model(i, 2);
  vp = model(i, 3);
  vs = model(i, 4);
  qk = model(i, 5);
  qm = model(i, 6);
%   rcsq=(vs./vp).^2;
%   qp = ((1-rcsq).*qk.^(-1) + rcsq.*qm.^(-1)).^(-1);
%   qs = qm;
  disp([num2str(i+1),' 1 ',sprintf(fmt, rho),'  ',sprintf(fmt, vp),' ',sprintf(fmt, vs),' 0 0 ',sprintf(fmt, qk),' ',sprintf(fmt, qm),' 0 0 0 0 0 0']);
%   3 1 2830  6500. 3710. 0 0 2704. 2850. 0 0 0 0 0 0 # middle crust: 16560-31140 (Qp=2750, Qs=2850)
%   4 1 2920  6900. 3930. 0 0 2954. 3100. 0 0 0 0 0 0 # lower crust: 31140-44100 (Qp=3000, Qs=3100)
end

% Produce geo file.
% wo_0__w033L_1__w3L_2__w033Lheight_3 = 0;
% wo_0__w033L_1__w3L_2__w033Lheight_3 = 1;
% wo_0__w033L_1__w3L_2__w033Lheight_3 = 2;
% wo_0__w033L_1__w3L_2__w033Lheight_3 = 3;
% wo_0__w033L_1__w3L_2__w033Lheight_3 = 4; % realistic mountains
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
      outputFile = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_realistic/EXTMSH/extMesh.geo';
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
      load('/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_realistic/cross_section.mat');
      x = rq-0.5*range(rq);
      xstop = max(x)+1e3;
      xminmax = [-1,1]*xstop;
      z = elvq;
      
      buf=5e3;
      apol = 0.5*(1-cos((x-min(x))*2*pi/(2*buf))) .* (x<=min(x)+buf) + (x>min(x)+buf); %plot(x, apol);
      apor = 0.5*(1-cos((max(x)-x)*2*pi/(2*buf))) .* (x>=max(x)-buf) + (x<max(x)-buf); %plot(x, apor);
      apo = apol.*apor; %plot(x, apo);
      z = z .* apo;
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

  % Build line loops, plane surfaces, and physical surfaces.
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
    
    materialnumber = NMATERIALS-gill+1;
    printLineLoopPlaneSurfPhySurf(fid, gill, [botlist, right, toplist, left], materialnumber);
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
  
  nvertiline = size(lines_list_verti, 1)/2;
  left = lines_list_verti(1:nvertiline, 1);
  right = lines_list_verti((nvertiline+1):end, 1);
  
  printPhyLine(fid, 'Top', lines_list_horiz(end, 1));
  printPhyLine(fid, 'Bottom', lines_list_horiz(1, 1));
  printPhyLine(fid, 'Left', left);
  printPhyLine(fid, 'Right', right);
  
  printRecombine(fid, [1:gill-1]);
  
  
%   Physical Line("Top") = {20};
%   Physical Line("Left") = {21, 22,7,4};
%   Physical Line("Bottom") = {1};
%   Physical Line("Right") = {2,5,8, 19};
%   Recombine Surface {1,1,2,3,4};

  fclose(fid);
end

function printPhyLine(fid, name, lineList)
  list = sprintf('%d,', lineList);
  list(end)=[];
  fprintf(fid, ['Physical Line("',name,'") = {',list,'};\n']);
end
function printRecombine(fid, surfList)
  list = sprintf('%d,', surfList);
  list(end)=[];
  fprintf(fid, ['Recombine Surface {',list,'};\n']);
end
function printCharacteristicLength(fid, ptslist, var)
  list = sprintf('%d,', ptslist);
  list(end)=[];
  fprintf(fid, ['Characteristic Length {',list,'} = 2*',var,';\n']);
end
function printLineLoopPlaneSurfPhySurf(fid, i, lines, materialnumber)
  list = sprintf('%d,', lines);
  list(end)=[];
  fprintf(fid, ['Line Loop(',num2str(i),') = {',list,'};\n']);
  fprintf(fid, 'Plane Surface(%d) = {%d};\n', i, i);
  fprintf(fid, 'Physical Surface("M%d") = {%d};\n', materialnumber, i);
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
  displaySurfaceEdge = 1;
  displaySurfaceFaces = 1;
  displayLineNumber = 1;
  displayPointsNumber = 0;
  displayMeshPoints = 0;
  displayMeshPointsNumber = 0;
  displayMeshLineNumber = 0;
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