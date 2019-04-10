% Author:        Léo Martire.
% Description:   Prepares stations for SPECFEM2D runs.
% Notes:         TODO.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

clear all;
close all;
clc;
addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run-specific.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simulationfolder = input(['[',mfilename,'] Path to simulation folder > '],'s');
if(not(simulationfolder(end)==filesep)); simulationfolder=[simulationfolder,filesep]; end;

typeeee=-1;
while(not(numel(typeeee)==1 & ismember(typeeee,[0,1,101])))
  disp(['[',mfilename,'] Type of simulation?']);
  disp([blanks(length(mfilename)+2),'     0 AboveAndBelowGround_Periodic_WithTilt']);
  disp([blanks(length(mfilename)+2),'     1 AboveGroundSimple']);
  disp([blanks(length(mfilename)+2),'   101 AboveGroundSimple_Microbaroms']);
  typeeee=input([blanks(length(mfilename)+2),' > ']);
end
switch (typeeee)
  case 0
    [xminmax, zminmax, interface, Xsource, debfin, d, name] = AboveAndBelowGround_Periodic_WithTilt(simulationfolder);
  case 1
    [xminmax, zminmax, interface, Xsource, debfin, d, name] = AboveGroundSimple(simulationfolder);
  case 101
    [xminmax, zminmax, interface, Xsource, debfin, d, name] = AboveGroundSimple_Microbaroms(simulationfolder);
  otherwise
    error('kuk');
end
% [xminmax, zminmax, interface, Xsource, debfin, d, name] = mars_insight([0,500]);
% [xminmax, zminmax, interface, Xsource, debfin, d, name] = tir_de_mine();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Automatic.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = get_domain_and_stations(xminmax, zminmax, interface, Xsource, debfin, d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ask if ok.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ok = -1;
while(not(numel(ok) == 1 & ismember(ok, [0, 1])))
  ok = input(['[', mfilename, '] Are those ', num2str(sum(N)), ' stations ok (0 for no, 1 for yes)? > ']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display copy-pastable.      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ok)
  clc;
  disp(['nreceiversets = ', num2str(numel(d)), ' # (total number of stations: ', num2str(sum(N)), ')']);
  disp(['# Orientation']);
  disp(['anglerec              = 0.0d0   # Angle to rotate components at receivers.']);
  disp(['rec_normal_to_surface = .false. # Base anglerec normal to surface (external mesh and curve file needed).']);
  for i = 1:numel(d)
    disp(['# ', name{i}]);
    disp(['nrec = ', num2str(N(i))]);
    disp(['xdeb = ', char(sprintf('%17.16g', debfin(i, 1, 1)))]);
    disp(['zdeb = ', char(sprintf('%17.16g', debfin(i, 2, 1)))]);
    disp(['xfin = ', char(sprintf('%17.16g', debfin(i, 1, 2)))]);
    disp(['zfin = ', char(sprintf('%17.16g', debfin(i, 2, 2)))]);
    disp(['record_at_surface_same_vertical = .false.']);
  end

  disp(' ');
else
  disp(['[', mfilename, '] Stations not ok, stopping script.']);
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run configurations.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xminmax, zminmax, interface, Xsource, debfin, d, name] = AboveGroundSimple_Microbaroms(simulationfolder)
  parfile        = [simulationfolder, 'parfile_input'];
  sourcefile     = [simulationfolder, 'source_input'];
  interfacesfile = [simulationfolder, 'interfaces_input'];
  [xminmax, zminmax, Xsource] = getRunParameters(parfile, sourcefile, interfacesfile);
  STAT_DX = 15e3;
  STAT_Z = [1, 15, 30]*1e3; % altitude of stations.
  groundclearance = 3;
  interface = [-1e9, 1e9;0, 0];
  isitok=-1;
  while(not(ismember(isitok,[0,1])))
    isitok=input(['[', mfilename, '] Stations'' DX planned to be ',num2str(STAT_DX),' [m]. Is that ok (0 for no, 1 for yes)? > ']);
  end
  if(isitok==0)
    error(['[',mfilename,', ERROR] Parameter was not ok, re-choose parametrisation in script.']);
  end
  isitok=-1;
  while(not(ismember(isitok,[0,1])))
    isitok=input(['[', mfilename, '] Stations'' Z planned to be ',num2str(STAT_Z),' [m]. Is that ok (0 for no, 1 for yes)? > ']);
  end
  if(isitok==0)
    error(['[',mfilename,', ERROR] Parameter was not ok, re-choose parametrisation in script.']);
  end
  lid = 0;
%   lid = lid+1;
%   d(lid) = 0;
%   debfin(lid, 1, :) = Xsource(1)*[1, 1]; % xdeb xfin
%   debfin(lid, 2, :) = Xsource(2)*[1, 1]; % zdeb zfin
%   name{lid} = ['source'];
  for i=1:numel(STAT_Z)
    lid = lid+1;
    d(lid) = STAT_DX;
    debfin(lid, 1, :) = [xminmax(1)+d(lid), Xsource(1)]; % xdeb xfin
    debfin(lid, 2, :) = [1, 1]*STAT_Z(i); % zdeb zfin
    name{lid} = ['horizontal altitude ',num2str(STAT_Z),' [m]'];
  end
  lid = lid+1;
  d(lid) = 1000;
  debfin(lid, 1, :) = [-1, 1]*8e3; % xdeb xfin
  debfin(lid, 2, :) = [1, 1]*groundclearance; % zdeb zfin
  name{lid} = ['microbarom monitoring'];
end

function [xminmax, zminmax, interface, Xsource, debfin, d, name] = AboveGroundSimple(simulationfolder)
  parfile        = [simulationfolder, 'parfile_input'];
  sourcefile     = [simulationfolder, 'source_input'];
  interfacesfile = [simulationfolder, 'interfaces_input'];
  [xminmax, zminmax, Xsource] = getRunParameters(parfile, sourcefile, interfacesfile);
  STAT_DX = 15e3;
  STAT_Z = [1, 15, 30]*1e3; % altitude of stations.
  interface = [-1e9, 1e9;0, 0];
  isitok=-1;
  while(not(ismember(isitok,[0,1])))
    isitok=input(['[', mfilename, '] Stations'' DX planned to be ',num2str(STAT_DX),' [m]. Is that ok (0 for no, 1 for yes)? > ']);
  end
  if(isitok==0)
    error(['[',mfilename,', ERROR] Parameter was not ok, re-choose parametrisation in script.']);
  end
  isitok=-1;
  while(not(ismember(isitok,[0,1])))
    isitok=input(['[', mfilename, '] Stations'' Z planned to be ',num2str(STAT_Z),' [m]. Is that ok (0 for no, 1 for yes)? > ']);
  end
  if(isitok==0)
    error(['[',mfilename,', ERROR] Parameter was not ok, re-choose parametrisation in script.']);
  end
  lid = 0;
  lid = lid+1;
  d(lid) = 0;
  debfin(lid, 1, :) = Xsource(1)*[1, 1]; % xdeb xfin
  debfin(lid, 2, :) = Xsource(2)*[1, 1]; % zdeb zfin
  name{lid} = ['source'];
  for i=1:numel(STAT_Z)
    lid = lid+1;
    d(lid) = STAT_DX;
    debfin(lid, 1, :) = [xminmax(1)+d(lid), xminmax(2)-d(lid)]; % xdeb xfin
    debfin(lid, 2, :) = [1, 1]*STAT_Z(i); % zdeb zfin
    name{lid} = ['horizontal altitude ',num2str(STAT_Z),' [m]'];
  end
end

function [xminmax, zminmax, interface, Xsource, debfin, d, name] = AboveAndBelowGround_Periodic_WithTilt(simulationfolder)
  parfile        = [simulationfolder, 'parfile_input'];
  sourcefile     = [simulationfolder, 'source_input'];
  interfacesfile = [simulationfolder, 'interfaces_input'];
%   xminmax = [extractParamFromInputFile(parfile, 'xmin', 'float'), extractParamFromInputFile(parfile, 'xmax', 'float')];
%   Xsource = [extractParamFromInputFile(sourcefile, 'xs', 'float'), extractParamFromInputFile(sourcefile, 'zs', 'float')];
%   [zmin, zmax] = extractZminZmaxFromInterfacesFile(interfacesfile);
%   zminmax = [zmin,zmax];
  [xminmax, zminmax, Xsource] = getRunParameters(parfile, sourcefile, interfacesfile);
  spacingstations = 300;
  isitok=-1;
  while(not(ismember(isitok,[0,1])))
    isitok=input(['[', mfilename, '] Stations planned to be spaced by ',num2str(spacingstations),' [m]. Is that ok (0 for no, 1 for yes)? > ']);
  end
  if(isitok==0)
    error(['[',mfilename,', ERROR] Parameter was not ok, re-choose parametrisation in script.']);
  end
  
  interface = [-1e9, 1e9;0, 0];
  ground_clearance = 5; % altitude/depth of the ground stations.
  isitok=-1;
  while(not(ismember(isitok,[0,1])))
    isitok=input(['[', mfilename, '] Ground clearance (above and below) planned to be ',num2str(ground_clearance),' [m]. Is that ok (0 for no, 1 for yes)? > ']);
  end
  if(isitok==0)
    error(['[',mfilename,', ERROR] Parameter was not ok, re-choose parametrisation in script.']);
  end
  shift_for_tiltcomputation = 5; % horizontal shift for stations used for tilt computation
  disp(['[', mfilename, '] Horizontal shift for tilt computation planned to be ',num2str(shift_for_tiltcomputation),' [m].']);

  lid = 0;
  
  lid = lid+1; idhoriz = lid;
  d(lid) = 0;
  debfin(lid, 1, :) = Xsource(1)*[1, 1]; % xdeb xfin
  debfin(lid, 2, :) = Xsource(2)*[1, 1]; % zdeb zfin
  name{lid} = ['source'];
  
  lid = lid+1; idhoriz = lid;
  d(lid) = spacingstations;
  debfin(lid, 1, :) = [xminmax(1)+d(lid), xminmax(2)-d(lid)]; % xdeb xfin
  debfin(lid, 2, :) = [1, 1]*ground_clearance; % zdeb zfin
  name{lid} = ['horizontal over ground'];
  lid = lid+1; idhorizvz = lid;
  d(lid) = d(idhoriz);
  debfin(lid, 1, :) = debfin(idhoriz, 1, :); % xdeb xfin same as horizontal over ground
  debfin(lid, 2, :) = -1*[1, 1]*ground_clearance; % zdeb zfin
  name{lid} = ['horizontal under ground'];
  lid = lid+1;
  d(lid) = d(idhoriz);
  debfin(lid, :, :) = debfin(idhorizvz, :, :); % all same as horizontal under ground
  debfin(lid, 1, :) = debfin(lid, 1, :)-shift_for_tiltcomputation; % but shifted to the left
  name{lid} = ['horizontal under ground shifted -', num2str(shift_for_tiltcomputation), 'm for tilt computation'];
  lid = lid+1;
  d(lid) = d(idhoriz);
  debfin(lid, :, :) = debfin(idhorizvz, :, :); % all same as horizontal under ground
  debfin(lid, 1, :) = debfin(lid, 1, :)+shift_for_tiltcomputation; % but shifted to the right
  name{lid} = ['horizontal under ground shifted ', num2str(shift_for_tiltcomputation), 'm for tilt computation'];
  
  % on periodic boundary condition
  lid = lid+1; d(lid) = 2*ground_clearance;
  debfin(lid, 1, :) = [xminmax(2), xminmax(2)]; debfin(lid, 2, :) = [1, -1]*ground_clearance;
  name{lid} = ['over & under ground @ right boundary'];
  lid = lid+1; d(lid) = diff(xminmax)-2*shift_for_tiltcomputation;
  debfin(lid, 1, :) = [xminmax(1)+shift_for_tiltcomputation,xminmax(2)-shift_for_tiltcomputation]; debfin(lid, 2, :) = -1*[1, 1]*ground_clearance;
  name{lid} = ['left&right of right boundary for tilt computation'];
end

function [xminmax, zminmax, interface, Xsource, debfin, d, name] = tir_de_mine()
  Xsource = [0, -25];
  xminmax = [-440, 240];
  zminmax = [-190, 150];
  ISAEIS2_z = 65;
  ISAEIS3_z = 104;
  to_gla08_z = (195.-167.);
  to_gla08_r = 187.;
  to_gla08_x = (to_gla08_r^2.-to_gla08_z^2.)^0.5;
  xr0 = -to_gla08_x/to_gla08_r; zr0 = -to_gla08_z/to_gla08_r; zx0 = -to_gla08_z/to_gla08_x;
  interface = [-400, 10;400*zx0, -10*zx0];

  lid = 0;
  lid = lid+1; d(lid) = 0; r = 0;   name{lid} = ['above source (@r = 0)', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 100; name{lid} = ['GLA04 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 193; name{lid} = ['GLA05 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 150; name{lid} = ['GLA06 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 5;   name{lid} = ['GLA07 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['GLA08 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0; idGLA08 = lid;
  lid = lid+1; d(lid) = 0; r = 185; name{lid} = ['GLA09 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 50;  name{lid} = ['GLA11 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  debfin(:, 2, 1) = debfin(:, 2, 1)+1; % beginning above ground
  debfin(:, 2, 2) = debfin(:, 2, 2)-1; % end under ground
  d = d+2; % separate correctly
  lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['ISAE IS Sensor 2, above GLA08 (@r = ', num2str(r), '), z = +', num2str(ISAEIS2_z), ' (ISAEIS ground is GLA08 z = +1m)']; debfin(lid, :, :) = debfin(idGLA08, :, :); debfin(lid, 2, :) = debfin(lid, 2, :)+ISAEIS2_z;
  lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['ISAE IS Sensor 3, above GLA08 (@r = ', num2str(r), '), z = +', num2str(ISAEIS3_z), ' (ISAEIS ground is GLA08 z = +1m)']; debfin(lid, :, :) = debfin(idGLA08, :, :); debfin(lid, 2, :) = debfin(lid, 2, :)+ISAEIS3_z;
  lid = lid+1; d(lid) = 1; r = 0; name{lid} = ['Monitor Source']; debfin(lid, 1, :) = Xsource(1)*[1, 1]; debfin(lid, 2, :) = (Xsource(2)+d(lid))*[1, 1];
end






















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N] = get_domain_and_stations(xminmax, zminmax, interface, Xsource, debfin, delta)
  if(not( numel(xminmax) == numel(zminmax) & numel(zminmax) == 2))
    error('kek');
  end
  if(not(size(debfin, 1) == numel(delta)))
    error('kek');
  end
  lid = numel(delta);
  figure('units', 'normalized', 'outerposition', [0 0 1 1]);
  set(gca, 'Color', 'k');
  set(gca, 'GridColor', 'white');
  rectangle('Position', [xminmax(1), zminmax(1), diff(xminmax), diff(zminmax)], 'edgecolor', 'white', 'linewidth', 1); hold on;
  line(interface(1, :), interface(2, :), 'color', 'white', 'linewidth', 1);
  pct = 0.05; margin = min(pct*diff(xminmax), pct*diff(zminmax));
  axis([xminmax(1)-margin, xminmax(2)+margin, zminmax(1)-margin, zminmax(2)+margin]);
  grid on; box on;
  plot(Xsource(1), Xsource(2), 'r+');
  for i = 1:lid
    Lline = norm(diff(debfin(i, :, :), 1, 3));
    if(Lline == 0 || delta(i) == 0)
      N(i) = 1;
    else
      N(i) = ceil(Lline/delta(i))+1;
    end
    N(i)
    xz = [linspace(debfin(i, 1, 1), debfin(i, 1, 2), N(i)); ...
        linspace(debfin(i, 2, 1), debfin(i, 2, 2), N(i))]
    plot(xz(1, :), xz(2, :), 'g.');
  end
  daspect([1 1 1]);
  xlabel(['$x$ [m]']); ylabel(['$z$ [m]']);
  set(gca, 'ticklabelinterpreter', 'latex');
end