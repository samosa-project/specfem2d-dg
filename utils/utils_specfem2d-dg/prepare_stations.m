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

[SPCFMEXloc] = setup_overall();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run-specific.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simulationfolder = input(['[',mfilename,'] Path to simulation folder > '],'s');
if(not(simulationfolder(end)==filesep)); simulationfolder=[simulationfolder,filesep]; end;

typeeee=-1;
while(not(numel(typeeee)==1 & ismember(typeeee,[0,1,101,102,103,201,3])))
  disp(['[',mfilename,'] Type of simulation?']);
  disp([blanks(length(mfilename)+2),'   1   AboveAndBelowGround_Periodic_WithTilt']);
  disp([blanks(length(mfilename)+2),'   101 AboveAndBelowGround_Periodic_WithTilt__MarsInSIGHT']);
  disp([blanks(length(mfilename)+2),'   102 AboveAndBelowGround_Periodic_WithTilt__MarsInSIGHT_Impact']);
  disp([blanks(length(mfilename)+2),'   103 AboveAndBelowGround_Periodic_WithTilt__MarsInSIGHT_Waveguide']);
  disp([blanks(length(mfilename)+2),'   2   AboveGroundSimple']);
  disp([blanks(length(mfilename)+2),'   201 AboveGroundSimple_Microbaroms']);
  disp([blanks(length(mfilename)+2),'   3   tir_de_mine']);
  typeeee=input([blanks(length(mfilename)+2),' > ']);
end
switch (typeeee)
  case 1
    [xminmax, zminmax, interface, Xsource, debfin, d, name] = AboveAndBelowGround_Periodic_WithTilt(simulationfolder);
  case 101
    [xminmax, zminmax, interface, Xsource, debfin, d, name] = AboveAndBelowGround_Periodic_WithTilt__MarsInSIGHT(simulationfolder);
  case 102
    [xminmax, zminmax, interface, Xsource, debfin, d, name] = AboveAndBelowGround_Periodic_WithTilt__MarsInSIGHT_Impact(simulationfolder);
  case 103
    [xminmax, zminmax, interface, Xsource, debfin, d, name] = AboveAndBelowGround_Periodic_WithTilt__MarsInSIGHT_Waveguide(simulationfolder);
  case 2
    [xminmax, zminmax, interface, Xsource, debfin, d, name] = AboveGroundSimple(simulationfolder);
  case 201
    [xminmax, zminmax, interface, Xsource, debfin, d, name] = AboveGroundSimple_Microbaroms(simulationfolder);
  case 3
    [xminmax, zminmax, interface, Xsource, debfin, d, name] = tir_de_mine(simulationfolder);
  otherwise
    error('kuk');
end
% [xminmax, zminmax, interface, Xsource, debfin, d, name] = mars_insight([0,500]);
% [xminmax, zminmax, interface, Xsource, debfin, d, name] = tir_de_mine();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Automatic.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = plot_domain_and_stations(xminmax, zminmax, interface, Xsource, debfin, d);

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
% Functions.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N] = plot_domain_and_stations(xminmax, zminmax, interface, Xsource, debfin, delta)
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