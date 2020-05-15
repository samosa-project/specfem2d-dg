clear all;
close all;
clc;

OF_rangeDep = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/LNS_GWB/OUTPUT_FILES_4236576';
OF_rangeIndep = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/LNS_GWB_rangeIndep/OUTPUT_FILES_4236578';
OFs = {OF_rangeDep, OF_rangeIndep};
IT = 18000;
% IT = 25000;
NX = 100*5;
NZ = 500*5;
maxpercenterr = 10;
XLAB = '$x$ [km]';
YLAB = '$z$ [km]';

% Start.
DT = 0.005;
C = 371;
time = IT*DT;
D_wavefront = C*time;
timestr=['at $t=',sprintf('%.0f',time),'$~s'];

thisFolder = [regexprep(mfilename('fullpath'),mfilename,'')];
output_figures_folder = [thisFolder, filesep, 'FIGURES'];
figfieldpath = [output_figures_folder, filesep, 'field'];
figdiffpath = [output_figures_folder, filesep, 'field_difference'];
extToSave = {'eps'};

% Read.
disp(['[',mfilename,'] Reading dumps (may be long).']);
Vi = {};
for io = 1:numel(OFs)
  OF = OFs{io};
  [X, Y, V] = readDumpsUnique(OF, IT, 0);
  [Xi, Yi, Vi{io}] = interpDumps(X, Y, V, NX, NZ, 0);
end

% Plot.
fig_field = {};
for io = 1:numel(OFs)
  if(io==1)
    maxabsval = max(abs(Vi{io}.pre(Yi>=D_wavefront-3e3 & Yi<=(D_wavefront+1e3)))); % grab amplitude of wavefront
    addendum = '_rD';
    TIT_add = ', Range Dependent';
  else
    addendum = '_rI';
    TIT_add = ', Range Independent';
  end

  fig_field{io} = figure('units','normalized','outerposition',[0,0,1,1]);
  pcolor(Xi/1e3, Yi/1e3, Vi{io}.pre); hold on;
  shading interp; daspect([1,1,1]);
  colormaps_fromPython('seismic', 1);
  colorbar;
  caxis([-1,1]*maxabsval);
  xlabel(XLAB); ylabel(YLAB); title({['Field ',timestr, TIT_add]});
  customSaveFig(fig_field{io}, [figfieldpath, addendum], extToSave, 9999);
end

difference_percent = 100*((Vi{1}.pre-Vi{2}.pre)./maxabsval); % reference to rangedep (io=1)
fig_diff = figure('units','normalized','outerposition',[0,0,1,1]);
pcolor(Xi/1e3, Yi/1e3, difference_percent); hold on;
shading interp; daspect([1,1,1]);
colormaps_fromPython('bwr', 1);
colorbar;
caxis([-1,1]*maxpercenterr);
xlabel(XLAB); ylabel(YLAB); title({['Field Difference ',timestr,' [\%]']});
customSaveFig(fig_diff, [figdiffpath], extToSave, 9999);