clear all;
close all;
clc;

addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new'));

OF_rangeDep = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/LNS_GWB/OUTPUT_FILES_4236576';
OF_rangeIndep = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/LNS_GWB_rangeIndep/OUTPUT_FILES_4236578';
OFs = {OF_rangeIndep, OF_rangeDep};
% IT = 18000;
IT = 25000;
NX = 100*5;
NZ = 500*5;
maxpercenterr = 10;
XLAB = '$x$ [km]';
YLAB = '$z$ [km]';

% Start.
parfile = [OF_rangeIndep, filesep, 'input_parfile'];
bgmodel_indep = [OF_rangeIndep, filesep, 'input_background_model.bin'];
bgmodel_indep_h = [OF_rangeIndep, filesep, 'input_background_model_header.dat'];
bgmodel_rngdep = [OF_rangeDep, filesep, 'input_background_model.bin'];
% bgmodel_rngdep_h = [OF_rangeDep, filesep, 'input_background_model_header.dat'];

DT = readExampleFiles_extractParam(parfile, 'DT', 'float');
BGMODEL_indep = try_make_bg_model_matrix(load_bg_model(bgmodel_indep, bgmodel_indep_h));
% BGMODEL_rngdep = try_make_bg_model_matrix(load_bg_model(bgmodel_rngdep, bgmodel_rngdep_h));
C = mean(sqrt(BGMODEL_indep.ga .* BGMODEL_indep.pr ./ BGMODEL_indep.rh), 'all');
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
[prefix, factor] = prefix_factor_values({Vi{io}.pre});

% Plot.
fig_field = figure('units','normalized','outerposition',[0,0,1,1]);
tightAxes = tight_subplot(2, 2, [0.08,0.12], [0.1, 0.06], [0.055, 0.075]); % gap, marg_h, marg_w
for io = 1:numel(OFs)
  if(io==1)
    maxabsval = max(abs(Vi{io}.pre(Yi>=D_wavefront-3e3 & Yi<=(D_wavefront+1e3)))); % grab amplitude of wavefront
    TIT_add = 'Range Independent Model';
  else
    TIT_add = 'Range Dependent Model';
  end
  
  axes(tightAxes(io));
  pcolor(Xi/1e3, Yi/1e3, factor * Vi{io}.pre); hold on;
  shading interp;
%   daspect([1,1,1]);
  colormaps_fromPython('seismic', 1);
  caxis([-1,1]*maxabsval*factor);
  title({[TIT_add]});
  
  if(io==1)
    ylabel(YLAB);
  else
%     xlabel(XLAB);
    hcb = colorbar;
    ylabel(hcb, ['pressure [',prefix,'Pa]'], 'interpreter', 'latex', 'fontsize', 26);
  end
end

difference_percent = 100*((Vi{2}.pre-Vi{1}.pre)./maxabsval); % reference to rangedep (io=1)
% fig_diff = figure('units','normalized','outerposition',[0,0,1,1]);
axes(tightAxes(3));
pcolor(Xi/1e3, Yi/1e3, difference_percent); hold on;
shading interp;
% daspect([1,1,1]);
colormaps_fromPython('bwr', 1);
hcb = colorbar;
caxis([-1,1]*maxpercenterr);
ylabel(hcb, ['\% of direct wavefront'], 'interpreter', 'latex', 'fontsize', 26);
xlabel(XLAB); ylabel(YLAB);
title({['Pressure Field Difference']});

axes(tightAxes(4));
title({['Ground Pressure Difference']});

% set(tightAxes(1:3), 'xlim', XLIM, 'ylim', YLIM);
set(tightAxes(1:2), 'xticklabel', {});
set(tightAxes(2), 'yticklabel', {});
tightAxes(2).Position = [tightAxes(2).Position(1), tightAxes(1).Position(2), tightAxes(1).Position(3:4)];
tightAxes(3).Position = [tightAxes(1).Position(1), tightAxes(3).Position(2), tightAxes(1).Position(3:4)];
tightAxes(4).Position = [tightAxes(2).Position(1), tightAxes(3).Position(2), tightAxes(1).Position(3:4)];

customSaveFig(fig_field, [figfieldpath], extToSave, 9999);

