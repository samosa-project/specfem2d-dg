clear all;
close all;
clc;

% Choices.
xminmax = [-1,1]*23e3;
depthmax = -15e3; altitudemax = 15e3;
casesToPrepare = 0:4; % 0 = without, 1 = short high, 2 = long high, 3 = short low, 4 = realistic

% Parameters.
addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new'));
thisFolder = [regexprep(mfilename('fullpath'),mfilename,'')];
% NOTE: the script 'produce_cross_section.m'
% (/home/l.martire/Documents/software/MotleyToolbox/meshing/produce_cross_section.m)
% needs to be used to produce a cross-section for the realistic case.
path_to_LITHO = '/home/l.martire/Documents/software/LITHO1.0/access_litho';
EXDIR = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/';
parfile = [thisFolder, filesep, 'parfile_input'];
sourcefile = [thisFolder, filesep, 'source_input'];

% EXAMPLE folders.
prfx_0without = 'mountain_scattering_without_mountains';
prfx_033_high = 'mountain_scattering_with_mountains_0.33L0';
prfx_033_loww = 'mountain_scattering_with_mountains_0.33L0_lower';
prfx_300_high = 'mountain_scattering_with_mountains_3.00L0';
prfx_real = 'mountain_scattering_with_realistic';
source_subfold = thisFolder; % from where copy parfile and sourcefile
dest_subfolds = {prfx_0without, prfx_033_high, prfx_033_loww, prfx_300_high, prfx_real}; % to where send them

% Deduce choice of lat/lon for ground model.
crossSection = load([EXDIR, prfx_real, filesep, 'cross_section.mat']);
lon = mean(crossSection.xq); lat = mean(crossSection.yq);

% Prepare ground model.
f0 = readExampleFiles_extractParam(sourcefile, 'f0', 'float');
prepare_ground

% Deduce L0.
z0 = readExampleFiles_extractParam(sourcefile, 'zs', 'float');
id_up = find(model(:,1)>z0);
id_lo = find(model(:,1)<z0);
if(isempty(id_lo))
  layer_which_has_source = model(id_up, :);
else
  layer_which_has_source = model(id_up+1:id_lo,:);
end
vp_at_source = layer_which_has_source(3);
L0 = vp_at_source/f0;
disp(' ');
disp(['Found $L_0=\SI{',sprintf('%.3f',L0/1e3),'}{\km}$.']);
% disp(['Domain = ']);
disp(['DX at interface asked to be at most ',num2str(dx(end-1)),'.']);

% Produce geo file.
prepare_geofile

propagate_input_files = 0;
if(propagate_input_files)
  % Copy parfile and sourcefile to all subfolders.
  
  % NOTE:
  % for 3L0, dt=3e-4 (cfl 0.429 cf 408826), nstep=140000, runs in ~45 minutes
  % for 0.33, dt=1e-4 (dt=3e-4 causes cfl>1 cf 40882{4,5}), nstep=420000, runs in ??
  % for no topo, dt=3e-4 (cfl 0.29 cf 408827), nstep=140000, should run ~45 minutes
  
  files_to_copy = {'parfile_input', 'source_input'};
  for i=1:numel(dest_subfolds)
    for j=1:numel(files_to_copy)
      command = ['cp ',source_subfold,'/',files_to_copy{j},' ',EXDIR,dest_subfolds{i},filesep,''];
      system(command);
    end
  end
end