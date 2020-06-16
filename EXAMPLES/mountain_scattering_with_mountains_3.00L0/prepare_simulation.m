clear all;
close all;
clc;

% Choices.
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

% Produce geo file.
prepare_geofile

% TODO: copy parfile and sourcefile to all subfolders.
for i=1:numel(dest_subfolds)
  command = ['cp ',source_subfold,'/parfile_input ',EXDIR,dest_subfolds{i},filesep,''];
  system(command);
end