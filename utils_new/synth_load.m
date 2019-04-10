% Author:        Léo Martire.
% Description:   Loads SPECFEM2D OUTPUT_FILES synthetics, and does a bunch of thing with them.
% Notes:         TODO.
%
% Usage:
%  1) Specify script parameters in the 'Parameters' section.
%  2) Configure fig_title, rootd, and OFd in the 'OUTPUT_FILES' section.
%     Eventually re-specify the script's parameters for each run.
%  3) Answer input prompts as they go.

clc;
% clear all;
clear('Zamp','Ztime'); disp(['[',mfilename,'] Cleared Zamp and Ztime variables.']);
% close all;
format compact;
set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'defaultTextFontSize', 12); set(0, 'defaultAxesFontSize', 12);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/Atmospheric_Models');
addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools'); % truncToShortest, readAndSubsampleSynth
addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/standalone');
SPCFMEXloc = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% You can re-set those parameters inline for each OUTPUT_FILES directory (see "OUTPUT_FILES location" section below) individually.
rescale_factor = 1; % Rescaling factor, applied on all stations (confirmation will be asked).
convert_to_relative_coords = 0; pos_interface = 0; % Convert to relative coordinates: x = 0 above source, z = 0 on surface (defined by pos_interface).
plot_amplitude = 0; % Plot amplitude (0 for no, 1 for yes sorted by x, 2 for yes sorted by z, 3 for yes sorted by d)?
subsample = 0; wanted_dt = 1; % Sub-sample? Useful for lengthy seismograms. If set to 1, sub-sample so that final time sampling is as parametrised by wanted_dt.
type_display = 2; % Quantity to display (should be the same as the seismotype variable in parfile). 1 = {displacement for non-DG, velocity for DG}. 2 = {velocity for non-DG, pressure perturbation [Pa] for DG}.
% Unknown (for direct plots only). Note that for type_display == 2 and stations in DG zones, pressure perturbation [Pa] is saved both in BXX and BXZ files.
% unknown = 'BXX'; % _x.
unknown = 'BXZ'; % _z.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT_FILES location.       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mars AGW.
fig_title = strcat('Mars Coupling');
% rootd = strcat(SPCFMEXloc,'mars_insight_incidence/'); OFd = strcat(rootd, 'OUTPUT_FILES_151319_20h_f3_larger/'); subsample = 1; wanted_dt = 0.01;
% rootd = strcat(SPCFMEXloc,'mars_insight/'); OFd = strcat(rootd, 'OUTPUT_FILES_151120_z830_f0p1_crashed_but_later/'); subsample = 1; wanted_dt = 0.01;
% rootd = strcat(SPCFMEXloc,'mars_insight_incidence/'); OFd = strcat(rootd, 'OUTPUT_FILES_150395_20h_f3/'); subsample = 1; wanted_dt = 0.01;
rootd = strcat(SPCFMEXloc,'mars_insight/'); OFd = strcat(rootd, 'OUTPUT_FILES_1633618_z800/'); subsample = 1; wanted_dt = 0.01;
% rootd = strcat(SPCFMEXloc,'mars_insight/'); OFd = strcat(rootd, 'OUTPUT_FILES_1601166_z12k/'); subsample = 1; wanted_dt = 0.01;
% rootd = strcat(SPCFMEXloc,'mars_insight/'); OFd = strcat(rootd, 'OUTPUT_FILES_1538139_22h/'); subsample = 1; wanted_dt = 0.01;
% rootd = strcat(SPCFMEXloc,'mars_insight/'); OFd = strcat(rootd, 'OUTPUT_FILES_1529789_20h_cleanusable/'); subsample = 1; wanted_dt = 0.01;
% rootd = strcat(SPCFMEXloc,'mars_insight/'); OFd = strcat(rootd, 'OUTPUT_FILES_1529411_interrupted/');
% rootd = strcat(SPCFMEXloc,'mars_insight/'); OFd = strcat(rootd, 'OUTPUT_FILES_1479218_clean/');
% rootd = strcat(SPCFMEXloc,'mars_insight_cut/'); OFd = strcat(rootd, 'OUTPUT_FILES_1484867/');
% rootd = strcat(SPCFMEXloc,'mars_insight_cut/'); OFd = strcat(rootd, 'OUTPUT_FILES_1486113/');

% TNTGlanes
% fig_title = strcat('Tirs de Mine Glanes');
% rootd = strcat(SPCFMEXloc,'tir_de_mine/'); OFd = strcat(rootd, 'OUTPUT_FILES_1508049_40hz/');
% rootd = strcat(SPCFMEXloc,'tir_de_mine/'); OFd = strcat(rootd, 'OUTPUT_FILES_1479543_clean/');
% rootd = strcat(SPCFMEXloc,'tirdemine_75040_redone/'); OFd = strcat(rootd, 'OUTPUT_FILES_1240682_nospurious_butbadgeom/');
% rootd = strcat(SPCFMEXloc,'tir_de_mine/'); OFd = strcat(rootd, 'OUTPUT_FILES_1216334/');
% rootd = strcat(SPCFMEXloc,'tir_de_mine/'); OFd = strcat(rootd, 'OUTPUT_FILES_1204916_spuriousreflexions/');
% rootd = strcat(SPCFMEXloc,'tir_de_mine/'); OFd = strcat(rootd, 'OUTPUT_FILES_75040/');
% rootd = strcat(SPCFMEXloc,'tir_de_mine/'); OFd = strcat(rootd, 'OUTPUT_FILES_74752/');
% rootd = strcat(SPCFMEXloc,'tntglanes_10/'); OFd = strcat(rootd, 'OUTPUT_FILES_full/');
% rootd = strcat(SPCFMEXloc,'tntglanes_10/'); OFd = strcat(rootd, 'OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'tntglanes_10/'); OFd = strcat(rootd, 'OUTPUT_FILES_long600hz/');
% rootd = strcat(SPCFMEXloc,'tntglanes_10/'); OFd = strcat(rootd, 'OUTPUT_FILES_long300hz/');

% Microbaroms ULDB.
% fig_title = strcat('Microbaroms, (49N, 178W), 6:00 UT');
% rootd = strcat(SPCFMEXloc,'mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_1560350_crash@20600/'); subsample = 1; wanted_dt = 0.01; % same as 1560350 but fns
% rootd = strcat(SPCFMEXloc,'mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_1560541_crashfns/'); % same as 1560350 but fns
% rootd = strcat(SPCFMEXloc,'mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_1560350_crash/'); % same as 1560541 but lns
% rootd = strcat(SPCFMEXloc,'mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_1447857_longestyet/');
% rootd = strcat(SPCFMEXloc,'mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_1206217/'); % Basically same as 1205575.
% rootd = strcat(SPCFMEXloc,'mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_1205575/');
% rootd = strcat(SPCFMEXloc,'mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_1204148_LNS/');
% rootd = strcat(SPCFMEXloc,'mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_1203633_FNS/');
% rootd = strcat(SPCFMEXloc,'mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_74710/');
% rootd = strcat(SPCFMEXloc,'mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_74565/');
% rootd = strcat(SPCFMEXloc,'mb_huge/'); OFd = strcat(rootd, 'OUTPUT_FILES_672048/');
% rootd = strcat(SPCFMEXloc,'mb_huge/'); OFd = strcat(rootd, 'OUTPUT_FILES_642746/');

% PML.
% fig_title = strcat('LNS ABC');
% rootd = strcat(SPCFMEXloc,'test_realstretching_planewave/'); OFd = strcat(rootd, '/OUTPUT_FILES/'); fig_title=[fig_title,' planewave'];
% rootd = strcat(SPCFMEXloc,'test_realstretching_square_buffs_w000/'); OFd = strcat(rootd, '/OUTPUT_FILES_155014_e0p2/'); fig_title=[fig_title,' buffw0 e0.2'];
% rootd = strcat(SPCFMEXloc,'test_realstretching_square_buffs_w000/'); OFd = strcat(rootd, '/OUTPUT_FILES_154970/'); fig_title=[fig_title,' buffw0 e1e-4'];
% rootd = strcat(SPCFMEXloc,'test_realstretching_square_buffs_w100/'); OFd = strcat(rootd, '/OUTPUT_FILES_154971/'); fig_title=[fig_title,' buffw100 e1e-4'];
% rootd = strcat(SPCFMEXloc,'test_realstretching_square_nobuffs_w000/'); OFd = strcat(rootd, '/OUTPUT_FILES_154980/'); fig_title=[fig_title,' nobuffw0 e1e-4'];
% rootd = strcat(SPCFMEXloc,'test_realstretching_square_nobuffs_w100/'); OFd = strcat(rootd, '/OUTPUT_FILES_154979/'); fig_title=[fig_title,' nobuffw100 e1e-4'];

% rootd = strcat(SPCFMEXloc,'test_pml/'); OFd = strcat(rootd, 'OUTPUT_FILES_alright/');

% Validation LNS.
% fig_title = strcat('Validation LNS');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_rho_M0_dx1gmsh_wow/'); % without wind in the source term
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_rho_M.3_dx1gmsh_wow/'); % without wind in the source term
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_rho_M0_dx1gmsh/'); % bad
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_rho_M.3_dx1gmsh/'); % bad
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_E_M0_dx1/');

% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.3/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0_FNS/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.3_FNS/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.3_0.5cfl/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.3_FNS_0.5cfl/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.15/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.15_FNS/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0_refined/')
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.3_refined/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0_gmsh/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.3_gmsh/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0_gmsh_refined_121954/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.3_gmsh_refined_1217753/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_corrected/'); % FIRST RUN AFTER PATCH
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.3_corrected/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_gmsh/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.3_gmsh/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_gmshrefined_1218660/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.3_gmshrefined_1218665/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0_FNS_testlambda0.2/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.3_FNS_testlambda0.2/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0_FNS_testlambda1.2/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.3_FNS_testlambda1.2/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0_testlambda0.2/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.3_testlambda0.2/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0_testlambda1.2/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.3_testlambda1.2/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0_chglambda/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.3_chglambda/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0_chglambda2woperturb/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M = 0.3_chglambda2woperturb/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_1222576_M0_refined/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_1222565_M.3_refined/');

% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_dx1_cfl.245/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.3_dx1_cfl.245/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.6_dx1_cfl.245/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_dx.5_cfl.49/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.6_dx.5_cfl.49/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_dx1_cfl.245_FNS/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.3_dx1_cfl.245_FNS/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_dx1.33_cfl.184_FNS/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.3_dx1.33_cfl.184_FNS/'); % alright base result
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M-.3_dx1_cfl.245/');
% rootd = strcat(SPCFMEXloc,'stationspositionz/'); OFd = strcat(rootd, 'OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_dx1gmshstruct_cfl.440/'); % bad result, source was offset to the right
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_dx1gmshstruct_cfl.440_spreadssf/'); % alright result
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_dx.5gmshstruct_cfl.424_spreadssf/'); % very good result
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_dx1gmshstruct_cfl.440_spreadssf_FNS/'); % alright result
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.3_dx.5gmshstruct_cfl.424_spreadssf/'); % "bad" result, same as OUTPUT_FILES_M.3_dx1gmshstruct_cfl.440_spreadssf
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.3_dx1gmshstruct_cfl.440_spreadssf/'); % "bad" result, same as OUTPUT_FILES_M.3_dx.5gmshstruct_cfl.424_spreadssf
% rootd = strcat(SPCFMEXloc,'validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.3_dx1gmshstruct_cfl.440_spreadssf_FNS/'); % ?? result

% Tests.
% fig_title = 'test';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading.                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OFd = checkOFd(OFd); % Test if OUTPUT_FILES directory exists.
pos_sources = loadSources(OFd); % Load sources' positions.
[xstattab, ystattab, stations_data] = loadStations(OFd); % Load stations data (first try OUTPUT folder, then if not found, try parent DATA folder).
% Compute distance to sources.
dist_to_sources = zeros(size(xstattab, 1), size(pos_sources, 1));
for n_source = 1:size(pos_sources, 1)
  dist_to_sources(:, n_source) = sqrt((xstattab-pos_sources(n_source, 1)) .^ 2 + (ystattab-pos_sources(n_source, 2)) .^ 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ask for user input.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display stations' information.
format longg;
if (convert_to_relative_coords == 1)
  disp(['[',mfilename,'] [istattab, xstattab(istattab), ystattab(istattab), dist_to_sources(istattab)] for all stations (x relative to source, z relative to ground, d relative to source):']);
  %disp([(1:size(pos_stations, 1)); (pos_stations - [pos_sources(1, 1), pos_interface])'; dist_to_sources']);
  disp([(1:size(xstattab, 1)); ([xstattab, ystattab] - [pos_sources(1, 1), pos_interface])'; dist_to_sources']);
else
  disp(['[',mfilename,'] [istattab, xstattab(istattab), ystattab(istattab), dist_to_sources(istattab)] for all stations (absolute x and z, d relative to source):']);
%   disp([(1:size(pos_stations, 1)); pos_stations'; dist_to_sources']);
  disp([(1:size(xstattab, 1))', [xstattab, ystattab], dist_to_sources]);
end
format compact;

% Ask for behaviour.
behaviour = - 1;
while (not(length(behaviour) == 1 && ismember(behaviour, [0, 1, 2, 3])))
  behaviour = input(['[',mfilename,'] Load and plot separately (0), load only (1), load plot time-distance (2), or load and plot polarisation (3)? > ']);
end
% Ask for stations.
istattab = input(['[',mfilename,'] Stations (Matlab format)? > ']); istattab = reshape(istattab,[1,numel(istattab)]);
if(isempty(istattab))
  error(['[',mfilename,', ERROR] Empty vector of stations, retry.']);
end
disp(['[',mfilename,'] Loading [istattab, xstattab(istattab), ystattab(istattab), dist_to_sources(istattab)] (absolute x and z, d relative to source):']);
disp([istattab', xstattab(istattab), ystattab(istattab), dist_to_sources(istattab)]);
nstat = numel(istattab);
% Ask for geometric attenuation (relies on distance to source).
geometric_attenuation = - 1;
while (not(ismember(geometric_attenuation, [0, 1, 2, 3])))
  geometric_attenuation = input(['[',mfilename,'] Apply geometric attenuation factor to data? (0 for no, 1 for d, 2 for |x|, 3 for |z|) > ']);
end
% Ask if plot y-axis should be normalised to same scale.
normalise_ylims = 0; % Default value.
if (behaviour == 0 && nstat > 1)
  normalise_ylims = - 1;
  while (not(normalise_ylims == 0 || normalise_ylims == 1))
    normalise_ylims = input(['[',mfilename,'] Normalise y-scale? (0 for no, 1 for yes) > ']);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and eventually plot.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (convert_to_relative_coords == 1) % Eventually remove source components for display.
  xstattab = xstattab - pos_sources(1, 1); ystattab = ystattab - pos_interface;
end

% Empty arrays.
Ztime = []; Zamp = [];
% Loop on synthetics to be loaded.
for istat = 1:nstat
  istat_glob = istattab(istat); % Recover global number of station.
  factor = getScalings(istat_glob, geometric_attenuation, xstattab, ystattab, dist_to_sources, rescale_factor); % Get scaling factors.
  
  if(ismember(behaviour, [0, 1, 2])) % If direct plots, get the one unknown and proceed.
    [extension, ylabel_unknown] = getUnknowns(type_display, unknown);
    [data, nsamples] = readAndSubsampleSynth(OFd, istat_glob, unknown, extension, subsample, wanted_dt, istat);
    Ztime(istat, 1:nsamples) = data(:, 1)'; Zamp(istat, 1:nsamples) = data(:, 2)'; % Recover time/amplitude data.
    Zamp(istat, :) = factor * Zamp(istat, :); % Scale.
    
  elseif (behaviour == 3) % If polarisation plot, get the two unknowns and proceed.
    [extension_x, ~] = getUnknowns(type_display, 'BXX');
    [data_x] = readAndSubsampleSynth(OFd, istat_glob, 'BXX', extension_x, subsample, wanted_dt, istat);
    [extension_z, ~] = getUnknowns(type_display, 'BXZ');
    [data_z] = readAndSubsampleSynth(OFd, istat_glob, 'BXZ', extension_z, subsample, wanted_dt, istat);
%     sig_t_x = data_x(:, 1)'; sig_v_x = data_x(:, 2)';
    Ztime(istat, :) = data_z(:, 1)'; Zamp(istat,:) = data_z(:, 2)'; % useful to still load and save Z
    Xtime(istat, :) = data_x(:, 1)'; Xamp(istat,:) = data_x(:, 2)'; % useful to still load and save X
    if(not(all(size(data_x(:, 2)') == size(Zamp(istat, :)))))
      error(['size mismatch']);
    end
    data_x(:, 2) = factor * data_x(:, 2); Zamp(istat, :) = factor * Zamp(istat, :); % Scale.
    Xlab='$v_x$';
    Zlab='$v_z$';
    plot_polarisation(Ztime(istat, :), Xamp(istat, :)', Zamp(istat, :), Xlab, Zlab, ['S',num2str(istat_glob),' Polarisation']);
    
    figure(); axx=[];
    subplot(2,1,1); axx=[axx, gca()]; plot(Ztime(istat, :),Zamp(istat,:)); title(['S',num2str(istat_glob),' Time Series']); ylabel(Zlab);
    subplot(2,1,2); axx=[axx, gca()];  plot(Xtime(istat, :),Xamp(istat,:)); ylabel(Xlab);
    linkaxes(axx,'x'); xlabel('time [s]');
  
  else
    error(['[',mfilename,'] behaviour choice not implemented.']);
  end
end
[Ztime, Zamp] = truncToShortest(Ztime, Zamp); % If sythetics don't have the same length, truncate them to shortest.
% Display information.
disp([' ']);
disp(['[',mfilename,'] Data loaded. [matlab_id, istattab, xstattab(istattab), ystattab(istattab), d]:']);
% disp([(1:length(istattab))',istattab', pos_stations(istattab,:),dist_to_sources(istattab)]);
disp([(1:length(istattab))', istattab', xstattab(istattab), ystattab(istattab), dist_to_sources(istattab)]);
disp(strcat("  Example: Data of station ", num2str(istattab(1)), " are in         Zamp(", num2str(1), ", :)."));
disp(strcat("           Corresponding time values are in Ztime(", num2str(1), ", :)."));
disp([' ']);

global synth_load_was_ran
synth_load_was_ran = 1;

if (behaviour == 0) % Eventually, plot one by one in subplots.
  plotOneByOne(Ztime, Zamp, istattab, xstattab(istattab), ystattab(istattab), normalise_ylims, fig_title, ylabel_unknown);
  
elseif (behaviour == 2) % Eventually, plot as time-distance.
  distancechoice = - 1;
  while (~ ismember(distancechoice, [1, 2, 3, 4]))
    distancechoice = input(['[', mfilename, '] Distance choice? (1 for x, 2 for |x|, 3 for z, 4 for d) > ']);
  end
  switch distancechoice
    case 1
      distance = xstattab; dist_symbol = "x"; dist_name = "horizontal distance";
    case 2
      distance = abs(xstattab); dist_symbol = "|x|"; dist_name = "horizontal distance";
    case 3
      distance = ystattab; dist_symbol = "z"; dist_name = "altitude";
    case 4
      distance = dist_to_sources; dist_symbol = "d"; dist_name = "distance";
  end
%   addpath('/home/l.martire/Documents/Ongoing_Work/1811_glanes/treatment_leo'); % wtff?
  reducedtime = - 1;
  while (not(numel(reducedtime)==1 && reducedtime>=0))
    reducedtime = input(['[', mfilename, '] Reduced time? (0 for no, any other value for speed [m/s]) > ']);
  end
  
  plot_time_v_dist(Ztime, Zamp, distance(istattab), reducedtime, fig_title, dist_name);
  
end

if (ismember(plot_amplitude, [1, 2, 3])) % Eventually, plot amplitude as function of distance.
  % Plot amplitude.
  amp = zeros(size(Zamp));
  for ii = 1:length(Zamp(:, 1))
    amp(ii) = max(Zamp(ii, :)) - min(Zamp(ii, :));
  end
  switch(plot_amplitude)
    case 1
      [sorted_dist, isort] = sort(xstattab(istattab));
    case 2
      [sorted_dist, isort] = sort(ystattab(istattab));
    case 3
      [sorted_dist, isort] = sort(dist_to_sources(istattab));
  end
  sorted_amp = amp(isort);
  figure();
  loglog(sorted_dist, sorted_amp);
end


% StratoBaro, 66, June, 12:00
% fig_title = strcat('Microbaroms, lat66, June, 12:00');
% rootd = strcat(SPCFMEXloc,'microbaroms_patch'); OFd = strcat(rootd, '/OUTPUT_FILES_668482_disp7_isp6_full/');
% rootd = strcat(SPCFMEXloc,'microbaroms_patch'); OFd = strcat(rootd, '/OUTPUT_FILES_668482_disp7/');
% rootd = strcat(SPCFMEXloc,'microbaroms_patch'); OFd = strcat(rootd, '/OUTPUT_FILES_668446_disp7_wrongstations/');
% rootd = strcat(SPCFMEXloc,'microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_668354_testlarger_str_1e-1mps_isp6/');
% rootd = strcat(SPCFMEXloc,'microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_656744_straight_1mps_isp6/');
% rootd = strcat(SPCFMEXloc,'microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_656505_straight_1e-2mps_isp6/');
% rootd = strcat(SPCFMEXloc,'microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_656465_straight_1mps_test/');
% rootd = strcat(SPCFMEXloc,'microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_655513_analytic/');
% rootd = strcat(SPCFMEXloc,'microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_655494_alright/');
% rootd = strcat(SPCFMEXloc,'microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_655487/');
% rootd = strcat(SPCFMEXloc,'microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_655369_unstable/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_650851/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_644923_EBF_ispread3/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_641616_EBF_ispread1.5/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_641395_testEBF_stopped/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_639014_long_betterEBF_crash/');
% rootd = strcat(SPCFMEXloc,'stratobaro_test_EBF/'); OFd = strcat(rootd, '/OUTPUT_FILES/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'stratobaro_test_EBF/'); OFd = strcat(rootd, '/OUTPUT_FILES_test1/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_637450_long_EBF_crash/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_634307_testexternalforcing/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_624650_long/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_624515_rpw_spatially_fixed_s0.2/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_624478_apo+rpw0.2_10.5p/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_624436_apo+rpw0.15/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_623945_apodised/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_622147_final_test/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_622125/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_622037/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_621952/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_621860/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_621802/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_599638_testmicrobarom/');

% StratoExplo, 66, June, 12:00
% fig_title = strcat('Atmospheric Explosions, lat66, June, 12:00');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_597316/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_597250/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_597099/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_595500_crash40k1it/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRdfsdfATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_594736_crash27kit/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_594361_dt1e-3_cancelled/');

% Seismic Hammer, soft soil.
% fig_title = strcat('Seismic Hammer Simulation (Soft Soil)'); coord_units = 'm'; convert_to_relative_coords = 1; pos_interface = 308;
% rootd = strcat(SPCFMEXloc,'SH_final'); OFd = strcat(rootd, '/OUTPUT_FILES_669168_fullretweaked/'); rescale_factor = 236; % Same as 593959 but test with first layers changed.
% rootd = strcat(SPCFMEXloc,'SH_final'); OFd = strcat(rootd, '/OUTPUT_FILES_668888_stopped_12kit/'); rescale_factor = 236; % Same as 593959 but test with first layers changed.
% rootd = strcat(SPCFMEXloc,'SH_axisym'); OFd = strcat(rootd, '/OUTPUT_FILES_660223_full_dec1m/'); % Same as 593959 but axisymmetric.
% rootd = strcat(SPCFMEXloc,'SH_final'); OFd = strcat(rootd, '/OUTPUT_FILES_627577_qk4sls_truefreesurf/');
% rootd = strcat(SPCFMEXloc,'SH_final'); OFd = strcat(rootd, '/OUTPUT_FILES_623195_qk_4sls_freesurf/');
% rootd = strcat(SPCFMEXloc,'SH_final/SH_soft_final_redone'); OFd = strcat(rootd, '/OUTPUT_FILES_610770/'); % With additionnal stations for comparison with data.
% rootd = strcat(SPCFMEXloc,'SH_final'); OFd = strcat(rootd, '/SH_soft_final_redone_Qkappa_616368/'); % Same as 593959 only with Qp converted to Qk and additionnal stations.
% rootd = strcat(SPCFMEXloc,'SH_final'); OFd = strcat(rootd, '/SH_soft_final_redone_Qkappa+f = 0_618645/');
% rootd = strcat(SPCFMEXloc,'SH_final'); OFd = strcat(rootd, '/SH_soft_final_redone_Qkappa+f = f0_618882/');
% rootd = strcat(SPCFMEXloc,'SH_final'); OFd = strcat(rootd, '/SH_soft_final_redone_qk_noatt_619264');
% rootd = strcat(SPCFMEXloc,'SH_final'); OFd = strcat(rootd, '/SH_soft_final_redone_qk_att4sls_620294');
% rootd = strcat(SPCFMEXloc,'SH_final/SH_soft_final'); OFd = strcat(rootd, '/OUTPUT_FILES_593959/'); % Original (used in paper).
% rootd = strcat(SPCFMEXloc,'ON_EOS__seismic_hammer_new_model'); OFd = strcat(rootd, '/OUTPUT_FILES_551980_seismic_potential_with_memvars_solid/');

% Seismic Hammer, hard soil.
% fig_title = strcat('Seismic Hammer Simulation (Hard Soil)'); coord_units = 'm'; convert_to_relative_coords = 1; pos_interface = 308;
% rootd = strcat(SPCFMEXloc,'SH_final/SH_hard_final'); OFd = strcat(rootd, '/OUTPUT_FILES_593960/'); % Original (used in paper).
% rootd = strcat(SPCFMEXloc,'SH_hard_axisym'); OFd = strcat(rootd, '/OUTPUT_FILES_661601_full_dec1m/'); % Same as 593960 but axisymmetric.
% rootd = strcat(SPCFMEXloc,'SH_hard_axisym'); OFd = strcat(rootd, '/OUTPUT_FILES_661609_full_onlypress/'); type_display = 4; unknown = 'PRE'; % Same as 661601 but only recording above ground.

% Quake, 45.
% fig_title = strcat('Quake Simulation (45$^\circ$ dip)');
% rootd = strcat(SPCFMEXloc,'OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_1811221612_local');
% rootd = strcat(SPCFMEXloc,'OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_1811221544_local');
% rootd = strcat(SPCFMEXloc,'test_impedance'); OFd = strcat(rootd, '/OUTPUT_FILES');
% rootd = strcat(SPCFMEXloc,'OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_75118_isoth_d6_savedvdg');
% rootd = strcat(SPCFMEXloc,'OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_71984_isoth_d6');
% rootd = strcat(SPCFMEXloc,'OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_71980_isothermal');
% rootd = strcat(SPCFMEXloc,'OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_pot');
% rootd = strcat(SPCFMEXloc,'OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_71931_d9');
% rootd = strcat(SPCFMEXloc,'OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_71920_force_instead_of_moment_d9');
% rootd = strcat(SPCFMEXloc,'OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_71913');
% rootd = strcat(SPCFMEXloc,'OKQ'); OFd = strcat(rootd, '/OUTPUT_FILES_668844_OKQ45_redone'); rescale_factor = 1e-3;
% rootd = strcat(SPCFMEXloc,'OKQ/ON_EOS_quake_ok_45'); OFd = strcat(rootd, '/OUTPUT_FILES_583041_long');

% Quake, 0.
% fig_title = strcat('Quake Simulation (0$^\circ$ dip)');
% rootd = strcat(SPCFMEXloc,'OKQ'); OFd = strcat(rootd, '/OUTPUT_FILES_668833_OKQ0_redone'); rescale_factor = 1e-3;
% rootd = strcat(SPCFMEXloc,'OKQ/ON_EOS_quake_ok_0'); OFd = strcat(rootd, '/OUTPUT_FILES_586984_full');

% Seismic Hammer, hard soil.
% fig_title = strcat('Seismic Hammer Simulation (Hard Soil)'); coord_units = 'm'; convert_to_relative_coords = 1;
% rootd = strcat(SPCFMEXloc,'SH_final/SH_hard_final'); OFd = strcat(rootd, '/OUTPUT_FILES_593960/');
% rootd = strcat(SPCFMEXloc,'ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580457_full/'); rescale_factor = 8.840811261618920e-04;
% rootd = strcat(SPCFMEXloc,'ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580113/');
% rootd = strcat(SPCFMEXloc,'ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580185/');
% rootd = strcat(SPCFMEXloc,'ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580228/');
% rootd = strcat(SPCFMEXloc,'ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580333/');
% rootd = strcat(SPCFMEXloc,'ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580712/');

% Seismic Hammer, soft soil.
% fig_title = strcat('Seismic Hammer Simulation (Soft Soil)'); coord_units = 'm'; convert_to_relative_coords = 1;
% rootd = strcat(SPCFMEXloc,'SH_final/SH_soft_final'); OFd = strcat(rootd, '/OUTPUT_FILES_593959/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9113508_seismic_DG_with_memvars_solid/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9048100_seismic_DG/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9081476_seismic_potential/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9091088_seismic_DG_new_coupling/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9102702_seismic_potential_rem_forcing/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_551980_seismic_potential_with_memvars_solid/');

% Tests.
% fig_title = 'test';
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf4/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf10_1hz/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf2_1hz/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf5_jpeguz/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf31hz_homo_otherstation/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf31hz_homo/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf4_homogenous/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_sharpstf3/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf5/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf3/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf2/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'test_external_modelDG_only'); OFd = strcat(rootd, '/OUTPUT_FILES_1e0mu/');
% rootd = strcat(SPCFMEXloc,'test_external_modelDG_only'); OFd = strcat(rootd, '/OUTPUT_FILES_1e2mu/');
% rootd = strcat(SPCFMEXloc,'test_external_modelDG_only'); OFd = strcat(rootd, '/OUTPUT_FILES_1e4mu/');
% rootd = strcat(SPCFMEXloc,'test_external_modelDG_only'); OFd = strcat(rootd, '/OUTPUT_FILES_1e5mu/');
% rootd = strcat(SPCFMEXloc,'full_DG_square'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'test_stretching'); OFd = strcat(rootd, '/OUTPUT_FILES_long/');
% rootd = strcat(SPCFMEXloc,'test_FTS'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'test_coupling'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'test_stretching'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'test_stretching_wind'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'test_stretching_FFcounterpart'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'ON_EOS_test_atmo'); OFd = strcat(rootd, '/OUTPUT_FILES_TEST');
% rootd = strcat(SPCFMEXloc,'ON_EOS_test_densitysource'); OFd = strcat(rootd, '/OUTPUT_FILES_TEST');
% rootd = strcat(SPCFMEXloc,'test_quake'); OFd = strcat(rootd, '/OUTPUT_FILES');
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES__long_working');
% rootd = strcat(SPCFMEXloc,'quake_oklahoma0'); OFd = strcat(rootd, '/OUTPUT_FILES_test_vz');
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_test_vz');
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_test_dz'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_ok_45'); OFd = strcat(rootd, '/OUTPUT_FILES_narrow_okdx');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_583123_100km_3sources_nospread');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_583128_100km_3sources_nospread');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_583138_100km_3sources_spread');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_test_scale_sources');
% rootd = strcat(SPCFMEXloc,'ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_591778_66j1200_regmukap0_softground_nocrash');
% rootd = strcat(SPCFMEXloc,'seismic_hammer_zooms/soft/'); OFd = strcat(rootd, 'OUTPUT_FILES_583180');
% rootd = strcat(SPCFMEXloc,'seismic_hammer_zooms/hard/'); OFd = strcat(rootd, 'OUTPUT_FILES_583194');
% rootd = strcat(SPCFMEXloc,'seismic_hammer_zooms/hard/'); OFd = strcat(rootd, 'OUTPUT_FILES_586795_d4_source_flipped'); rescale_factor = 8.840811261618920e-04;
% rootd = strcat(SPCFMEXloc,'seismic_hammer_zooms/soft/'); OFd = strcat(rootd, 'OUTPUT_FILES_591776_additional_stations');fig_title = 'test_soft';
% rootd = strcat(SPCFMEXloc,'seismic_hammer_zooms/hard/'); OFd = strcat(rootd, 'OUTPUT_FILES_591777_additional_stations');fig_title = 'test_hard';

% Mars Gravity Wave.
% fig_title = strcat('Mars Gravity Wave Simulation');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_533937_vNEW_full/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_534758_long_instab/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_535011_with_FTS/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_535489_removed_discontinuity_long/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave/test_RAPHAEL'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_540064_FTS_no_disc_long/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_9078210_spread_source/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_9081352_spread_cut_source/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_9091089_new_coupling/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_9103256_same_as_previous_but_factor_1/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_552471_atmo_only/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_552455_factor0p1/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_557219_long/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_558183_Gderiv/');
% rootd = strcat(SPCFMloc, 'Ongoing_Work/SPECFEM-DG_Mars_AGW_runs/explo_mars_sub'); OFd = strcat(rootd, '/OUTPUT_FILES_KappaON/');

% DAG
% fig_title = strcat('DAG');
% rootd = strcat(SPCFMEXloc,'DAG_testlocal/'); OFd = strcat(rootd, 'OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'DAG_testlocal/'); OFd = strcat(rootd, 'OUTPUT_FILES_testlocal/');
% rootd = strcat(SPCFMEXloc,'DAG/'); OFd = strcat(rootd, 'OUTPUT_FILES_1303372/');
% rootd = strcat(SPCFMEXloc,'DAG/'); OFd = strcat(rootd, 'OUTPUT_FILES_103088_uglybutnicerefrac/');

% Tests
% fig_title = strcat('test');
% rootd = strcat(SPCFMEXloc,'test_pml'); OFd = strcat(rootd, '/OUTPUT_FILES_d = 0_kmax = 2/');
% rootd = strcat(SPCFMEXloc,'test_pml'); OFd = strcat(rootd, '/OUTPUT_FILES_d = d_kmax = 2/');
% rootd = strcat(SPCFMEXloc,'test_pml'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'test_plot_perio'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'tir_mars'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'demo_pot'); OFd = strcat(rootd, '/OUTPUT_FILES_826234/');
% rootd = strcat(SPCFMEXloc,'demo_fns'); OFd = strcat(rootd, '/OUTPUT_FILES_826226/');
% rootd = strcat(SPCFMEXloc,'demo_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_826213/');
% rootd = strcat(SPCFMEXloc,'demo_fns'); OFd = strcat(rootd, '/OUTPUT_FILES_fnsf2s_local/');
% rootd = strcat(SPCFMEXloc,'demo_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_lnsf2s_local/');
% rootd = strcat(SPCFMEXloc,'demo_fns'); OFd = strcat(rootd, '/OUTPUT_FILES_fnsf2s_local_butd7/');
% rootd = strcat(SPCFMEXloc,'demo_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_lnsf2s_local_butd7/');
% rootd = strcat(SPCFMEXloc,'demo_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_lns_t = 105s/');
% rootd = strcat(SPCFMEXloc,'demo_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_fns_t = 195s/');
% rootd = strcat(SPCFMEXloc,'mb_gmsh'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'test_impedance'); OFd = strcat(rootd, '/OUTPUT_FILES_extatm+oksoil+lowdt_BUTONCALMIP/');
% rootd = strcat(SPCFMEXloc,'test_impedance'); OFd = strcat(rootd, '/OUTPUT_FILES_extatm+oksoil+lowdt/');
% rootd = strcat(SPCFMEXloc,'test_impedance'); OFd = strcat(rootd, '/OUTPUT_FILES_lns/');
% rootd = strcat(SPCFMEXloc,'test_lns'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'test_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_fts_lns_19s/');
% rootd = strcat(SPCFMEXloc,'test_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_fts_fns_26s/');
% rootd = strcat(SPCFMEXloc,'test_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_f_lns_44s/');
% rootd = strcat(SPCFMEXloc,'test_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_f_fns_62s/');
% rootd = strcat(SPCFMEXloc,'test_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_stf_lns_27s/');
% rootd = strcat(SPCFMEXloc,'test_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_stf_fns_45s/');
% rootd = strcat(SPCFMEXloc,'axisym_test'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'test_external_forcing'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd = strcat(SPCFMEXloc,'test_stretching_wind'); OFd = strcat(rootd, '/OUTPUT_FILES/'); coord_units = 'm'; convert_to_relative_coords = 0; pos_interface = 0;
% rootd = strcat(SPCFMEXloc,'test_stretching_wind'); OFd = strcat(rootd, '/OUTPUT_FILES_observesignalinbuffer_cstrhdr/'); coord_units = 'm'; convert_to_relative_coords = 0; pos_interface = 0;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf4/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf10_1hz/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf2_1hz/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf5_jpeguz/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf31hz_homo_otherstation/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf31hz_homo/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf4_homogenous/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_sharpstf3/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf5/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf3/'); type_display = 1;
% rootd = strcat(SPCFMEXloc,'quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf2/'); type_display = 1;

function dir = checkOFd(dir)
  if (not(strcmp(dir(end), '/')))
    dir = [dir, '/'];
  end
  if (not(exist(dir, 'dir')))
    error(['[',mfilename,', ERROR] OUTPUT_FILES directory does not exist (', dir, ').']);
  end
end

function pos_sources = loadSources(OFdir)
  pos_sources = [inf, inf]; % Allocate a row for the first source's position.
  % fid = fopen([rootd, '/DATA/SOURCE']);
  fid = fopen([OFdir, 'SOURCE']);
  if (fid == - 1)
    fid = fopen([OFdir, 'input_source']);
    if (fid == - 1)
      error(['[',mfilename,', ERROR] Cannot open SOURCE file (', OFdir, 'SOURCE', ').']);
    end
  end
  line = 0; xfound = 0; zfound = 0; stop = 0;
  while (stop == 0)
    % TODO: Loop on source number.
    line = fgetl(fid);
    if length(line) > 0
      if (line == - 1)
        stop = 1;
      end
      line = regexprep(regexprep(line, ' +', ' '), '^ ', ''); % Remove multiple spaces, and then eventually remove space if it there is one as first character.
      if strcmp(line(1:2), 'xs')
        xfound = 1; pos_sources(1, 1) = str2num(regexprep(regexprep(line(3:end), ' *#.*', ''), ' * = * *', '')); % Remove comments (everything after a '#'), remove the equals sign and spaces around it, and cast it as source position.
      end
      if strcmp(line(1:2), 'zs')
        zfound = 1; pos_sources(1, 2) = str2num(regexprep(regexprep(line(3:end), ' *#.*', ''), ' * = * *', ''));
      end
    end
    if (xfound && zfound)
      break
    end
  end
  fclose('all');
end

function [ext, unknown] = getUnknowns(type_displ, unknown)
  % Switch on type of display.
  switch type_displ
%   if (type_display == 1) % Original SPECFEM2D's synthetic is displacement.
    case 1 % Original SPECFEM2D's synthetic is displacement.
      ext = 'semd'; % Because original SPECFEM2D's synthetic is displacement.
      % For stations in solid zones it's displacement. For stations in DG zones it's velocity.
      if (strcmp(unknown, 'BXZ'))
        unknown = 'vertical {$u_z$ (m), $v_z$ [m/s]}';
      elseif (strcmp(unknown, 'BXX'))
        unknown = 'horizontal {$u_x$ (m), $v_x$ [m/s]}';
      else
        error(['[',mfilename,', ERROR] The variable ''unknown'' has a non-standard value.']);
      end
%   elseif (type_display == 2) % Original SPECFEM2D's synthetic is velocity.
    case 2 % Original SPECFEM2D's synthetic is velocity.
      ext = 'semv'; % Because original SPECFEM2D's synthetic is velocity.
      if (strcmp(unknown, 'BXZ'))
        unknown = '{$v_z$ [m/s], $\delta P$ [Pa]}';
      elseif (strcmp(unknown, 'BXX'))
        unknown = '{$v_x$ [m/s], $\delta P$ [Pa]}';
      else
        error(['[',mfilename,', ERROR] The variable ''unknown'' has a non-standard value.']);
      end
%   elseif (type_display == 4) % Original SPECFEM2D's synthetic is pressure.
    case 4 % Original SPECFEM2D's synthetic is pressure.
      ext = 'semp'; % Because original SPECFEM2D's synthetic is pressure.
      if (strcmp(unknown, 'PRE'))
        unknown = '$\delta P$ [Pa]';
      else
        error(['[',mfilename,', ERROR] The variable ''unknown'' has a non-standard value.']);
      end
    otherwise
      error(['[',mfilename,', ERROR] This type_display (',num2str(type_displ),') is not implemented for this script.']);
  end
end

function factor = getScalings(stat_number, geomAtt, x_stat, z_stat, d_stat, rescale_fact)
  % Renormalisation (global, and geometric).
  factor = 1; % Reset to default value for each station.
  if (geomAtt ~= 0)
    % If geometric_attenuation was asked by user.
    switch geomAtt
      case 1
        geom_att_fact = d_stat(stat_number) ^ 0.5; % Geometric attenuation respective to raw distance to source.
      case 2
        geom_att_fact = abs(x_stat(stat_number)) ^ 0.5; % Geometric attenuation respective to horizontal distance to source.
      case 3
        geom_att_fact = abs(z_stat(stat_number)) ^ 0.5; % Geometric attenuation respective to vertical distance to source.
      otherwise
        error(['[',mfilename,', ERROR] Geometric attenuation parameter not implemented.']);
    end
    if (geom_att_fact == 0)
      % If exactly at zero distance, consider no geometric rescaling.
      geom_att_fact = 1;
    end
    factor = factor / geom_att_fact;
  end
  renorm_statbystat = - 1;
  if (rescale_fact ~= 1)
    % Rescaling was asked. Check again with user.
    renorm_statbystat = - 1;
    disp(['[',mfilename,'] Specified rescale factor is ', num2str(rescale_fact), '.']);
    inputtxt = ['[',mfilename,'] Rescale data for station ', num2str(stat_number), '? (0 for no, 1 for yes) > '];
    while (not(ismember(renorm_statbystat, [0, 1])))
      renorm_statbystat = input(inputtxt);
      if (isempty(renorm_statbystat))
        renorm_statbystat = 1;
      end
    end
  end
  if (renorm_statbystat == 1)
    % If rescaling is actually wanted by user, do it.
    factor = factor * rescale_fact;
  end
end

function plotOneByOne(Ztime, Zamp, istattab, xstat, zstat, norm_ylims, figtitle, ylab)
  if(not(all(size(Ztime)==size(Zamp))))
    error('must have same sizes');
  end
  if(not(exist('figtitle')))
    figtitle='';
  end
  if(not(exist('ylab')))
    ylab='SIGNAL';
  end
  nstat = size(Ztime, 1);
  figure();
  hold on;
  axes = {};
  for istat = 1:nstat
    istat_glob = istattab(istat); % Recover global number of station.
    subplot(nstat, 1, istat);
    axes{istat} = gca;
    legtext{istat} = ['S', num2str(istat_glob), ', $(x,z) = (', num2str(xstat(istat)), ',', num2str(zstat(istat)), ')$ [m]'];
    plot(Ztime(istat, :), Zamp(istat, :), 'displayname', legtext{istat}); hold on;
    % Cosmetics.
    if (istat == 1); title(figtitle); end; % Put title only on first subplot.
    if (istat == nstat); xlabel('time [s]'); end; % Put xlabel only on last subplot.
    if (istat ~= nstat); set(gca, 'xticklabel', []); end; % Remove xticks for all subplots except last one.
    if (istat == round(nstat / 2)); ylabel(ylab); end; % Put one ylabel at the middle subplot.
    xlim([min(min(Ztime)), max(max(Ztime))]);
    legend('Location', 'northeast');
    set(gca, 'TickLabelInterpreter', 'latex'); grid on; box on;
  end
  axess=[];
  for i=1:numel(axes); axess(i)=axes{i}; end;
  if (norm_ylims)
    linkaxes(axess); % Link both x and y.
  else
    linkaxes(axess, 'x'); % Link only x.
  end
end