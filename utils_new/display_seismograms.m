% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

clear all;
% close all
clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 8); % Default at 6.
set(0, 'defaultTextFontSize', 12);
set(0, 'defaultAxesFontSize', 12); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

%%%%%%%%%%%%%%%%%
% INITIALIZATION
addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/Atmospheric_Models');
renorm_factor=1; SPCFMloc='/home/l.martire/Documents/SPECFEM/';
% Direct waves parameters (??)
%windsurf = - 10.5; % m/s
%csurf = 233.5; % m/s
%tstart = 9.0;
% window size for Amplitude estimate
%tws = 20.0; % s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% StratoExplo, 66, June, 12:00
fig_title = strcat('Stratospheric Explosions, lat66, June, 12:00');
rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_594361_dt1e-3_cancelled/');

% Mars AGW.
% rootd=strcat(SPCFMloc, 'Ongoing_Work/SPECFEM-DG_Mars_AGW_runs/explo_mars_sub'); OFd = strcat(rootd, '/OUTPUT_FILES_KappaON/');

% Mars Gravity Wave.
% fig_title = strcat('Mars Gravity Wave Simulation');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_533937_vNEW_full/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_534758_long_instab/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_535011_with_FTS/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_535489_removed_discontinuity_long/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave/test_RAPHAEL'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_540064_FTS_no_disc_long/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_9078210_spread_source/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_9081352_spread_cut_source/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_9091089_new_coupling/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_9103256_same_as_previous_but_factor_1/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_552471_atmo_only/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_552455_factor0p1/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_557219_long/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_558183_Gderiv/');

% Seismic Hammer, soft soil.
% fig_title = strcat('Seismic Hammer Simulation (Soft Soil)');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_soft_final'); OFd = strcat(rootd, '/OUTPUT_FILES_593959/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9113508_seismic_DG_with_memvars_solid/');

% rootd=strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9048100_seismic_DG/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9081476_seismic_potential/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9091088_seismic_DG_new_coupling/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9102702_seismic_potential_rem_forcing/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_551980_seismic_potential_with_memvars_solid/');

% Seismic Hammer, hard soil.
% fig_title = strcat('Seismic Hammer Simulation (Hard Soil)');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_hard_final'); OFd = strcat(rootd, '/OUTPUT_FILES_593960/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580457_full/'); renorm_factor=8.840811261618920e-04;

% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580113/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580185/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580228/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580333/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580712/');

% Quake, 45.
% fig_title = strcat('Quake Simulation (45d dip)');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_quake_ok_45'); OFd = strcat(rootd, '/OUTPUT_FILES_583041_long');

% Quake, 0.
% fig_title = strcat('Quake Simulation (0d dip)');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_quake_ok_0'); OFd = strcat(rootd, '/OUTPUT_FILES_586984_full');

% Tests.
% fig_title = 'test';
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/full_DG_square'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_stretching'); OFd = strcat(rootd, '/OUTPUT_FILES_long/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_FTS'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_coupling'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_stretching'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_stretching_wind'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_stretching_FFcounterpart'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_test_atmo'); OFd = strcat(rootd, '/OUTPUT_FILES_TEST');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_test_densitysource'); OFd = strcat(rootd, '/OUTPUT_FILES_TEST');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_quake'); OFd = strcat(rootd, '/OUTPUT_FILES');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES__long_working');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma0'); OFd = strcat(rootd, '/OUTPUT_FILES_test_vz');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_test_vz');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_ok_45'); OFd = strcat(rootd, '/OUTPUT_FILES_narrow_okdx');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_583123_100km_3sources_nospread');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_583128_100km_3sources_nospread');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_583138_100km_3sources_spread');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_test_scale_sources');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_591778_66j1200_regmukap0_softground_nocrash');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/seismic_hammer_zooms/soft/'); OFd = strcat(rootd, 'OUTPUT_FILES_583180');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/seismic_hammer_zooms/hard/'); OFd = strcat(rootd, 'OUTPUT_FILES_583194');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/seismic_hammer_zooms/hard/'); OFd = strcat(rootd, 'OUTPUT_FILES_586795_d4_source_flipped'); renorm_factor=8.840811261618920e-04;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/seismic_hammer_zooms/soft/'); OFd = strcat(rootd, 'OUTPUT_FILES_591776_additional_stations');fig_title = 'test_soft';
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/seismic_hammer_zooms/hard/'); OFd = strcat(rootd, 'OUTPUT_FILES_591777_additional_stations');fig_title = 'test_hard';

% Sub-sample? Useful for lengthy seismograms.
subsample = 0;

% Quantity to display:
%   1 = displacement for non-DG and velocity for DG,
%   2 = velocity for non-DG and pressure perturbation (Pa) for DG.
type_display = 2; % Should be the same as the seismotype variable in parfile.

% Unknown:
% For type_display==2 and stations in DG zones, pressure perturbation (Pa) is saved both in BXX and BXZ files.
% unknown = 'BXX';
unknown = 'BXZ';

% close all; % Close all figure. Comment this to keep them.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
renorm=0; % Do not renormalise by default.
if(not(strcmp(OFd(end),'/'))); OFd=[OFd,'/']; end;
if(not(exist(OFd))); error(strcat("OUTPUT_FILES directory does not exist (",OFd,').')); end;
% Load sources' positions.
pos_sources = [inf, inf]; % Allocate a row for the first source's position.
% fid = fopen([rootd, '/DATA/SOURCE']);
fid = fopen([OFd, 'SOURCE']);
if(fid==-1)
  error(strcat("Cannot open SOURCE file (",[OFd, 'SOURCE'],').'));
end
line = 0; xfound = 0; zfound = 0; stop=0;
while(stop==0)
  % TODO: Loop on source number.
  line = fgetl(fid);
  if length(line)>0
    if(line==-1)
      stop=1;
    end
    line = regexprep(regexprep(line, ' +', ' '),'^ ',''); % Remove multiple spaces, and then eventually remove space if it there is one as first character.
    if strcmp(line(1:2), 'xs')
      xfound = 1; pos_sources(1, 1) = str2num(regexprep(regexprep(line(3:end), ' *#.*', ''), ' *=* *','')); % Remove comments (everything after a '#'), remove the equals sign and spaces around it, and cast it as source position.
    end
    if strcmp(line(1:2), 'zs')
      zfound = 1; pos_sources(1, 2) = str2num(regexprep(regexprep(line(3:end), ' *#.*', ''), ' *=* *',''));
    end
  end
  if(xfound && zfound)
    break
  end
end
fclose('all');

% Load stations data (first try OUTPUT folder, then if not found, try parent DATA folder).
try
  A = importdata(strcat(OFd, 'STATIONS'));
catch
  try
    A = importdata(strcat(rootd, '/DATA/STATIONS'));
  catch
    error('Cannot find STATIONS file.');
  end
end
pos_stations = [A.data(:, 1) A.data(:, 2)];
xstattab = pos_stations(:, 1);
ystattab = pos_stations(:, 2);
% Compute distance to sources.
dist_to_sources = zeros(size(pos_stations, 1), size(pos_sources, 1));
for n_source = 1:size(pos_sources, 1)
  dist_to_sources(:, n_source) = sqrt((pos_stations(:, 1) - pos_sources(n_source, 1)) .^ 2 + (pos_stations(:, 2) - pos_sources(n_source, 2)) .^ 2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stations to display.        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display stations' informations.
format shortG; disp("  [station_id; x; z; d] for all stations:"); disp([(1:size(pos_stations, 1)); pos_stations';dist_to_sources']); format compact;

% Ask user for various inputs.
display_or_load=-1;
while(not(display_or_load==0 || display_or_load==1)); display_or_load=input('  Load and display (0) or load only (1)?\n  > '); end
istattab = input(['  Stations (Matlab format, e.g. [1, 4, 7] or 1:20)?\n  > ']);
nstat = size(pos_stations(istattab, 1), 1);
geometric_attenuation=-1; inputtxt=char(strcat("  Apply geometric attenuation (1/sqrt(d) factor) to data? (0 for no, 1 for yes)\n  > "));
while(not(geometric_attenuation==0 || geometric_attenuation==1)); geometric_attenuation=input(inputtxt); end
if(display_or_load==0 && nstat>1)
  normalise_ylims=-1;
  while(not(normalise_ylims==0 || normalise_ylims==1)); normalise_ylims=input('  Normalise y-scale? (0 for no, 1 for yes)\n  > '); end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ??                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scale of recordsection (??)
% lbar = 10;
% if (type_display == 1)
%   recscale = 50;
% else
%   % recscale = 300.0;
%   recscale = 1;
% end

% reduction velocity (??)
% redvel = 250; % m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and eventually plot.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_or_load==0)
  figure(); hold on;
end

% Loop on sismograms.
max_ylim_plus=-Inf;
min_ylim_minus=+Inf;
for istat = 1 : nstat
  istat_glob = istattab(istat); % Recover global numer of station.

  % Unused?
  %   if (istat_glob > 9)
  %     cmplt = '000';
  %   else
  %     cmplt = '0000';
  %   end

  % Switch on type of display.
  if (type_display == 1)
    % Original SPECFEM2D's sismogram is displacement.
    extension = "semd"; % Because original SPECFEM2D's sismogram is displacement.
    % For stations in solid zones it's displacement. For stations in DG zones it's velocity.
    if(strcmp(unknown,'BXZ'))
      unknown_name = 'vertical {$u_z$ (m), $v_z$ (m/s)}';
    elseif(strcmp(unknown,'BXX'))
      unknown_name = 'horizontal {$u_x$ (m), $v_x$ (m/s)}';
    else
      error("The variable 'unknown' has a non-standard value.");
    end
  elseif (type_display == 2)
    % Original SPECFEM2D's sismogram is velocity.
    extension = "semv"; % Because original SPECFEM2D's sismogram is velocity.
    if(strcmp(unknown,'BXZ'))
      unknown_name = '{$v_z$ (m/s), $\delta P$ (Pa)}';
    elseif(strcmp(unknown,'BXX'))
      unknown_name = '{$v_x$ (m/s), $\delta P$ (Pa)}';
    else
      error("The variable 'unknown' has a non-standard value.");
    end
  end

  % Read the sismogram.
  file = strcat(OFd, 'AA.', A.textdata(istat_glob, 1), '.', unknown, '.', extension);
  data = load(file{1});
  nt = max(size(data));
  if(subsample==1)
    % Sub-sample of records.
    nsub = ceil(nt/1000);
    nd = max(size(data(1:nsub:nt, 1)));
  else
    nsub=1;
    nd=nt;
  end
  % Recover time/amplitude data.
  Ztime(istat, 1:nd) = data(1:nsub:nt, 1)';
  %Ztime(istat, 1:nd) = Ztime(istat, 1:nd) - Ztime(istat, 1); % Make time values start at zero.
  % Correct for sign of source function for N wave pattern (??).
  Zamp(istat, 1:nd) = data(1:nsub:nt, 2)';

  % Display.
  if(display_or_load==0) % If display.
    ax(istat) = subplot(nstat, 1, istat);
    legtext=strcat('S', num2str(istat_glob), ', (x,z,d)=(', num2str(xstattab(istat_glob) / 1000), ',', num2str(ystattab(istat_glob) / 1000), ',', num2str(dist_to_sources(istat_glob) / 1000), ') km');
    
    factor=1;
    if(geometric_attenuation==1)
      factor=factor/(dist_to_sources(istat_glob)^0.5);
    end
    if(renorm_factor~=1)
      renorm=-1; disp(strcat("    Specified renormalisation factor is ", num2str(renorm_factor), ".")); inputtxt=char(strcat("    Renormalise data for station ", legtext,"? (0 for no, 1 for yes) > "));
      while(not(renorm==0 || renorm==1)); renorm=input(inputtxt); end;
    end
    if(renorm==1)
      factor=factor*renorm_factor;
    end
    plot(Ztime(1, :), factor*Zamp(istat, :));

    % Cosmetics.
    if (istat == 1)
      title(fig_title)
    end
    if (istat == nstat)
      xlabel('time (s)')
    end
    if (istat ~= nstat)
      set(gca, 'xticklabel', []);
    end
    if (istat == round(nstat / 2))
      ylabel(unknown_name);
    end
    xlim([Ztime(1, 1), Ztime(1, end)]);

    legend(legtext, 'Location', 'northeast');

    hold on;

    ax=gca;
    if(ax.YLim(1)<min_ylim_minus)
      min_ylim_minus=ax.YLim(1);
    end
    if(ax.YLim(2)>max_ylim_plus)
      max_ylim_plus=ax.YLim(2);
    end
  else % Load only.
    % Nothing to do.
  end
end
if(display_or_load==0)
  linkaxes(ax, 'x');
  %tightfig;

  if(nstat>1)
    if(normalise_ylims)
      f=gcf;
      for i=1:length(f.Children)
        if(strcmp(f.Children(i).Type,'axes'))
          f.Children(i).YLim=[min_ylim_minus, max_ylim_plus];
        end
      end
    end
  end
  f=gcf; figure(f.Number);
end
disp("  Data loaded. [matlab_id, station_id]:");
disp([(1:length(istattab))',istattab']);
disp(strcat("  Example: Data of station ",num2str(istattab(1))," are in         Zamp(",num2str(1),", :)."));
disp(strcat("           Corresponding time values are in Ztime(",num2str(1),", :)."));