% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

clear all;
clc;
format compact;
set(0, 'DefaultLineLineWidth', 1.5); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 6); % Default at 6.
set(0, 'defaultTextFontSize', 16);
set(0, 'defaultAxesFontSize', 14); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

% scale of recordsection

%%%%%%%%%%%%%%%%%
% INITIALIZATION

% Direct waves parameters (??)
%windsurf = - 10.5; % m/s
%csurf = 233.5; % m/s
%tstart = 9.0;
% window size for Amplitude estimate
%tws = 20.0; % s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mars AGW.
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/SPECFEM-DG_Mars_AGW_runs/explo_mars_sub'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_KappaON/');

% Mars Gravity Wave.
% fig_title = strcat('Mars Gravity Wave Simulation');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Mars_Gravity_Wave'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_533937_vNEW_full/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Mars_Gravity_Wave'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_534758_long_instab/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Mars_Gravity_Wave'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_535011_with_FTS/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Mars_Gravity_Wave'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_535489_removed_discontinuity_long/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Mars_Gravity_Wave/test_RAPHAEL'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Mars_Gravity_Wave'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_540064_FTS_no_disc_long/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Mars_Gravity_Wave'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_9078210_spread_source/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Mars_Gravity_Wave'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_9081352_spread_cut_source/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Mars_Gravity_Wave'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_9091089_new_coupling/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Mars_Gravity_Wave'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_9103256_same_as_previous_but_factor_1/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Mars_Gravity_Wave'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_552471_atmo_only/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Mars_Gravity_Wave'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_552455_factor0p1/');

% Seismic Hammer.
% fig_title = strcat('Seismic Hammer Simulation (New Model)');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Balloons/simulations'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_9048100_seismic_DG/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Balloons/simulations'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_9081476_seismic_potential/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Balloons/simulations'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_9091088_seismic_DG_new_coupling/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Balloons/simulations'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_9102702_seismic_potential_rem_forcing/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Balloons/simulations'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_9113508_seismic_DG_with_memvars_solid/');
% root_dir = '/home/l.martire/Documents/SPECFEM/Ongoing_Work/Balloons/simulations'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES_551980_seismic_potential_with_memvars_solid/');

% Tests.
fig_title = 'test';
% root_dir = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/full_DG_square'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES/');
root_dir = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_stretching'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES/');
% root_dir = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_FTS'; output_files_dir = strcat(root_dir, '/OUTPUT_FILES/');

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
if(not(strcmp(output_files_dir(end),'/'))); folder=[folder,'/']; end;
% Load sources' positions.
pos_sources = [inf, inf]; % Allocate a row for the first source's position.
% fid = fopen([root_dir, '/DATA/SOURCE']);
fid = fopen([output_files_dir, 'SOURCE']);
if(fid==-1)
  error(strcat("Cannot open SOURCE file in ", output_files_dir, '(',[output_files_dir, 'SOURCE'],').'))
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
  A = importdata(strcat(output_files_dir, 'STATIONS'));
catch
  try
    A = importdata(strcat(root_dir, '/DATA/STATIONS'));
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
% Stations to display. %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display stations informations, and ask user for input.
%disp([num2str(size(pos_stations,1)) ' stations found.']);
%[(1:size(pos_stations,1));[[0;0],(diff(pos_stations)~=0)']]
format shortG;
[(1:size(pos_stations, 1)); pos_stations';dist_to_sources']
format compact;
display_or_load=-1;
while(not(display_or_load==0 || display_or_load==1))
  display_or_load=input('  Display (0) or load only (1)? > ');
end
istattab = input(['  Stations (Matlab format, eg. [1, 4, 7] or 1:20)? > ']);
nstat = size(pos_stations(istattab, 1), 1);
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
% Plot.                       %
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
      unknown_name = 'vertical {displacement (m), velocity (m/s)}';
    elseif(strcmp(unknown,'BXX'))
      unknown_name = 'horizontal {displacement (m), velocity (m/s)}';
    else
      error("The variable 'unknown' has a non-standard value.");
    end
  elseif (type_display == 2)
    % Original SPECFEM2D's sismogram is velocity.
    extension = "semv"; % Because original SPECFEM2D's sismogram is velocity.
    if(strcmp(unknown,'BXZ'))
      unknown_name = '{vertical velocity (m/s), pressure perturbation (Pa)}';
    elseif(strcmp(unknown,'BXX'))
      unknown_name = '{horizontal velocity (m/s), pressure perturbation (Pa)}';
    else
      error("The variable 'unknown' has a non-standard value.");
    end
  end

  % Read the sismogram.
  file = strcat(output_files_dir, 'AA.', A.textdata(istat_glob, 1), '.', unknown, '.', extension);
  data = load(file{1});
  nt = max(size(data));
  % Sub-sample of records.
  nsub = ceil(nt/1000);
  nd = max(size(data(1:nsub:nt, 1)));
  % Recover time/amplitude data.
  % Ztime(istat,1:nt) = data(1:nsub:nt,1)';
  % Zamp(istat,1:nt)  = data(1:nsub:nt,2)';
  Ztime(istat, 1:nd) = data(1:nsub:nt, 1)';
  %Ztime(istat, 1:nd) = Ztime(istat, 1:nd) - Ztime(istat, 1); % Make time values start at zero.
  % Correct for sign of source function for N wave pattern (??)
  Zamp(istat, 1:nd) = data(1:nsub:nt, 2)';

  % Display.
  if(display_or_load==0) % If display.
    ax(istat) = subplot(nstat, 1, istat); plot(Ztime(1, :), Zamp(istat, :));

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

    legend(strcat('S', num2str(istat_glob), ', (x,z,d)=(', num2str(xstattab(istat_glob) / 1000), ',', num2str(ystattab(istat_glob) / 1000), ',', num2str(dist_to_sources(istat_glob) / 1000), ') km'), 'Location', 'west');

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
    normalise_ylims=-1;
    while(not(normalise_ylims==0 || normalise_ylims==1))
      normalise_ylims=input('  Normalise y-scale? (0 for no, 1 for yes) > ');
    end
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
else
  disp("  Data loaded. [id, station]:");
  [(1:length(istattab))',istattab']
  disp(strcat("  Example: Data of station ",num2str(istattab(1))," are in         Zamp(",num2str(1),", :)."));
  disp(strcat("           Corresponding time values are in Ztime(",num2str(1),", :)."));
end