% Author:        Léo Martire.
% Description:   Computes PSDs and NPSDs of synthetics.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         /utils_new/synth_load.m should have been ran before.

% clear all;
% close all;
% clc;
format compact;
set(0, 'DefaultLineLineWidth', 3); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'defaultTextFontSize', 12); set(0, 'defaultAxesFontSize', 18);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: ask for user input.
raw_t = Ztime; raw_s = Zamp; cd(OFd);
% raw_t = stf_from_run(:, 1)'; raw_s = stf_from_run(:, 2)';
% raw_t = data_leo_t(1,:);raw_s = data_leo_v(1,:);
% OKQ0
% raw_t = data_voon_t{1};raw_s = data_voon_v{1}'; fig_tit = 'Voon 53m'; nstat = 1;
% raw_t = data_voon_t{12};raw_s = data_voon_v{12}'; fig_tit = 'Voon 102m'; nstat = 1;
% raw_t = data_voon_t{14};raw_s = data_voon_v{14}'; fig_tit = 'Voon 151m'; nstat = 1;
% raw_t = data_voon_t{17};raw_s = data_voon_v{17}'; fig_tit = 'Voon 297m'; nstat = 1;
% raw_t = data_voon_t{4};raw_s = data_voon_v{4}'; fig_tit = 'Voon 15';
% raw_t = data_voon_t{5};raw_s = data_voon_v{5}'; fig_tit = 'Voon 30';
% raw_t = data_voon_t{6};raw_s = data_voon_v{6}'; fig_tit = 'Voon 45';
% raw_t = Ztime(1,:); raw_s = Zamp(1,:); fig_tit = fig_title;
% raw_t = Ztime(2,:); raw_s = Zamp(2,:); fig_tit = fig_title;
% raw_t = Ztime(3,:); raw_s = Zamp(3,:); fig_tit = fig_title;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters and              %
% pre-treatment.              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(not(all(size(raw_t) == size(raw_s))))
  error(['[(',mfilename,', ERROR] time data and amplitude should have the same size, but right now do not.']);
end

% doWindowSignal = -1;
% while (not((numel(doWindowSignal)==1 && doWindowSignal==0) || (numel(doWindowSignal)==2)))
%   doWindowSignal = input(['[', mfilename, '] Window signal ([t1 t2] vector for yes, 0 for no)? > ']);
% end
% if(any(doWindowSignal))
%   sel = (raw_t(1,:)>=min(doWindowSignal) & raw_t(1,:)<=max(doWindowSignal));
%   raw_t = raw_t(:,sel);
%   raw_s = raw_s(:,sel);
% end

nstat = min(size(Ztime));

NWINDOWZZZZ = -1;
while (not(length(NWINDOWZZZZ) == 1 && NWINDOWZZZZ>0))
  NWINDOWZZZZ = input(['[', mfilename, '] Number of windows for PSD (integer, >0, higher values <=> smoother curve & higher lowest frequency)? > ']);
end

normalise_wpsd = - 1;
while (not(length(normalise_wpsd) == 1 && ismember(normalise_wpsd, [0, 1])))
  normalise_wpsd = input(['[', mfilename, '] Normalise PSD (0 for no, 1 for yes)? > ']);
end
if (normalise_wpsd)
  normalise_wpsd_txt = " (normalised)";
else
  normalise_wpsd_txt = "";
end

% Eventually, compute integral or derivative.
% TODO: Ask for user input.
% signal = raw_s; signal_name = "displacement"; unit = "m";
% signal = cumtrapz(raw_t, raw_s); signal_name = "displacement"; unit = "m";

% TODO: Ask for user input.
signal_name = "SIGNAL"; signal_unit = "UNIT";
% signal_name = "$\delta P$"; signal_unit = "Pa";

avgwpsds = - 1;
if (nstat > 1)
  while (not(length(avgwpsds) == 1 && ismember(avgwpsds, [0, 1, 2, 3, 4])))
    avgwpsds = input(['[', mfilename, '] Multiple data found. Choose first (0), average WPSDs (1), compute WPSD of average signal (2), plot every WPSD on top of each other (3), or plot every PSD as surf (4)? > ']);
  end
  switch(avgwpsds)
    case 0
      disp(['[', mfilename, '] Will compute PSD of first data.']);
      WPSD_txt = strcat("PSD of ", signal_name, " [", signal_unit, "/Hz$^{.5}$]", normalise_wpsd_txt);
      IDs_to_process = 1;
    case 1
      disp(['[', mfilename, '] Averaging PSDs. Be wary of the stations you use.']);
      WPSD_txt = strcat("Average of PSDs of ", signal_name, " [", signal_unit, "/Hz$^{.5}$]", normalise_wpsd_txt);
      IDs_to_process = 1:nstat;
    case 2
      disp(['[', mfilename, '] Computing PSD of average signal. Be wary of the stations you use.']);
      WPSD_txt = strcat("PSD of averaged ", signal_name, " [", signal_unit, "/Hz$^{.5}$]", normalise_wpsd_txt);
      IDs_to_process = 1;
    case 3
      disp(['[', mfilename, '] Plotting every PSD on top of each other.']);
      WPSD_txt = strcat("PSDs of ", signal_name, " [", signal_unit, "/Hz$^{.5}$]", normalise_wpsd_txt);
      IDs_to_process = 1:nstat;
    case 4
      disp(['[', mfilename, '] Plotting every PSD as surf.']);
      WPSD_txt = strcat("PSDs of ", signal_name, " [", signal_unit, "/Hz$^{.5}$]", normalise_wpsd_txt);
      IDs_to_process = 1:nstat;
  end
else
  avgwpsds = 0;
  WPSD_txt = strcat("PSD of ", signal_name, " [", signal_unit, "/Hz$^{.5}$]", normalise_wpsd_txt);
  IDs_to_process = 1;
end

if (ismember(avgwpsds, [1, 2]))
  if (avgwpsds == 2)
    raw_s = mean(raw_s, 1);
  end
  timeseries_txt = strcat("Averaged ", signal_name, " (", signal_unit, ")");
else
  timeseries_txt = strcat(signal_name, " (", signal_unit, ")");
end

WPSD_tab = [];
stopAskingWindowing = 0;
for i = IDs_to_process
  % Set signal to be used.
  signal = raw_s(i, :);
  
  if(not(stopAskingWindowing))
    % Select time frame.
    doWindowSignal = -Inf;
    while (not((numel(doWindowSignal)==1 && ismember(doWindowSignal,[0,-1])) || (numel(doWindowSignal)==2)))
      disp(['[', mfilename, '] Window signal for S',num2str(istattab(i)),' @[',num2str(xstattab(istattab(i))),',',num2str(ystattab(istattab(i))),'] ([t1 t2] vector for yes, 0 for no, -1 for no for all)?']);
      disp(['[', mfilename, ', INFO] You may use ''xstattab(istattab(i))'' to get the value of the the station''s x position.']);
      doWindowSignal = input(['[', mfilename, '] > ']);
    end
    if(numel(doWindowSignal)==1 && doWindowSignal==-1)
      stopAskingWindowing = 1;
      doWindowSignal = 0;
    end
  end
  if(any(doWindowSignal))
    sel = (raw_t(1,:)>=min(doWindowSignal) & raw_t(1,:)<=max(doWindowSignal));
%     raw_t = raw_t(:,sel);
%     raw_s = raw_s(:,sel);
    disp(['[',mfilename,', INFO] Windowing: [',num2str(doWindowSignal),'] [s].']);
    select_time_l = min(doWindowSignal); select_time_u = max(doWindowSignal); % no window
  else
    select_time_l = raw_t(i, 1); select_time_u = raw_t(i, end); % no window
  end
  % select_time_l = raw_t(1); select_time_u = 48;
  % select_time_l = 0; select_time_u = 2.9;

  % Select.
  % time = raw_t; time = time(select_time_l<= time); time = time(time<= select_time_u); signal = signal(time<= select_time_u); % Cut signal and time around selection window.
%   time = raw_t(i, :); signal(select_time_l >=  time) = 0; signal(time >=  select_time_u) = 0; % Zero signal value around selection window.
  time = raw_t(i, :); signal(time <= select_time_l) = 0; signal(time >=  select_time_u) = 0; % Zero signal value around selection window.
  % signal = [signal, zeros(1,size(time,2))]; time = [time-time(1), time(end)-2*time(1)+time]; % Zero signal value around selection window and add zeros at the end.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PSD and plot.               %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [~, ~, ~, ~, PSD, PSD_f] = PSD_Spectrogram(signal, mean(1./diff(time)), NWINDOWZZZZ); % UNIT/HZ^.5
%   [PSD, PSD_f] = custom_psd(time, signal); % UNIT^2/HZ

  WPSD_tab(i, :) = PSD;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Stack".                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (ismember(avgwpsds, [0, 1, 2]))
  timeseries_to_plot = mean(raw_s(IDs_to_process, :), 1); % Will be the first time series loaded if stack was deactivated, or the only time series if only one station was loaded.
  WPSD_to_plot = mean(WPSD_tab(IDs_to_process, :), 1); % Will be the PSD of the first time series if stack was deactivated, or the PSD of the only time series if only one station was loaded.

  if (normalise_wpsd == 1)
    WPSD_to_plot = WPSD_to_plot / max(WPSD_to_plot);
  end

  figure();
  plot(time, timeseries_to_plot);
  xlim([time(1), time(end)]);
  xlabel("$t$ (s)"); ylabel(timeseries_txt);
  title(timeseries_txt);
  set(gca, 'TickLabelInterpreter', 'latex');

  figure();
  loglog(PSD_f, WPSD_to_plot);
  xlim([PSD_f(1), PSD_f(end)]);
  xlabel("$f$ [Hz]"); ylabel(WPSD_txt);
  title(WPSD_txt);
  set(gca, 'TickLabelInterpreter', 'latex');
  grid;

  disp(['[', mfilename, '] Amplitude of signal: ',sprintf('%1.6e',max(timeseries_to_plot)-min(timeseries_to_plot))]);
elseif (avgwpsds == 3)
  figure();
  colours = jet(numel(IDs_to_process));
  for i = IDs_to_process
%     if (strcmp(coord_units, 'km'))
%       PSDName = strcat('S', num2str(istattab(i)), ', $(x,z) = (', num2str(xstattab(istattab(i)) / 1000), ',', num2str(ystattab(istattab(i)) / 1000), "$) ", coord_units);
%     elseif (strcmp(coord_units, 'm'))
    PSDName = ['S', num2str(istattab(i)), ', $(x,z) = (', num2str(xstattab(istattab(i))), ',', num2str(ystattab(istattab(i))), ')$ [m]'];
%     else
%       error(['[', mfilename, ', ERROR] coord_units = ', coord_units, 'not implemented.']);
%     end
    loglog(PSD_f, WPSD_tab(i, :), 'displayname', PSDName, 'color', colours(i, :));
    hold on;
  end
  legend('location', 'best');
  xlim([PSD_f(1), PSD_f(end)]);
  xlabel("$f$ [Hz]"); ylabel(WPSD_txt);
  title(WPSD_txt);
%   set(gca, 'TickLabelInterpreter', 'latex'); grid on; box on;
  prettyAxes(gcf);
elseif (avgwpsds == 4)
  distancechoice = - 1;
  while (~ ismember(distancechoice, [1, 2, 3, 4]))
    distancechoice = input(['[', mfilename, '] Distance choice? (1 for x, 2 for |x|, 3 for z, 4 for d) > ']);
  end
  figure();
%   renorm_for_unit = - 1;
%   if (strcmp(coord_units, 'km'))
%     renorm_for_unit = 1000;
%   elseif (strcmp(coord_units, 'm'))
  renorm_for_unit = 1;
%   else
%     error(['[', mfilename, ', ERROR] coord_units = ', coord_units, 'not implemented.']);
%   end
  switch distancechoice
    case 1
      SURFx = xstattab(istattab(IDs_to_process)) / renorm_for_unit; dist_symbol = 'x';
    case 2
      SURFx = abs(xstattab(istattab(IDs_to_process))) / renorm_for_unit; dist_symbol = '|x|';
    case 3
      SURFx = ystattab(istattab(IDs_to_process)) / renorm_for_unit; dist_symbol = 'z';
    case 4
      SURFx = dist_to_sources(istattab(IDs_to_process)) / renorm_for_unit; dist_symbol = 'd';
  end
  SURFy = PSD_f;
  [SURFX, SURFY] = meshgrid(SURFx, SURFy);
  surf(SURFX, SURFY, log10(WPSD_tab'));
  shading interp;
  xlim([min(SURFx), max(SURFx)]);
  ylim([min(SURFy), max(SURFy)]);
  set(gca, 'yscale', 'log');
  xlabel(['$', dist_symbol, '$ [m]']);
  ylabel(['$f$ [Hz]']);
  %   title(strcat("$\log($",WPSD_txt,"$)$"));
  title(WPSD_txt);
  set(gca, 'TickLabelInterpreter', 'latex');
  view([0, 0, 1]);
  grid;
  cb = colorbar;
  set(cb, 'TickLabelInterpreter', 'latex');
  set(cb, 'ticklabels', split(sprintf('$10^{%g}$ ', get(cb, 'ticks'))));
else
  error(['[', mfilename, ', ERROR] [', mfilename, ', ERROR] bad value for avgwpsds.']);
end

% Plot ratios.
if(nstat > 1)
  refID = 1;
  % if(0)
  IDs_to_process_RATIO = IDs_to_process;
  IDs_to_process_RATIO(IDs_to_process_RATIO == refID) = []; % remove ref.
  figure();
  for i = IDs_to_process_RATIO
    loglog(PSD_f, WPSD_tab(i, :)./WPSD_tab(refID, :), 'displayname', ['S', num2str(istattab(i)),'/S', num2str(istattab(refID)), '@$(x,z) = (', num2str(xstattab(istattab(i))), ',', num2str(ystattab(istattab(i))), ')$ [m]'], 'color', colours(i, :));
    hold on;
  end
  legend('location', 'best');
  xlim([PSD_f(1), PSD_f(end)]);
  xlabel("$f$ [Hz]"); ylabel('PSD Ratio [1]');
  title(['PSD Ratio [1] w.r.t. S', num2str(istattab(refID))]);
  % set(gca, 'TickLabelInterpreter', 'latex'); grid on; box on;
  prettyAxes(gcf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear variables.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('select_time_l', 'select_time_u', 'signal_name', 'unit');


function prettyAxes(f)
  children=f.Children;
  for i=1:numel(children)
    child=children(i);
    if(strcmp(child.Type,'axes'))
      axes(child);
      set(gca, 'Color','k');
      set(gca, 'GridColor','white');
      set(gca, 'TickLabelInterpreter', 'latex');
      set(gca, 'TickDir','both');
      set(gca, 'TickLabelInterpreter', 'latex');
      grid on;
      box on;
    elseif(strcmp(children(i).Type,'legend'))
      set(child,'fontsize', 25);
      set(child, 'Color',[1,1,1]*0.25);
      set(child, 'textColor','w');
    end
  end
end