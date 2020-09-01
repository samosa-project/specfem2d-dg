% Author:        Léo Martire.
% Description:   Computes ASDs and NASDs of synthetics.
% Usage:         N. A.
% Notes:         synth_load.m should have been ran before.

% clear all;
% close all;
% clc;

[~] = setup_overall();

% prepare window: Cb = polyfit([4e3, 124e3],[-0.4, 451.3],1); Ca = polyfit([0, 120e3],[20, 488],1);
% window: dsorted(i)*[Cb(1), Ca(1)] + [Cb(2), Ca(2)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dsorted, isort] = sort(distance(istattab));
raw_t = Ztime(isort,:); raw_s = Zamp(isort,:); cd(OFd);

% TODO: ask for user input.
% raw_t = Ztime; raw_s = Zamp; cd(OFd);
% raw_t = Ztime(1:2:end,:); raw_s = Zamp(1:2:end,:); cd(OFd);
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
coord_units = 'km';
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


nstat = min(size(raw_t));

NWINDOWZZZZ = -1;
while (not(length(NWINDOWZZZZ) == 1 && NWINDOWZZZZ>0))
  NWINDOWZZZZ = input(['[', mfilename, '] Number of windows for ASD (integer, >0, higher values <=> smoother curve & higher lowest frequency)? > ']);
end

normalise_wpsd = - 1;
while (not(length(normalise_wpsd) == 1 && ismember(normalise_wpsd, [0, 1])))
  normalise_wpsd = input(['[', mfilename, '] Normalise ASD (0 for no, 1 for yes)? > ']);
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
    avgwpsds = input(['[', mfilename, '] Multiple data found. Choose first (0), average WASDs (1), compute WASD of average signal (2), plot every WASD on top of each other (3), or plot every ASD as surf (4)? > ']);
  end
  switch(avgwpsds)
    case 0
      disp(['[', mfilename, '] Will compute ASD of first data.']);
      WASD_txt = strcat("ASD of ", signal_name, " [", signal_unit, "/$\sqrt{\mathrm{Hz}}$]", normalise_wpsd_txt);
      IDs_to_process = 1;
    case 1
      disp(['[', mfilename, '] Averaging ASDs. Be wary of the stations you use.']);
      WASD_txt = strcat("Average of ASDs of ", signal_name, " [", signal_unit, "/$\sqrt{\mathrm{Hz}}$]", normalise_wpsd_txt);
      IDs_to_process = 1:nstat;
    case 2
      disp(['[', mfilename, '] Computing ASD of average signal. Be wary of the stations you use.']);
      WASD_txt = strcat("ASD of averaged ", signal_name, " [", signal_unit, "/$\sqrt{\mathrm{Hz}}$]", normalise_wpsd_txt);
      IDs_to_process = 1;
    case 3
      disp(['[', mfilename, '] Plotting every ASD on top of each other.']);
      WASD_txt = strcat("ASDs of ", signal_name, " [", signal_unit, "/$\sqrt{\mathrm{Hz}}$]", normalise_wpsd_txt);
      IDs_to_process = 1:nstat;
    case 4
      disp(['[', mfilename, '] Plotting every ASD as surf.']);
      WASD_txt = strcat("ASDs of ", signal_name, " [", signal_unit, "/$\sqrt{\mathrm{Hz}}$]", normalise_wpsd_txt);
      IDs_to_process = 1:nstat;
  end
else
  avgwpsds = 0;
  WASD_txt = strcat("ASD of ", signal_name, " [", signal_unit, "/$\sqrt{\mathrm{Hz}}$]", normalise_wpsd_txt);
  IDs_to_process = 1;
end

if (ismember(avgwpsds, [1, 2]))
  if (avgwpsds == 2)
    raw_s = mean(raw_s, 1);
  end
  timeseries_txt = strcat("Averaged ", signal_name, " (", signal_unit, "$^2$)");
else
  timeseries_txt = strcat(signal_name, " (", signal_unit, "$^2$)");
end

WASD_tab = [];
stopAskingWindowing = 0;
for i = IDs_to_process
  % Set signal to be used.
  signal = raw_s(i, :);
  
  if(not(stopAskingWindowing))
    % Select time frame.
    doWindowSignal = -Inf;
    while (not((numel(doWindowSignal)==1 && ismember(doWindowSignal,[0,-1])) || (numel(doWindowSignal)==2)))
      disp(['[', mfilename, '] Window signal for S',num2str(istattab(i)),' @[',num2str(xstattab(istattab(i))),',',num2str(ystattab(istattab(i))),'] ([t1 t2] vector for yes, 0 for no, -1 for no for all)?']);
%       disp(['[', mfilename, ', INFO] You may use ''xstattab(istattab(i))'' to get the value of the the station''s x position.']);
      disp(['[', mfilename, ', INFO] You may use ''dsorted(i)'' to get the value of the the station''s x position.']);
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
  % ASD and plot.               %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [~, ~, ~, ~, ASD, ASD_f] = PSD_Spectrogram(signal, mean(1./diff(time)), NWINDOWZZZZ); % UNIT/HZ^.5
  [ASD_f, ASD] = asd(signal, mean(1./diff(time)), 'nwindow', NWINDOWZZZZ); % UNIT/HZ^.5
%   [ASD, ASD_f] = custom_psd(time, signal); % UNIT^2/HZ

  WASD_tab(i, :) = ASD;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Stack".                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (ismember(avgwpsds, [0, 1, 2]))
  timeseries_to_plot = mean(raw_s(IDs_to_process, :), 1); % Will be the first time series loaded if stack was deactivated, or the only time series if only one station was loaded.
  WASD_to_plot = mean(WASD_tab(IDs_to_process, :), 1); % Will be the ASD of the first time series if stack was deactivated, or the ASD of the only time series if only one station was loaded.

  if (normalise_wpsd == 1)
    WASD_to_plot = WASD_to_plot / max(WASD_to_plot);
  end

  figure();
  plot(time, timeseries_to_plot);
  xlim([time(1), time(end)]);
  xlabel("$t$ (s)"); ylabel(timeseries_txt);
  title(timeseries_txt);
  set(gca, 'TickLabelInterpreter', 'latex');

  figure();
  loglog(ASD_f, WASD_to_plot);
  xlim([ASD_f(1), ASD_f(end)]);
  xlabel("$f$ [Hz]"); ylabel(WASD_txt);
  title(WASD_txt);
  set(gca, 'TickLabelInterpreter', 'latex');
  grid;

  disp(['[', mfilename, '] Amplitude of signal: ',sprintf('%1.6e',max(timeseries_to_plot)-min(timeseries_to_plot))]);
elseif (avgwpsds == 3)
  figure();
  colours = jet(numel(IDs_to_process));
  for i = IDs_to_process
    if (strcmp(coord_units, 'km'))
      dFactor = 1e-3;
    else
      dFactor = 1;
    end
%       ASDName = strcat('S', num2str(istattab(i)), ', $(x,z) = (', num2str(xstattab(istattab(i)) / 1000), ',', num2str(ystattab(istattab(i)) / 1000), "$) ", coord_units);
%     elseif (strcmp(coord_units, 'm'))
%     ASDName = ['S', num2str(istattab(i)), ', $(x,z) = (', num2str(xstattab(istattab(i))), ',', num2str(ystattab(istattab(i))), ')$ [m]'];
    ASDName = ['S', num2str(istattab(i)), '@$d=', num2str(dsorted(i)*dFactor), '$ [',coord_units,']'];
%     else
%       error(['[', mfilename, ', ERROR] coord_units = ', coord_units, 'not implemented.']);
%     end
    loglog(ASD_f, WASD_tab(i, :), 'displayname', ASDName, 'color', colours(i, :));
    hold on;
  end
  legend('location', 'best');
  xlim([ASD_f(1), ASD_f(end)]);
  xlabel("$f$ [Hz]"); ylabel(WASD_txt);
  title(WASD_txt);
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
  SURFy = ASD_f;
  [SURFX, SURFY] = meshgrid(SURFx, SURFy);
  surf(SURFX, SURFY, log10(WASD_tab'));
  shading interp;
  xlim([min(SURFx), max(SURFx)]);
  ylim([min(SURFy), max(SURFy)]);
  set(gca, 'yscale', 'log');
  xlabel(['$', dist_symbol, '$ [m]']);
  ylabel(['$f$ [Hz]']);
  %   title(strcat("$\log($",WASD_txt,"$)$"));
  title(WASD_txt);
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
    loglog(ASD_f, WASD_tab(i, :)./WASD_tab(refID, :), 'displayname', ['S', num2str(istattab(i)),'/S', num2str(istattab(refID)), '@$(x,z) = (', num2str(xstattab(istattab(i))), ',', num2str(ystattab(istattab(i))), ')$ [m]'], 'color', colours(i, :));
    hold on;
  end
  legend('location', 'best');
  xlim([ASD_f(1), ASD_f(end)]);
  xlabel("$f$ [Hz]"); ylabel('ASD Ratio [1]');
  title(['ASD Ratio [1] w.r.t. S', num2str(istattab(refID))]);
  % set(gca, 'TickLabelInterpreter', 'latex'); grid on; box on;
  prettyAxes(gcf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear variables.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('select_time_l', 'select_time_u', 'signal_name', 'unit');