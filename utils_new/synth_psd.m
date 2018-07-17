% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Computes PSDs and NPSDs of synthetics.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         /utils_new/synth_load.m should have been ran
%                before.

% clear all;
% close all;
% clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'defaultTextFontSize', 12); set(0, 'defaultAxesFontSize', 12);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: ask for user input.
raw_t=Ztime;
raw_s=Zamp;
% raw_t=data_leo_t(1,:);raw_s=data_leo_v(1,:);
% OKQ0
% raw_t=data_voon_t{4};raw_s=data_voon_v{4}'; fig_tit='Voon 15';
% raw_t=data_voon_t{5};raw_s=data_voon_v{5}'; fig_tit='Voon 30';
% raw_t=data_voon_t{6};raw_s=data_voon_v{6}'; fig_tit='Voon 45';
% raw_t=Ztime(1,:); raw_s=Zamp(1,:); fig_tit=fig_title;
% raw_t=Ztime(2,:); raw_s=Zamp(2,:); fig_tit=fig_title;
% raw_t=Ztime(3,:); raw_s=Zamp(3,:); fig_tit=fig_title;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters and              %
% pre-treatment.              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normalise_wpsd=-1;
while(not(length(normalise_wpsd)==1 && ismember(normalise_wpsd,[0,1])))
  normalise_wpsd=input('  Normalise Welch PSD (0 for no, 1 for yes)? > ');
end
if(normalise_wpsd)
  normalise_wpsd_txt=" (normalised)";
else
  normalise_wpsd_txt="";
end

% Eventually, compute integral or derivative.
% TODO: Ask for user input.
% signal = raw_s; signal_name = "displacement"; unit="m";
% signal = cumtrapz(raw_t, raw_s); signal_name = "displacement"; unit="m";

% TODO: Ask for user input.
% signal_name = "SIGNAL"; signal_unit="UNIT";
signal_name = "$\delta P$"; signal_unit="Pa";

if(nstat>1)
  avgwpsds=-1;
  while(not(length(avgwpsds)==1 && ismember(avgwpsds,[0,1,2])))
    avgwpsds=input('  Multiple data found. Choose first (0), average WPSDs (1), or compute WPSD of average signal (2)? > ');
  end
  switch(avgwpsds)
    case 0
      disp('  Will compute Welch PSD of first data.');
      WPSD_txt=strcat("Welch PSD of " ,signal_name, " [",signal_unit,"$^2$/Hz]",normalise_wpsd_txt);
      IDs_to_process=1;
    case 1
      disp('  Averaging Welch PSDs. Be wary of the stations you use.');
      WPSD_txt=strcat("Averaged Welch PSD of " ,signal_name, " [",signal_unit,"$^2$/Hz]",normalise_wpsd_txt);
      IDs_to_process=1:nstat;
    case 2
      disp('  Computing Welch PSD of average signal. Be wary of the stations you use.');
      WPSD_txt=strcat("Welch PSD of averaged " ,signal_name, " [",signal_unit,"$^2$/Hz]",normalise_wpsd_txt);
      IDs_to_process=1;
  end
else
  WPSD_txt=strcat("Welch PSD of " ,signal_name, " [",signal_unit,"$^2$/Hz]",normalise_wpsd_txt);
  IDs_to_process=1;
end

if(ismember(avgwpsds,[1,2]))
  if(avgwpsds==2)
    raw_s=mean(raw_s,1);
  end
  timeseries_txt=strcat("Averaged ", signal_name, " (",signal_unit,")");
else
  timeseries_txt=strcat(signal_name, " (",signal_unit,")");
end

WPSD_tab=[];
for i=IDs_to_process
  % Set signal to be used.
  signal = raw_s(i,:);

  % Select time frame.
  % TODO: ask for user input.
  select_time_l=raw_t(i,1); select_time_u=raw_t(i,end);
  % select_time_l=raw_t(1); select_time_u=48;
  % select_time_l=0; select_time_u=2.9;

  % Select.
  % time=raw_t; time = time(select_time_l<=time); time=time(time<=select_time_u); signal = signal(time<=select_time_u); % Cut signal and time around selection window.
  time=raw_t(i,:); signal(select_time_l>=time)=0; signal(time>=select_time_u)=0; % Zero signal value around selection window.
  % signal=[signal, zeros(1,size(time,2))]; time=[time-time(1), time(end)-2*time(1)+time]; % Zero signal value around selection window and add zeros at the end.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PSD and plot.               %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [WPSD,WPSD_f]=custom_psd(time, signal);

  WPSD_tab(i,:)=WPSD;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Stack".                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeseries_to_plot = mean(raw_s(IDs_to_process,:),1); % Will be the first time series loaded if stack was deactivated, or the only time series if only one station was loaded.
WPSD_to_plot = mean(WPSD_tab(IDs_to_process,:),1); % Will be the PSD of the first time series if stack was deactivated, or the PSD of the only time series if only one station was loaded.

if(normalise_wpsd==1)
  WPSD_to_plot=WPSD_to_plot/max(WPSD_to_plot);
end

figure();
plot(time, timeseries_to_plot);
xlim([time(1), time(end)]);
xlabel("$t$ (s)"); ylabel(timeseries_txt);
title(timeseries_txt);

figure();
loglog(WPSD_f, WPSD_to_plot);
xlim([WPSD_f(1), WPSD_f(end)]);
xlabel("$f$ (Hz)"); ylabel(WPSD_txt);
title(WPSD_txt);

disp(sprintf("Amplitude of signal: %1.6e",max(timeseries_to_plot)-min(timeseries_to_plot)));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear variables.             %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('select_time_l', 'select_time_u', 'signal_name', 'unit');