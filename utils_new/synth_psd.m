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
raw_t=Ztime(1,:);
raw_s=Zamp(1,:);
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
% Set signal to be used.
% TODO: ask for user input.
signal = raw_s; signal_name = "vertical velocity"; unit="(m/s)";
% signal = raw_s; signal_name = "displacement"; unit="m";
% signal = cumtrapz(raw_t, raw_s); signal_name = "displacement"; unit="m";

% Select time frame.
% TODO: ask for user input.
select_time_l=raw_t(1); select_time_u=raw_t(end);
% select_time_l=raw_t(1); select_time_u=48;
% select_time_l=0; select_time_u=2.9;

% Select.
% time=raw_t; time = time(select_time_l<=time); time=time(time<=select_time_u); signal = signal(time<=select_time_u); % Cut signal and time around selection window.
time=raw_t; signal(select_time_l>=time)=0; signal(time>=select_time_u)=0; % Zero signal value around selection window.
% signal=[signal, zeros(1,size(time,2))]; time=[time-time(1), time(end)-2*time(1)+time]; % Zero signal value around selection window and add zeros at the end.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ask for user input.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normalise_psd=-1;
while(not(length(normalise_psd)==1 && ismember(normalise_psd,[0,1])))
  normalise_psd=input('  Normalise PSD (0 for no, 1 for yes)? > ');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSD and plot.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Powerf,Freqf]=custom_psd(time, signal);
if(normalise_psd==1)
  Powerf=Powerf/max(Powerf);
end

figure();
plot(time, signal);
xlim([time(1), time(end)]);
xlabel("$t$ (s)"); ylabel(strcat(signal_name, " (",unit,")"));
title({strcat(signal_name, "")});

figure();
loglog(Freqf, Powerf);
xlim([Freqf(1), Freqf(end)]);
xlabel("$f$ (Hz)"); ylabel(strcat(signal_name, " PSD (",unit,"$^2$/Hz)"));
title({strcat(signal_name, " PSD")});

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear variables.             %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('raw_s','raw_t', 'select_time_l', 'select_time_u', 'signal_name', 'unit');