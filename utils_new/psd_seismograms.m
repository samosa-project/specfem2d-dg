% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         display_seismograms.m should have been ran before.

% clear all
% close all
% clc
format compact;
set(0, 'DefaultLineLineWidth', 2); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 8); % Default at 6.
set(0, 'defaultTextFontSize', 12);
set(0, 'defaultAxesFontSize', 12); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

% Load raw seismogram.
% TODO: ask for user input.
% raw_t=Ztime(1,:);
% raw_s=Zamp(1,:);
% raw_t=data_leo_t(1,:);raw_s=data_leo_v(1,:);

% OKQ0
% raw_t=data_voon_t{4};raw_s=data_voon_v{4}'; fig_tit='Voon 15';
% raw_t=data_voon_t{5};raw_s=data_voon_v{5}'; fig_tit='Voon 30';
% raw_t=data_voon_t{6};raw_s=data_voon_v{6}'; fig_tit='Voon 45';
% raw_t=Ztime(1,:); raw_s=Zamp(1,:); fig_tit=fig_title;
% raw_t=Ztime(2,:); raw_s=Zamp(2,:); fig_tit=fig_title;
raw_t=Ztime(3,:); raw_s=Zamp(3,:); fig_tit=fig_title;

% Parameters.
% TODO: ask for user input.
% select_time_l=raw_t(1); select_time_u=raw_t(end);
select_time_l=raw_t(1); select_time_u=48;
% select_time_l=0; select_time_u=2.9;

% Set signal to be used.
% TODO: ask for user input.
signal = raw_s; signal_name = "vertical velocity"; unit="(m/s)";
% signal = raw_s; signal_name = "displacement"; unit="m";
% signal = cumtrapz(raw_t, raw_s); signal_name = "displacement"; unit="m";

% Select.
% time=raw_t; time = time(select_time_l<=time); time=time(time<=select_time_u); signal = signal(time<=select_time_u); % Cut signal and time around selection window.
time=raw_t; signal(select_time_l>=time)=0; signal(time>=select_time_u)=0; % Zero signal value around selection window.
% signal=[signal, zeros(1,size(time,2))]; time=[time-time(1), time(end)-2*time(1)+time]; % Zero signal value around selection window and add zeros at the end.

% PSD.
signal=detrend(signal); % Remove eventual linear trend.
dt=time(2)-time(1);
Fls=1/dt;
% nfft= 2^(floor(log2(2000)))/2;
nfft= 2^(8); nfft=min(length(raw_s), nfft);
window = hann(nfft);
noverlap= nfft/2;

[Powerf,Freqf]=pwelch(signal,window,noverlap,nfft,Fls);

% Plots.
figure();
plot(time, signal);
xlim([time(1), time(end)]);
xlabel("$t$ (s)"); ylabel(strcat(signal_name, " (",unit,")"));
title({fig_tit,strcat(signal_name, "")});

figure();
loglog(Freqf, Powerf);
xlim([Freqf(1), Freqf(end)]);
xlabel("$f$ (Hz)"); ylabel(strcat(signal_name, " PSD (",unit,"$^2$/Hz)"));
title({fig_tit,strcat(signal_name, " PSD")});