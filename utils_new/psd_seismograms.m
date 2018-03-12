% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

% clear all
% close all
% clc
format compact;
set(0, 'DefaultLineLineWidth', 1.5); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 6); % Default at 6.
set(0, 'defaultTextFontSize', 16);
set(0, 'defaultAxesFontSize', 14); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

% Load raw seismograms.
raw_t=Ztime(1,:);
raw_s=Zamp(1,:);

% Parameters.
select_time_l=0;
select_time_u=45;

% Set signal to be used.
signal = cumtrapz(raw_t, raw_s);
signal_name = "displacement"; unit="m";

% Select.
% time=raw_t; time = time(select_time_l<=time); time=time(time<=select_time_u); signal = signal(time<=select_time_u); % Cut signal and time around selection window.
time=raw_t; signal(select_time_l>=time)=0; signal(time>=select_time_u)=0; % Zero signal value around selection window.
% signal=[signal, zeros(1,size(time,2))]; time=[time-time(1), time(end)-2*time(1)+time]; % Zero signal value around selection window and add zeros at the end.

% PSD.
signal=detrend(signal); % Remove eventual linear trend.
dt=time(2)-time(1);
Fls=1/dt;
% nfft= 2^(floor(log2(2000)))/2;
nfft= 2^(12);
window = hann(nfft);
noverlap= nfft/2;

[Powerf,Freqf]=pwelch(signal,window,noverlap,nfft,Fls);

% Plots.
figure();
plot(time, signal);
xlim([time(1), time(end)]);
xlabel("$t$ (s)"); ylabel(strcat(signal_name, " (",unit,")"));
title({fig_title,strcat(signal_name, "")});

figure();
loglog(Freqf, Powerf);
xlim([Freqf(1), Freqf(end)]);
xlabel("$f$ (Hz)"); ylabel(strcat(signal_name, " PSD (",unit,"$^2$/Hz)"));
title({fig_title,strcat(signal_name, " PSD")});