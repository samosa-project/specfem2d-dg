% Author:        Léo Martire.
% Description:   Filter synthetics.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         /utils_new/synth_load.m has to have been ran
%                before.

% clear all;
% close all;
% clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'defaultTextFontSize', 12); set(0, 'defaultAxesFontSize', 12);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

if(not(exist('synth_load_was_ran') && synth_load_was_ran==1))
  error("[ERROR] synth_load was not ran before.");
end

% Ask user for filter, low, and/or high cut-off frequencies.
% TODO.
f_low=1;
f_high=20;

% Low-pass.
% TODO.

% High-pass.
% TODO.

% Band-pass.
for i=1:nstat
  [~, signal_HP]=custom_filter(Ztime(i,:), Zamp(i,:), f_low);
  [signal_HP_LP, ~]=custom_filter(Ztime(i,:), signal_HP, f_high);
  Zamp(i,:)=signal_HP_LP;
end