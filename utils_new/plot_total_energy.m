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

folder = input('  Folder containing the energy.dat file? > ', 's');
if(not(strcmp(folder(end),'/'))); folder=[folder,'/']; end;
logscale=-1;
while(not(logscale==0 || logscale==1))
  logscale=input('  Normal y-scale (0) or log y-scale (1)? > ');
end

% Load the data.
a=importdata(strcat(folder, 'energy.dat'));
t=a.data(:,1); ke=a.data(:,2); pe=a.data(:,3); te=a.data(:,4);

figure();
if(logscale)
  semilogy(t, te);
else
  plot(t, te);
end
xlabel("t (s)"); ylabel("total energy (J)");
title("Total simulation energy");