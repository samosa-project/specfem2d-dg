% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Convolves synthetics by a function.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOT FUNCTIONNAL FOR NOW.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'defaultTextFontSize', 12); set(0, 'defaultAxesFontSize', 12);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

SPCFMloc='/home/l.martire/Documents/SPECFEM/';

% Quantity to display:
%   1 = displacement for non-DG and velocity for DG,
%   2 = velocity for non-DG and pressure perturbation (Pa) for DG.
type_display = 2; % Should be the same as the seismotype variable in parfile.

% Unknown:
% For type_display==2 and stations in DG zones, pressure perturbation (Pa) is saved both in BXX and BXZ files.
% unknown = 'BXX';
unknown = 'BXZ';

% Root directory.
rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45');

% Directory of base synthetics (Dirac? Heaviside?).
synOFd = strcat(rootd, '/OUTPUT_FILES_stf4/'); type_display = 1; % Dirac.
% synOFd = strcat(rootd, '/OUTPUT_FILES_sharpstf3/'); type_display = 1; % Sharp Gaussian.
% synOFd = strcat(rootd, '/OUTPUT_FILES_stf5/'); type_display = 1; % Heaviside

% Directory where is save the wanted STF.
% TODO: Ask user for custom STF.
stfOFd = strcat(rootd, '/OUTPUT_FILES_stf3/'); % Gaussian.
% stfOFd = strcat(rootd, '/OUTPUT_FILES_stf2/'); % Gaussian derivative.

station=input('  Station number?\n  > ');

if (type_display == 1)
  extension = "semd";
elseif (type_display == 2)
  extension = "semv";
end

stfd=importdata(strcat(stfOFd,"plot_source_time_function.txt"));
stftime=stfd(stfd(:,1)<=-stfd(1,1),1);
stfvals=stfd(stfd(:,1)<=-stfd(1,1),2);

figure();
plot(stftime,stfvals);
xlim([min(stftime),max(stftime)]);
title("STF");

STATIONSFILE = importdata(strcat(synOFd, 'STATIONS'));
data = load(strcat(synOFd, 'AA.', STATIONSFILE.textdata(station, 1), '.', unknown, '.', extension));
syntime = data(:, 1)';
synvals = data(:, 2)';

figure();
plot(syntime,synvals);
xlim([min(syntime),max(syntime)]);
title("Synthetic");

convv=conv(synvals,stfvals,"same");

figure();
plot(syntime,convv);
xlim([min(syntime),max(syntime)]);
title("Synthetic convolved with STF");