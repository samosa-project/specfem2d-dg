% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

clear all;
% close all;
clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'defaultTextFontSize', 12); set(0, 'defaultAxesFontSize', 12);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

% User inputs.
disp(['User inputs.']);
% folder = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200/';
folder = input('  Folder containing the "external_bottom_forcing.dat" file? > ', 's');

if (not(strcmp(folder(end), '/')))
  folder = [folder, '/'];
end
gst = input('  Time grid size (s)? > ');
gsx = input('  Space grid size (m)? > ');
maxt = input('  Maximum t (s, Inf for all)? > ');
minx = input('  Minimum x (m, -Inf for all)? > ');
maxx = input('  Maximum x (m, Inf for all)? > ');

% Load the data.
EBFFILE = [folder, 'external_bottom_forcing.dat'];
disp(['Loading (can be very long for large files).']);
a = importdata(EBFFILE);
t = a(:, 1);
x = a(:, 2);
d = a(:, 3);
clear('a');
% selt=(t<maxt); t=t(selt); selx=(x>=minx & x<=maxx); x=x(selx);
% sel=(t<maxt & x>=minx & x<=maxx); t=t(sel); x=x(sel); d=d(sel);

% Interpolate the point cloud.
disp(['Interpolating (can be very long for large files).']);
F = scatteredInterpolant(t, x, d);

% Plot.
disp(['Plotting.']);
% tt = min(t(:)):gst:max(t(:));
% tx = min(x(:)):gsx:max(x(:));
tt = min(t(:)):gst:min(max(t(:)),maxt);
tx = max(min(x(:)),minx):gsx:min(max(x(:)),maxx);
[X, Z] = meshgrid(tt, tx);
D = F(X, Z);
figure();
surf(X, Z, D, 'edgecolor', 'interp', 'facecolor', 'interp'); view([0, 0, 1]); axis([min(tt), max(tt), min(tx), max(tx)]);
xlabel('$t$'); ylabel('$x$'); title("External Bottom Forcing ");
colorbar;
set(gca, 'TickLabelInterpreter','latex');