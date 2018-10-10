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
EBFFILE = input(['[',mfilename,'] Path to EBF file? > '], 's');
tmin = input(['[',mfilename,'] Minimum t (s, -Inf for all)? > ']);
tmax = input(['[',mfilename,'] Maximum t (s, Inf for all)? > ']);
xmin = input(['[',mfilename,'] Minimum x (m, -Inf for all)? > ']);
xmax = input(['[',mfilename,'] Maximum x (m, Inf for all)? > ']);

% Load the data.
disp(['[',mfilename,'] Loading (can be very long for large files).']);
if(tmin==-Inf && tmax==+Inf && xmin==-Inf && xmax==Inf)
  % Brutal but simple way.
  a = importdata(EBFFILE);
  t = a(:, 1);
  x = a(:, 2);
  d = a(:, 3);
  clear('a');
else
  % Smooth but complex way.
  f=fopen(EBFFILE,'r');
  if(f==-1)
    error(['[',mfilename,', ERROR] Cannot open file ''', EBFFILE, '''.']);
  end
  t=[]; x=[]; d=[];
  while(not(feof(f)))
    txv=str2num(fgetl(f)); % Vector [t, x, v].
    if(txv(1)>=tmin && txv(1)<=tmax)
      if(txv(2)>=xmin && txv(2)<=xmax)
        t=[t;txv(1)];
        x=[x;txv(2)];
        d=[d;txv(3)];
      end
    elseif(txv(1)>tmax)
      break; % We have all we need, get out.
    end
  end
  fclose('all');
end
dt=unique(diff(t))';
dx=unique(diff(x))';

% selt=(t<maxt); t=t(selt); selx=(x>=minx & x<=maxx); x=x(selx);
% sel=(t<maxt & x>=minx & x<=maxx); t=t(sel); x=x(sel); d=d(sel);

% More user inputs.
gst = input(['[',mfilename,'] Plotting time grid size (s)? File dt: [',num2str(dt),'] > ']);
gsx = input(['[',mfilename,'] Plotting space grid size (m)? File dx: [',num2str(dx),'] > ']);

% Interpolate the point cloud.
disp(['[',mfilename,'] Interpolating (can be very long for large files).']);
F = scatteredInterpolant(t, x, d);

% Plot.
disp(['[',mfilename,'] Plotting.']);
% tt = min(t(:)):gst:max(t(:));
% tx = min(x(:)):gsx:max(x(:));
tt = min(t(:)):gst:min(max(t(:)),tmax);
tx = max(min(x(:)),xmin):gsx:min(max(x(:)),xmax);
[X, Z] = meshgrid(tt, tx);
D = F(X, Z);
figure();
surf(X, Z, D, 'edgecolor', 'interp', 'facecolor', 'interp'); view([0, 0, 1]); axis([min(tt), max(tt), min(tx), max(tx)]);
xlabel('$t$'); ylabel('$x$'); title("External Bottom Forcing ");
colorbar;
set(gca, 'TickLabelInterpreter','latex');