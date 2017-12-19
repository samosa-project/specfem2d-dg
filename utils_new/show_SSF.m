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

folder = input('Folder containing the SSF files? > ', 's');
if(not(strcmp(folder(end),'/'))); folder=[folder,'/']; end;

% Load the data.
listSSF=dir(strcat(folder, 'SSF*'));
a=[];
for i=1:length(listSSF)
  a=[a;importdata(strcat(folder, listSSF(i).name))];
  if(mod(i,floor(length(listSSF)/10))==0 || i==length(listSSF)); disp(strcat('Loading data (', num2str(100*i/length(listSSF)), '%).')); end;
end
x=a(:,1); z=a(:,2); d=a(:,3);

% Interpolate the point cloud.
gs = input('Plotting grid size? > ');
disp(['Interpolating.'])
F = scatteredInterpolant(x,z,d);
tx = min(x(:)):gs:max(x(:));
tz = min(z(:)):gs:max(z(:));
[X,Z] = meshgrid(tx,tz);
D = F(X,Z);
ar=input('Axes aspect ratio (format [ar_x, ar_y, ar_z])? > ');
figure(1);
surf(X,Z,D,'FaceColor', 'interp', 'EdgeColor', 'white', 'LineStyle', ':'); set(gca,'DataAspectRatio',ar);
xlabel('x'); ylabel('z'); title('Source Spatial Function');
th = input('Threshold to recenter plot (show only where SSF>threshold)? > ');
xlim([min(x(d>th)), max(x(d>th))]); ylim([min(z(d>th)), max(z(d>th))]); zlim([th, max(d)]);
figure(1);