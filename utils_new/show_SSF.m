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

folder = input('  Folder containing the SSF files? > ', 's'); if(not(strcmp(folder(end),'/'))); folder=[folder,'/']; end;
listSSF=dir(strcat(folder, 'SSF*')); listSSF.name, nsources = input('  Number of sources? > ');
gs = input('  Plotting grid size? > ');
ar = input('  Axes aspect ratio (format [ar_x, ar_y, ar_z])? > ');
th = input('  Threshold to recenter plot (show only where SSF>threshold)? > ');

% Load the data.
disp(['Loading.']);
for s=1:nsources
  a=[];
  for i=1:length(listSSF)
    if(not(isempty(regexp(listSSF(i).name, strcat(num2str(s),'_')))))
%       data=importdata(strcat(folder, listSSF(i).name));
      a=unique([a;importdata(strcat(folder, listSSF(i).name))],'rows');
%       a(a(:,3)==0,:)=[]; % Remove lines where value is 0.
      if(mod(i,floor(length(listSSF)/10))==0 || i==length(listSSF)); disp(strcat('Loading data (', num2str(100*i/length(listSSF)), '%).')); end;
    end
  end
  x{s}=a(:,1); z{s}=a(:,2); d{s}=a(:,3);
end

% Interpolate the point cloud.
disp(['Interpolating.']);
for s=1:nsources
  figure(s);
  xx=x{s}; zz=z{s}; dd=d{s};
  F = scatteredInterpolant(x{s},z{s},d{s});
  tx = min(xx(:)):gs:max(xx(:));
  tz = min(zz(:)):gs:max(zz(:));
  [X,Z] = meshgrid(tx,tz);
  D = F(X,Z);
  surf(X,Z,D,'FaceColor', 'interp', 'EdgeColor', 'white', 'LineStyle', ':');
  set(gca,'DataAspectRatio',ar);
  xlabel('x'); ylabel('z'); title(strcat("Source Spatial Function ", num2str(s)));
  xlim([min(xx(d{s}>th)), max(xx(d{s}>th))]); ylim([min(zz(d{s}>th)), max(zz(d{s}>th))]); zlim([th, max(d{s})]);
end