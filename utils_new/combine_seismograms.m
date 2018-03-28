% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         display_seismograms.m should have been ran before.

% clear all;
% close all
clc;
format compact;
set(0, 'DefaultLineLineWidth', 3); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 8); % Default at 6.
set(0, 'defaultTextFontSize', 22);
set(0, 'defaultAxesFontSize', 22); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

colour='b';

distancechoice=-1;
while(~ismember(distancechoice,[1,2,3]))
  distancechoice=input('  Distance choice? (1 for x, 2 for z, 3 for d) > ');
end
switch distancechoice
  case 1
    distance=xstattab;
  case 2
    distance=zstattab;
  case 3
    distance=dist_to_sources;
end

[~,isort]=sort(distance(istattab(1:nstat)));
ptp_over_dist=max(peak2peak(Zamp(isort,:),2))/max(diff(distance(istattab(isort))));

satisfied=-1;
while(satisfied~=0)
  if(satisfied~=0 && satisfied~=-1)
    ptp_over_dist=satisfied;
  end
  s=sprintf("%.2e",ptp_over_dist);spls=split(s,"e");s=[char(spls(1)),'\cdot10^{',num2str(str2num(char(spls(2)))),'}'];
  figure();
  for i=1:size(isort,1)
    istat=isort(i);
    istat_glob = istattab(istat);
    vertical_shift=ptp_over_dist*distance(istat_glob);
    name{istat}=strcat('S', num2str(istat_glob));
    plot(Ztime(istat,:),vertical_shift+Zamp(istat,:),'displayname',name{istat},'color',colour);
    text(1.01*Ztime(istat,end),vertical_shift,name{istat},'HorizontalAlignment','left');
    hold on;
  end
  xlim([min(Ztime(1:nstat, 1)), max(Ztime(1:nstat, end))]);
  xlabel('time (s)');
  ylabel(strcat("$\left(",s,"\right)x_s$ + ", unknown_name));
  satisfied=input('  Rechoose coefficient? (0 for no, new value for yes)? > ');
end
% legend('location', 'best');

tmin=input('  t_min? > '); tmax=input('  t_min? > '); xlim([tmin, tmax]);