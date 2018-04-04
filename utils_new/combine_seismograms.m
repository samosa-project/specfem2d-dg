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
set(0, 'DefaultLineLineWidth', 2); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 8); % Default at 6.
set(0, 'defaultTextFontSize', 22);
set(0, 'defaultAxesFontSize', 22); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

fign=-1; fign=input('  Figure number? > ');

colour=-1; colour=input('  Colour choice? ([0, 0.447, 0.741] for nice blue, [0.851, 0.325, 0.098] for nice orange, or classical notations) > ');
% colour=[0, 0.447, 0.741];

distancechoice=-1;
while(~ismember(distancechoice,[1,2,3]))
  distancechoice=input('  Distance choice? (1 for x, 2 for |x|, 3 for z, 4 for d) > ');
end
switch distancechoice
  case 1
    distance=xstattab; dist_symbol="x"; dist_name="horizontal distance";
  case 2
    distance=abs(xstattab); dist_symbol="|x|"; dist_name="horizontal distance";
  case 3
    distance=ystattab; dist_symbol="z"; dist_name="altitude";
  case 4
    distance=dist_to_sources; dist_symbol="d"; dist_name="distance";
end

[~,isort]=sort(distance(istattab(1:nstat)));
ptp_over_dist=max(peak2peak(Zamp(isort,:),2))/max(diff(distance(istattab(isort))));

figure(fign);
name={};
yticklabel=[];
vertical_shift={};
scale=-1;
while(scale~=0)
  close(fign); figure(fign);
  if(scale~=0 && scale~=-1)
    ptp_over_dist=scale;
  end
  s=sprintf("%.2e",ptp_over_dist);spls=split(s,"e");s=[char(spls(1)),'\cdot10^{',num2str(str2num(char(spls(2)))),'}'];
  for i=1:size(isort,1)
    istat=isort(i);
    istat_glob = istattab(istat);
    vertical_shift{istat}=ptp_over_dist*distance(istat_glob);
    name{istat}=strcat('S', num2str(istat_glob));
    yticklabel=[yticklabel, sprintf("%.2f",distance(istat_glob))];
    plot(Ztime(istat,:),vertical_shift{istat}+Zamp(istat,:),'displayname',name{istat},'color',colour);
%     text(1.01*Ztime(istat,end),vertical_shift,name{istat},'HorizontalAlignment','left');
    hold on;
  end
  xlim([min(Ztime(1:nstat, 1)), max(Ztime(1:nstat, end))]);
  xlabel('time (s)');
  ylabel(strcat("$\left(",s,"\right)",dist_symbol,"$ + ", unknown_name));
  scale=input('  Rechoose coefficient? (0 for no, new value for yes)? > ');
end
% legend('location', 'best');
yticks(ptp_over_dist*distance(istattab(isort)));
yticklabels(yticklabel);
ylabel(strcat(dist_name," $",dist_symbol,"$ $\longrightarrow$"));

tmin=input('  t_min? > '); tmax=input('  t_max? > '); xlim([tmin, tmax]);
for i=1:size(isort,1)
  text(1.01*tmax,vertical_shift{isort(i)},name{isort(i)},'HorizontalAlignment','left');
end

scale=-1;
while(not(ismember(scale,[0,1])))
  scale=input('  Draw scale examples? (0 for no, 1 for yes)? > ');
end
if(scale==1)
  min_amplitude_log10=floor(log10(max(max(Zamp(:,:)')-min(Zamp(:,:)'))));
  for i=1:10
    abs=(tmax+tmin)*i/11;
    plot(abs*[1,1],ptp_over_dist*distance(istattab(isort(1)))*[1,1]+0.5*i*10^(min_amplitude_log10)*[-1,1],'DisplayName',sprintf("%.0e",i*10^(min_amplitude_log10)));
  end
end

ax=gca();
set(ax, 'Position', [1,1,0.9,1].*ax.Position);