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

data_t=Ztime;
data_v=Zamp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User input.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fign=-1;
fign=input('  Figure number? > ');
colour=-1;
colour=input('  Colour choice? ([0, 0.447, 0.741] for nice blue, [0.851, 0.325, 0.098] for nice orange, or classical notations) > ');
dist_unit="-1";
while(~ismember(dist_unit,["m", "km"]))
  dist_unit=input('  Distance unit (m, km)? > ', 's');
end
dist_factor=1;
distancechoice=-1;
while(~ismember(distancechoice,[1,2,3,4]))
  distancechoice=input('  Distance choice? (1 for x, 2 for |x|, 3 for z, 4 for d) > ');
end
normalise=-1;
while(~ismember(normalise,[0,1]))
  normalise=input('  Normalise data? (0 for no, 1 for yes) > ');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(dist_unit, "km"))
  dist_factor=1000;
end
dist_unit=strcat(" (",dist_unit,")");
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
distance=distance/dist_factor;
% Sort according to chosen distance.
[~,isort]=sort(distance(istattab(1:nstat)));

% Enventually normalise.
if(normalise==1)
  for i=1:size(isort,1)
    data_v(i,:)=(data_v(i,:)-min(data_v(i,:)))/(max(data_v(i,:))-min(data_v(i,:)));
    data_v(i,:)=data_v(i,:)-mean(data_v(i,:));
  end
end

% Plotting tools.
dist_over_ptp=max(diff(distance(istattab(isort))))/max(peak2peak(data_v(isort,:),2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure.                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(fign);
name={};
yticklabel=[];
vertical_shift={};
scale=-1;
while(scale~=0)
  close(fign); figure(fign);
  if(scale~=0 && scale~=-1)
    dist_over_ptp=scale;
  end
  coef_string=sprintf("%.2e",dist_over_ptp);spls=split(coef_string,"e");coef_string=[char(spls(1)),'\cdot10^{',num2str(str2num(char(spls(2)))),'}'];
  for i=1:size(isort,1)
    istat=isort(i);
    istat_glob = istattab(istat);
%     vertical_shift{istat}=*distance(istat_glob);
    name{istat}=strcat('S', num2str(istat_glob));
    yticklabel=[yticklabel, sprintf("%.2f",distance(istat_glob))];
    plot(data_t(istat,:),distance(istat_glob)+dist_over_ptp*data_v(istat,:),'displayname',name{istat},'color',colour);
    hold on;
  end
  xlim([min(data_t(1:nstat, 1)), max(data_t(1:nstat, end))]);
  xlabel('time (s)');
  ylabel(strcat("$",dist_symbol," + \left(",coef_string,"\right)\times$", unknown_name));
  scale=input('  Rechoose coefficient? (0 for no, new value for yes)? > ');
end

% Eventually label each plot.
labeleach=-1;
while(~ismember(labeleach,[0,1]))
  labeleach=input('  Label each line? (0 for no, 1 for yes) > ');
end
if(labeleach==1)
  yticks(distance(istattab(isort)));
end
% yticklabels(yticklabel);
ylabel(strcat(dist_name," $",dist_symbol,"$ ", dist_unit));%, " $\longrightarrow$"));

% Time axis limits.
tmin=input(['  t_min (',num2str(min(min(data_t))),' now)? > ']); tmax=input(['  t_max (',num2str(max(max(data_t))),' now)? > ']); xlim([tmin, tmax]);

% Stations' names labels.
for i=1:size(isort,1)
  text(1.01*tmax,distance(isort(i)),name{isort(i)},'HorizontalAlignment','left');
end

% Scales.
if(normalise==0)
  scale=-1;
  while(not(ismember(scale,[0,1])))
    scale=input('  Draw scale examples? (0 for no, 1 for yes)? > ');
  end
  if(scale==1)
    min_amplitude_log10=floor(log10(max(max(data_v(:,:)')-min(data_v(:,:)'))));
    for i=1:10
      abs=(tmax+tmin)*i/11;
      plot(abs*[1,1],dist_over_ptp*distance(istattab(isort(1)))*[1,1]+0.5*i*10^(min_amplitude_log10)*[-1,1],'DisplayName',sprintf("%.0e",i*10^(min_amplitude_log10)));
    end
  end
end

% Useful if stations' labels are outside plotting zone.
% ax=gca();
% set(ax, 'Position', [1,1,0.9,1].*ax.Position);