% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Combines synthetics under the classical one-panel plot
%                fashion.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_t=Ztime;
data_v=Zamp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ask for user input.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rescale=-1;
while(rescale==-1)
  rescale=input('  Rescale (0 for no, new value for yes)? > ');
end
if(rescale~=0)
  data_v=data_v*rescale;
  disp(['  [WARNING] Data was rescaled by a factor ',num2str(rescale),'.']);
end

filter_data=-1;
% while(not(ismember(filter_data,[0,1,2,3])))
%   filter_data=input('  Filter (0 for no, 1 for high-pass, 2 for low-pass, 3 for band-pass)? > ');
while(not(ismember(filter_data,[0,1])))
  filter_data=input('  Filter (0 for no, 1 for high-pass)? > ');
end
if(filter_data~=0)
  filter_fcutoff=-1;
  while(filter_fcutoff<=0)
    filter_fcutoff=input('  Filter cutoff frequency? > ');
  end
  switch filter_data
    case 1
      for(i=1:nstat)
        [~, data_v_HP] = custom_filter(data_t(i,:), data_v(i,:), filter_fcutoff);
        data_v(i,:)=data_v_HP;
      end
      disp(['  [WARNING] Data was high-pass filtered, with cutoff frequency ',num2str(filter_fcutoff),'.']);
      clear('data_v_HP');
    otherwise
      disp(['  [ERROR] Filtering type not implemented.']);
  end
end
fign=-1;
fign=input('  Figure number? > ');
colour=-1;
colour=input('  Colour choice? ([0, 0.447, 0.741] for nice blue, [0.851, 0.325, 0.098] for nice orange, or classical notations) > ','s');
if(length(colour)>1)
  eval(['colour=',colour]);
end
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

% Remove mean value.
for i=1:nstat
  % Enventually normalise.
  if(normalise==1)
    data_v(i,:)=(data_v(i,:)-min(data_v(i,:)))/(max(data_v(i,:))-min(data_v(i,:)));
  end
  data_v(i,:)=data_v(i,:)-mean(data_v(i,:));
end

% Plotting tools.
if(nstat==1)
  dist_over_ptp=1;
else
  dist_over_ptp=max(diff(distance(istattab(isort))))/max(peak2peak(data_v(isort,:),2));
end
if(dist_over_ptp>1e15)
  error("Variable dist_over_ptp is > 1e15, probably coming from the signal being very small everywhere.");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure.                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(fign);
name={};
yticklabel=[];
scale=-1;
while(scale~=0)
  close(fign);
  figure(fign);
  if(scale~=0 && scale~=-1)
    dist_over_ptp=scale;
  end
  % Prepare a string to give information about what is plotted.
  coef_string=sprintf("%.2e",dist_over_ptp);
  spls=split(coef_string,"e");
  coef_string=[char(spls(1)),'\cdot10^{',num2str(str2num(char(spls(2)))),'}'];
  clear('spls');
  for i=1:nstat
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
  scale=input('  Rechoose coefficient (0 for no, new value for yes)? > ');
end

ylabel(strcat(dist_name," $",dist_symbol,"$ ", dist_unit, ' $\longrightarrow$'));%, " $\longrightarrow$"));
set(gca, 'TickLabelInterpreter','latex');

% Time axis limits.
tmin=input(['  t_min (',num2str(min(min(data_t))),' now)? > ']);
tmax=input(['  t_max (',num2str(max(max(data_t))),' now)? > ']);
if(isempty(tmin))
  tmin=min(min(data_t));
end
if(isempty(tmax))
  tmax=max(max(data_t));
end
xlim([tmin, tmax]);

% Eventually label each plot.
labeleach=-1;
while(~ismember(labeleach,[0,1]))
  labeleach=input('  Label each line? (0 for no, 1 for yes) > ');
end
if(labeleach==1)
  % YTicks.
  yticks(distance(istattab(isort)));
  % Stations' names labels.
  for i=1:size(isort,1)
    text(1.01*tmax,distance(istattab(isort(i))),name{isort(i)},'HorizontalAlignment','left');
  end
end

% Scales.
if(normalise==0)
  scale=-1;
  while(not(ismember(scale,[0,1])))
    scale=input('  Draw scale examples? (0 for no, 1 for yes)? > ');
  end
  if(scale==1)
    min_amplitude_log10=floor(log10(max(max(data_v(:,:)')-min(data_v(:,:)'))));
    curxlim=get(gca,'xlim');
    for i=1:10
      absc=curxlim(1)+diff(curxlim)*i/11;
      plot(absc*[1,1],distance(istattab(isort(1)))*[1,1]+dist_over_ptp*0.5*i*10^(min_amplitude_log10)*[-1,1],'DisplayName',sprintf("%.0e",i*10^(min_amplitude_log10)));
    end
  end
end

% Useful if stations' labels are outside plotting zone.
% ax=gca();
% set(ax, 'Position', [1,1,0.9,1].*ax.Position);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear variables.             %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('coef_string', 'colour', 'data_t', 'data_v', 'dist_factor', 'dist_name', 'dist_over_ptp', 'distancechoice', 'fign', 'i', 'isort', 'istat', 'istat_glob', 'labeleach', 'name', 'normalise', 'scale', 'spls', 'tmax', 'tmin', 'yticklabel');