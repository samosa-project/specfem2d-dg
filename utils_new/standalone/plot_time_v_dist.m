% Author:        Léo Martire.
% Description:   Combines time series under the classical one-panel plot
%                fashion, with some distance as abscissas.
% Notes:         Needs:
%                a) .m scripts and functions (if not alongside this script, recover via Léo):
%                  1) bulkfilter.m
%
% Usage:
%   plot_time_v_dist(times, values, distances, reducedTime, figureTitle, distanceName, skipFilters)
% with:
%   TODO.
% yields:
%   TODO.

function fh = plot_time_v_dist(times, values, distances, reducedTime, figureTitle, distanceName, skipFilters, names, colours, linestyles)
  if(nargin<3)
    error(['[',mfilename,', ERROR] Not enough input arguments. Needs ''times,values,distances''.']);
  end
  if(not(exist('reducedTime')))
    % not found, make default
    reducedTime=0;
  end
  if(not(exist('figureTitle')))
    % not found, make default
    figureTitle='';
  end
  if(not(exist('distanceName')))
    distanceName=['distance $d$ '];
  end
  if(not(exist('skipFilters')))
    % not found, make default
    skipFilters=0;
  end
  if(not(exist('names')))
    namesGiven=0;
  else
    if(numel(names)==1 && names==0)
      namesGiven=0;
    else
      namesGiven=1;
    end
  end
  
  if(not(exist('colours')))
    coloursGiven=0;
  else
    if(numel(colours)==1 && colours==0)
      coloursGiven=0;
    else
      coloursGiven=1;
    end
  end
  
  if(not(exist('linestyles')))
    linestylesGiven=0;
  else
    if(numel(linestyles)==1 && linestyles==0)
      linestylesGiven=0;
    else
      linestylesGiven=1;
    end
  end
  
%   format compact;
%   set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
%   set(0, 'defaultTextFontSize', 12); set(0, 'defaultAxesFontSize', 12);
%   set(0, 'DefaultTextInterpreter', 'latex');
%   set(0, 'DefaultLegendInterpreter', 'latex');
  
  % Make sure data has right shape.
  % The following lines assume we have more time steps than stations to plot.
  times  = reshape(times, [min(size(times)), max(size(times))]);
  values = reshape(values, [min(size(values)), max(size(values))]);
  if(not(all(size(times)==size(values))))
    error(['[',mfilename,', ERROR] time data and amplitude should have the same size, but right now do not.']);
  end
  nbstat = size(times, 1);
  if(not(numel(distances)==nbstat))
    error(['[',mfilename,', ERROR] distance array should contain the same number of stations as the number of data series, but right now do not.']);
  end
  distances = reshape(distances, [nbstat, 1]);
  if(namesGiven)
    if(not(numel(names)==nbstat))
      error(['[',mfilename,', ERROR] If a name structure is given, it must contain same number of names as the number of data series.']);
    end
  end
  if(coloursGiven)
    if(not(numel(colours)==nbstat))
      error(['[',mfilename,', ERROR] If a colours structure is given, it must contain same number of names as the number of data series.']);
    end
  end
  if(linestylesGiven)
    if(not(numel(linestyles)==nbstat))
      error(['[',mfilename,', ERROR] If a linestyles structure is given, it must contain same number of names as the number of data series.']);
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Load.                       %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  data_t = times;
  data_v = values;
  unknown_name = 'v';
  xlab = 'time $t$ [s]';
  
  if(reducedTime~=0)
    data_t = data_t - distances/reducedTime;
    xlab = ['reduced time $t_r = t-d/',num2str(reducedTime),'$ [s]'];
  end
  
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Parse ignored data.         %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   notignored_IDs={};
%   for i = 1:nbstat
%     % ignore NANs.
% %     notignored_IDs{i} = find( not(isnan(data_v(i, :))) );
%     diftims=[Inf,diff(data_t(i,:))]; % store differential time, put Inf at the very beginning to always consider the very first one
%     notignored_IDs{i}=find(diftims>0); % ignore where differential time is negative or 0
%     % ie, if times=[-1,0,1,2,3,0,0,0], ignore the trailing [0,0,0] and only consider the [-1,0,1,2,3]
%   end
%   ttt=data_t(1,notignored_IDs{1}); [min(ttt), max(ttt)]
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Ask for user input.         %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   rescale = - 1;
%   while (rescale == - 1)
%     rescale = input(['[', mfilename, '] Rescale data (0 for no, new value for yes)? > ']);
%   end
%   if (rescale ~= 0)
%     data_v = data_v * rescale;
%     disp(['[', mfilename, ', WARNING] Data was rescaled by a factor ', num2str(rescale), '.']);
%   end
  
  if(not(skipFilters))
    data_v = bulkfilter(data_t, data_v);
  end
  
  fign = - 1;
  fign = input(['[', mfilename, '] Figure number? > ']);
%   colour = - 1;
%   colour = input(['[', mfilename, '] Colour choice? ([0, 0.447, 0.741] for nice blue, [0.851, 0.325, 0.098] for nice orange, or classical notations) > '], 's');
%   if (length(colour) > 1)
%     eval(['colour=', colour]);
%   end
%   colour='k';
%   dist_unit = "-1";
%   while (~ ismember(dist_unit, ["m", "km"]))
%     dist_unit = input(['[', mfilename, '] Distance unit (m, km)? > '], 's');
%   end
  dist_unit='m';
  dist_factor = 1;
%   distancechoice = - 1;
%   while (~ ismember(distancechoice, [1, 2, 3, 4]))
%     distancechoice = input(['[', mfilename, '] Distance choice? (1 for x, 2 for |x|, 3 for z, 4 for d) > ']);
%   end
  normalise = - 1;
  while (~ ismember(normalise, [0, 1]))
    normalise = input(['[', mfilename, '] Normalise all data? (0 for no, 1 for yes) > ']);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Treatment.                  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (strcmp(dist_unit, "km"))
    dist_factor = 1000;
  end
  dist_unit = [' [', dist_unit, ']'];
%   switch distancechoice
%     case 1
%       distance = xstattab; dist_symbol = "x"; dist_name = "horizontal distance";
%     case 2
%       distance = abs(xstattab); dist_symbol = "|x|"; dist_name = "horizontal distance";
%     case 3
%       distance = ystattab; dist_symbol = "z"; dist_name = "altitude";
%     case 4
%       distance = dist_to_sources; dist_symbol = "d"; dist_name = "distance";
%   end
%   dist_symbol='d';
%   dist_name = 'distance';
  distances = distances / dist_factor;
  % Sort according to chosen distance.
%   [~, isort] = sort(distance(istattab(1:nstat)));
  [~, isort] = sort(distances);

  % Remove mean value.
  for i = 1:nbstat
%     sel=notignored_IDs{i};
    gap2detrend = abs(data_v(i, :)-detrend(data_v(i, :)));
    maxpercentdetrend=100*[max(gap2detrend)]/(max(data_v(i, :)) - min(data_v(i, :))); % finds maximum gap ignoring eventual NaNs.
%     gap2detrend = abs(data_v(i, sel)-detrend(data_v(i, sel)));
%     maxpercentdetrend=100*[max(gap2detrend)]/(max(data_v(i, sel)) - min(data_v(i, sel))); % finds maximum gap ignoring eventual NaNs.
    if(max(maxpercentdetrend)>5)
      disp(['[',mfilename,'] For dataset n°',num2str(i),', detrend would shift data values by a quantity which is ',sprintf('%.2g',maxpercentdetrend),' % of signal amplitude. Discarding detrend.'])
    else
      data_v(i, :) = detrend(data_v(i, :));
%       data_v(i, sel) = detrend(data_v(i, sel)); % detrend only where data is not ignored
    end
%     data_v(i, :) = data_v(i, :) - mean(data_v(i, :));

    % Eventually normalise.
    if (normalise == 1)
      data_v(i, :) = data_v(i, :) / (max(data_v(i, :)) - min(data_v(i, :))); % Transform into signal with 1 p2p amplitude (without any shift).
%       data_v(i, sel) = data_v(i, sel) / (max(data_v(i, sel)) - min(data_v(i, sel))); % Transform into signal with 1 p2p amplitude (without any shift).
    end
  end

  % Plotting tools.
  if (nbstat == 1)
    dist_over_ptp = 1;
  else
    %dist_over_ptp = max(diff(distance(istattab(isort)))) / max(peak2peak(data_v(isort, :), 2));
    dist_over_ptp = max(diff(distances(isort))) / max(peak2peak(data_v(isort, :), 2)); % peaktopeak ignores NaN
  end
  if (dist_over_ptp > 1e15)
    error(['[', mfilename, ', ERROR] Variable dist_over_ptp is > 1e15, probably coming from the signal being very small everywhere for one of the signals.']);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Figure.                     %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fh = figure(fign);
%   name = {};
  yticklabel = [];
  scalez = - 1;
  while (scalez ~= 0)
    close(fign);
    fh = figure(fign);
    set(fh,'units','normalized','outerposition',[0 0 1 1])
    if (scalez ~= 0 && scalez ~= - 1)
      dist_over_ptp = scalez;
    end
    % Prepare a string to give information about what is plotted.
%     dist_over_ptp
    coef_string = sprintf(" %.2e",dist_over_ptp);
    spls = split(coef_string, "e");
    coef_string = [char(spls(1)), '\cdot10^{', num2str(str2num(char(spls(2)))), '}'];
    clear('spls');
    for i = 1:nbstat
%       sel = notignored_IDs{i};
      istat = isort(i);
%       istat_glob = istattab(istat);
      istat_glob = istat;
      %     vertical_shift{istat}=*distance(istat_glob);
      
      if(namesGiven)
        name = names{istat};
      else
        name = strcat('S', num2str(istat_glob));
      end
      
      if(coloursGiven)
        colour=colours{istat};
      else
        colour='k';
      end
      
      if(linestylesGiven)
        LS = linestyles{istat};
      else
        LS = '-';
      end
      
      yticklabel = [yticklabel, sprintf(" %.2f",distances(istat_glob))];
      
%       [i, dist_over_ptp, min(dist_over_ptp * data_v(istat, sel)), max(dist_over_ptp * data_v(istat, sel))]
      
%       plot(data_t(istat, sel), distances(istat_glob) + dist_over_ptp * data_v(istat, sel), 'displayname', name{istat}, 'color', colour);
      plot(data_t(istat, :), distances(istat_glob) + dist_over_ptp * data_v(istat, :), 'displayname', name, 'color', colour,'linestyle',LS);
      hold on;
    end
%     xlim([min(data_t(:, 1)), max(data_t(:, end))]); % max ignores nan
    xlim([min(min(data_t(:, :))), max(max(data_t(:, :)))]); % max ignores nan
%     xlim([min(min(data_t(:, sel))), max(max(data_t(:, sel)))]); % max ignores nan
    xlabel(xlab);
    ylabel(['$d + \left(', coef_string, '\right)\times ', unknown_name,'$']);
    scalez = input(['[', mfilename, '] Rechoose coefficient (0 for no, new value for yes)? > ']);
  end

  ylabel(strcat(distanceName, dist_unit)); %, " $\longrightarrow$"));
  set(gca, 'TickLabelInterpreter', 'latex');
  title(figureTitle);

  % Time axis limits.
  tmin = input(['[', mfilename, '] t_min (', num2str(min(min(data_t))), ' now)? > ']);
  tmax = input(['[', mfilename, '] t_max (', num2str(max(max(data_t))), ' now)? > ']);
  if (isempty(tmin))
    tmin = min(min(data_t));
  end
  if (isempty(tmax))
    tmax = max(max(data_t));
  end
  xlim([tmin, tmax]);

  % Eventually label each plot.
  labeleach = - 1;
  while (~ ismember(labeleach, [0, 1]))
    labeleach = input(['[', mfilename, '] Label each line with distance? (0 for no, 1 for yes) > ']);
  end
  if (labeleach == 1)
    % YTicks.
    yticks(unique(distances(isort)));
%     % Stations' names labels.
%     for i = 1:size(isort, 1)
%       text(1.01 * tmax, distance(istattab(isort(i))), name{isort(i)}, 'HorizontalAlignment', 'left');
%     end
  end

  % Scales.
  if (normalise == 0)
    scalez = - 1;
    while (not(ismember(scalez, [0, 1])))
      scalez = input(['[', mfilename, '] Draw scale examples? (0 for no, 1 for yes)? > ']);
    end
    if (scalez == 1)
      maxamp=max(max(data_v(:, :)')-min(data_v(:,:)'));
      maxpower=floor(log10(maxamp));
      maxscale=floor(maxamp/(10^maxpower))+1;
      scaa=linspace(maxscale-9,maxscale,10);
      scaa(scaa==0)=[];
      if(length(scaa)<10)
        scaa=[scaa(1)-1,scaa];
      end
      scaalz=((scaa<=0).*(10+scaa)+(scaa>0).*scaa).*((scaa<=0)*10^(maxpower-1)+(scaa>0)*10^maxpower);
%       min_amplitude_log10 = floor(log10(max(max(data_v(:, :)')-min(data_v(:,:)'))));
      curxlim = get(gca, 'xlim');
      for i = 1:10
        absc = curxlim(1) + diff(curxlim) * i / 11;
%         plot(absc * [1, 1], distance(isort(1)) * [1, 1] + dist_over_ptp * 0.5 * i * 10 ^ (min_amplitude_log10) * [- 1, 1], 'DisplayName', sprintf(" %.0e",i*10^(min_amplitude_log10)));
        plot(absc * [1, 1], distances(isort(1)) * [1, 1] + dist_over_ptp * 0.5 * scaalz(i) * [- 1, 1], 'DisplayName', sprintf(" %.0e",scaalz(i)));
      end
    end
  end

  % Useful if stations' labels are outside plotting zone.
  % ax=gca();
  % set(ax, 'Position', [1,1,0.9,1].*ax.Position);
  grid on;
end

