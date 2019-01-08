% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Combines time series under the classical one-panel plot
%                fashion, with some distance as abscissas.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [] = plot_time_v_dist(Ztime,Zamp,distance)

  % clear all;
  % close all;
  % clc;
  format compact;
  set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
  set(0, 'defaultTextFontSize', 12); set(0, 'defaultAxesFontSize', 12);
  set(0, 'DefaultTextInterpreter', 'latex');
  set(0, 'DefaultLegendInterpreter', 'latex');
  
  % Make sure data has right shape.
  % The following lines assume we have more time steps than stations to plot.
  Ztime=reshape(Ztime,[min(size(Ztime)), max(size(Ztime))]);
  Zamp=reshape(Zamp,[min(size(Zamp)), max(size(Zamp))]);
  if(not(all(size(Ztime)==size(Zamp))))
    error(['[',mfilename,', ERROR] time data and amplitude should have the same size, but right now do not.']);
  end
  nbstat=size(Ztime,1);
  if(not(numel(distance)==nbstat))
    error(['[',mfilename,', ERROR] distance array should contain the same number of stations as the number of data series, but right now do not.']);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Load.                       %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  data_t = Ztime;
  data_v = Zamp;
  unknown_name = 'v_z';

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

  filter_data = - 1;
  % while(not(ismember(filter_data,[0,1,2,3])))
  %   filter_data=input(['  Filter (0 for no, 1 for high-pass, 2 for low-pass, 3 for band-pass)? > ']);
  while (not(ismember(filter_data, [0, 1, 2, 3])))
    filter_data = input(['[', mfilename, '] Filter (0 for no, 1 for high-pass, 2 for low-pass, 3 for band-pass)? > ']);
  end
  if (filter_data ~= 0)
    filter_fchp = -1;
    filter_fclp = -1;
    filter_dhp = 0.65;
    filter_dlp = 0.65;
    filter_ord = 2;
    if(ismember(filter_data, [1, 3]))
      while (filter_fchp <= 0)
        filter_fchp = input(['[', mfilename, '] High-pass cutoff frequency? > ']);
      end
    end
    if(ismember(filter_data, [2, 3]))
      while (filter_fclp <= 0)
        filter_fclp = input(['[', mfilename, '] Low-pass cutoff frequency? > ']);
      end
    end
    
    for i = 1:nbstat
      switch filter_data
        case 1
          [filt_f,~,filt_hhp,~]=custom_Filter(mean(1./diff(data_t(i, :))),length(data_t(i, :)),filter_fclp,filter_fchp,filter_dlp,filter_dhp,filter_ord);
          data_v(i, :) = real(ifft(fft(detrend(data_v(i, :))) .* fftshift(filt_hhp)));
          disp(['[',mfilename,'] Highpass filtered data over ',num2str([filter_fchp]),' Hz with homemade order ',num2str(filter_ord),' filter.']);
        case 2
          [filt_f,filt_hlp,~,~]=custom_Filter(mean(1./diff(data_t(i, :))),length(data_t(i, :)),filter_fclp,filter_fchp,filter_dlp,filter_dhp,filter_ord);
          data_v(i, :) = real(ifft(fft(detrend(data_v(i, :))) .* fftshift(filt_hlp)));
          disp(['[',mfilename,'] Lowpass filtered data under ',num2str([filter_fclp]),' Hz with homemade order ',num2str(filter_ord),' filter.']);
        case 3
%           [~, data_v_HP] = custom_filter(data_t(i, :), data_v(i, :), filter_fcutoff);
          [filt_f,~,~,filt_hbp]=custom_Filter(mean(1./diff(data_t(i, :))),length(data_t(i, :)),filter_fclp,filter_fchp,filter_dlp,filter_dhp,filter_ord);
          data_v(i, :) = real(ifft(fft(detrend(data_v(i, :))) .* fftshift(filt_hbp)));
%           data_v(i, :) = data_v_HP;
          disp(['[',mfilename,'] Bandpass filtered data between (',num2str([filter_fchp,filter_fclp]),') Hz with homemade order ',num2str(filter_ord),' filter.']);
%           disp(['[', mfilename, ', WARNING] Data was high-pass filtered, with cutoff frequency ', num2str(filter_fchp), '.']);
%           clear('data_v_HP');
        otherwise
          error(['[', mfilename, ', ERROR] Filtering type not implemented.']);
      end
    end
  end
  fign = - 1;
  fign = input(['[', mfilename, '] Figure number? > ']);
%   colour = - 1;
%   colour = input(['[', mfilename, '] Colour choice? ([0, 0.447, 0.741] for nice blue, [0.851, 0.325, 0.098] for nice orange, or classical notations) > '], 's');
%   if (length(colour) > 1)
%     eval(['colour=', colour]);
%   end
  colour='k';
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
  dist_unit = strcat(" (", dist_unit, ")");
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
  dist_symbol='d';
  dist_name = 'distance';
  distance = distance / dist_factor;
  % Sort according to chosen distance.
%   [~, isort] = sort(distance(istattab(1:nstat)));
  [~, isort] = sort(distance);

  % Remove mean value.
  for i = 1:nbstat
    gap2detrend = abs(data_v(i, :)-detrend(data_v(i, :)));
    maxpercentdetrend=100*[max(gap2detrend)]/(max(data_v(i, :)) - min(data_v(i, :)));
    if(max(maxpercentdetrend)>5)
      disp(['[',mfilename,'] For dataset n°',num2str(i),', detrend would shift data values by a quantity which is ',sprintf('%.2g',maxpercentdetrend),' % of signal amplitude. Discarding detrend.'])
    else
      data_v(i, :) = detrend(data_v(i, :));
    end
%     data_v(i, :) = data_v(i, :) - mean(data_v(i, :));

    % Eventually normalise.
    if (normalise == 1)
      data_v(i, :) = (data_v(i, :) - min(data_v(i, :))) / (max(data_v(i, :)) - min(data_v(i, :)));
    end
  end

  % Plotting tools.
  if (nbstat == 1)
    dist_over_ptp = 1;
  else
    %dist_over_ptp = max(diff(distance(istattab(isort)))) / max(peak2peak(data_v(isort, :), 2));
    dist_over_ptp = max(diff(distance(isort))) / max(peak2peak(data_v(isort, :), 2));
  end
  if (dist_over_ptp > 1e15)
    error(['[', mfilename, ', ERROR] Variable dist_over_ptp is > 1e15, probably coming from the signal being very small everywhere for one of the signals.']);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Figure.                     %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(fign);
  name = {};
  yticklabel = [];
  scalez = - 1;
  while (scalez ~= 0)
    close(fign);
    figure(fign);
    if (scalez ~= 0 && scalez ~= - 1)
      dist_over_ptp = scalez;
    end
    % Prepare a string to give information about what is plotted.
    coef_string = sprintf(" %.2e",dist_over_ptp);
    spls = split(coef_string, "e");
    coef_string = [char(spls(1)), '\cdot10^{', num2str(str2num(char(spls(2)))), '}'];
    clear('spls');
    for i = 1:nbstat
      istat = isort(i);
%       istat_glob = istattab(istat);
      istat_glob = istat;
      %     vertical_shift{istat}=*distance(istat_glob);
      name{istat} = strcat('S', num2str(istat_glob));
      yticklabel = [yticklabel, sprintf(" %.2f",distance(istat_glob))];
      plot(data_t(istat, :), distance(istat_glob) + dist_over_ptp * data_v(istat, :), 'displayname', name{istat}, 'color', colour);
      hold on;
    end
    xlim([min(data_t(1:nbstat, 1)), max(data_t(1:nbstat, end))]);
    xlabel('time (s)');
    ylabel(strcat("$", dist_symbol, " + \left(", coef_string, "\right)\times ", unknown_name,'$'));
    scalez = input(['[', mfilename, '] Rechoose coefficient (0 for no, new value for yes)? > ']);
  end

  ylabel(strcat(dist_name, " $", dist_symbol, "$ ", dist_unit)); %, " $\longrightarrow$"));
  set(gca, 'TickLabelInterpreter', 'latex');

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
    yticks(distance);
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
        plot(absc * [1, 1], distance(isort(1)) * [1, 1] + dist_over_ptp * 0.5 * scaalz(i) * [- 1, 1], 'DisplayName', sprintf(" %.0e",scaalz(i)));
      end
    end
  end

  % Useful if stations' labels are outside plotting zone.
  % ax=gca();
  % set(ax, 'Position', [1,1,0.9,1].*ax.Position);
  grid on;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Clear variables.             %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  clear('coef_string', 'colour', 'data_t', 'data_v', 'dist_factor', 'dist_name', 'dist_over_ptp', 'distancechoice', 'fign', 'i', 'isort', 'istat', 'istat_glob', 'labeleach', 'name', 'normalise', 'scale', 'spls', 'tmax', 'tmin', 'yticklabel');
end

