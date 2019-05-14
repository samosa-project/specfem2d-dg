% Author:        Léo Martire.
% Description:   Plot a bunch of data into a single figure.
% Notes:         TODO.
%
% Usage:
%   plotOneByOne(Ztime, Zamp, istattab, xstat, zstat, norm_ylims)
% with:
%   TODO.
% yields:
%   TODO.

function fh=plotOneByOne(Ztime, Zamp, istattab, xstat, zstat, norm_ylims, figtitle, ylab, existingfignumber)
  if(not(exist('figtitle')))
    figtitle='';
  end
  if(not(exist('ylab')))
    ylab='SIGNAL';
  end
  if(not(exist('existingfignumber')))
    existingfignumber=-1;
  end
  
  if(not(all(size(Ztime)==size(Zamp))))
    error('must have same sizes');
  end
  nstat = size(Ztime, 1);
  
  if(existingfignumber>0)
    % if a figure is provided, plot on it
    fh=figure(existingfignumber);
    % retrieve axes, by checking children and removing everything which are not axes
    axxx=fh.Children;
    itodel=[];
    for i=1:numel(axxx)
      if(not(strcmp(axxx(i).Type,'axes')))
        itodel=[itodel,i];
      end
    end
    axxx(itodel)=[];
    axxx=flip(axxx);
  else
    % else, create new figure
    fh=figure('units','normalized','outerposition',[0 0 1 1]);
    [axxx] = tight_subplot(nstat, 1, 0.005, 0.05, [0.05, 0.01]);
  end
  
%   axes = {};
  for istat = 1:nstat
    istat_glob = istattab(istat); % Recover global number of station.
%     subplot(nstat, 1, istat);
%     axes{istat} = gca;
    axes(axxx(istat));
    legtext{istat} = ['S', num2str(istat_glob), ', $(x,z) = (', num2str(xstat(istat)), ',', num2str(zstat(istat)), ')$ [m]'];
    plot(Ztime(istat, :), Zamp(istat, :), 'displayname', legtext{istat}); hold on;
    % Cosmetics.
    if (istat == 1); title(figtitle); end; % Put title only on first subplot.
    if (istat == nstat); xlabel('time [s]'); end; % Put xlabel only on last subplot.
    if (istat ~= nstat); set(gca, 'xticklabel', []); end; % Remove xticks for all subplots except last one.
    if (istat == round(nstat / 2)); ylabel(ylab); end; % Put one ylabel at the middle subplot.
    xlim([min(min(Ztime)), max(max(Ztime))]);
    legend('Location', 'northeast');
    set(gca, 'TickLabelInterpreter', 'latex'); grid on; box on;
  end
%   axess=[];
%   for i=1:numel(axxx); axess(i)=axxx{i}; end;
  if (norm_ylims)
    linkaxes(axxx); % Link both x and y.
  else
    linkaxes(axxx, 'x'); % Link only x.
  end
end