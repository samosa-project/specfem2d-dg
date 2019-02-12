% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Simpler and more lax version of plot_time_v_dist. With an
%                accent on usability with periodised computation domains.
%                Cannot be really used 'in production' as is.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [] = plot_periodised_time_v_dist(inp_t, inp_v, inp_d, maxvelocity, duration, d_farthest2closest, nreps, scalez)
  
  if(nargin~=8)
    error(['[',mfilename,', ERROR] Not enough input arguments. Needs ''data_time, data_vals, distance, group_velocity, duration, d_farthest2closest, n_repetitions, plot_scaling''.']);
  end
  
  inp_t=reshape(inp_t,[min(size(inp_t)), max(size(inp_t))]);
  inp_v=reshape(inp_v,[min(size(inp_v)), max(size(inp_v))]);
  reshape(inp_d,[numel(inp_d),1]);

%   nreps=2;
  propagatingTo='right';
  
  Ztimeperio=repmat(inp_t,[nreps,1]);
  Zampperio=repmat(inp_v,[nreps,1]);
  
%   d_farthest2closest=10; % distance from farthest station to closest station w.r.t. periodisation
  
  switch propagatingTo
    case 'right'
      shift=-min(inp_d)+d_farthest2closest+max(inp_d);
      distperiodised=repmat(inp_d,[nreps,1]) + sort(repmat(shift * ((1:nreps)-1),[1,numel(inp_d)]))'; % ugly af
%       distperiodised=repmat(inp_d,[nreps,1]) + sort(repmat((d_farthest2closest+min(inp_d)) * ((1:nreps)-1),[1,numel(inp_d)]))'; % ugly af
    otherwise
      error();
  end
  
%   maxvelocity=60/(.1784-.00656); % speed of fastest signal
%   maxvelocity=130/(.3843-.00656); % speed of fastest signal
  
%   duration=.01952+.01576; % duration of signal to observe
  
  distperiodised;
  windowt=min(min(inp_t))+distperiodised/maxvelocity;
  
  sel=(Ztimeperio > windowt) & (Ztimeperio<windowt+duration);
  [windowt, windowt+duration];

  colour='k';
  
  data_t={};
  data_v={};
  maxp2p=-Inf;
  for istat = 1:size(Ztimeperio,1)
    data_t{istat}=Ztimeperio(istat,sel(istat,:));
    data_v{istat}=Zampperio(istat,sel(istat,:));
    if(length(data_v{istat})~=0)
      maxp2p=max(maxp2p,max(data_v{istat}));
    end
  end
  if(maxp2p==-Inf)
    error('kek');
  end
  dist_over_ptp=scalez * max(diff(distperiodised))/maxp2p;
  
  figure();
  hold on;
  for istat = 1:size(Ztimeperio,1)
    name{istat}=['$d=',num2str(distperiodised(istat)),'$'];
    plot(Ztimeperio(istat,:), distperiodised(istat) + dist_over_ptp * Zampperio(istat,:),'color',colour,'linestyle',':','linewidth',.5, 'displayname', name{istat});
    hold on;
  end
  for istat = 1:size(Ztimeperio,1)
    plot(data_t{istat}, distperiodised(istat) + dist_over_ptp * data_v{istat},'color',colour, 'displayname', name{istat});
    hold on;
  end
  xlim([min(inp_t(:, 1)), max(inp_t(:, end))]);
  ylim([min(distperiodised)+min(dist_over_ptp * Zampperio(1,:)),max(distperiodised)]);
  grid on;
  set(gca, 'TickLabelInterpreter', 'latex');
  xlabel('time [s]');
  ylabel(['$x$ [m]']);
  yticks(distperiodised);
end

