% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function plot_polarisation(time, sig_x, sig_z, titlefig, xlab, zlab, wavepropdirection)
  if(not(exist('wavepropdirection')))
    wavepropdirection='none';
  end
  if(not(exist('xlab')))
    xlab='$v_x$ [m/s]';
  end
  if(not(exist('zlab')))
    zlab='$v_z$ [m/s]';
  end
  if(not(strcmp(wavepropdirection,'none') || strcmp(wavepropdirection,'left') || strcmp(wavepropdirection,'right')))
    error('kek');
  end
  if(not(exist('titlefig')))
    titlefig='';
  end
  
  lw=1;
  
  figure('units','normalized','outerposition',[0 0 0.5 1]);
  
  plot(sig_x, sig_z, 'white', 'linewidth', 0.1*lw);
  hold on;
  
  c=time;
  colormap jet;
%   colormap hsv;
  scatter(sig_x, sig_z, 15, c, 'filled');
  h=colorbar;
  set(h, 'TickLabelInterpreter', 'latex');
  ticckks=unique(sort([h.Ticks, min(time), max(time)]));
  ticckks_lab=split(sprintf('%.0f ',ticckks),' '); ticckks_lab=ticckks_lab(~cellfun('isempty',ticckks_lab));
  set(h, 'Ticks', ticckks); set(h,'ticklabels',ticckks_lab);
  ylabel(h,'time $t$ [s]','interpreter','latex');
  
  daspect([1 1 1]);
  
  minmax=[min(sig_x), max(sig_x), ...
          min(sig_z), max(sig_z)];
  tenPow_x=10.^floor(log10(abs(minmax)));
  minmax_nextint=sign(minmax).*ceil(abs(minmax./tenPow_x)).*tenPow_x;
  
  axis(minmax_nextint);
  set(gca,'Color','k');
  set(gca,'GridColor','white');
  set(gca, 'TickLabelInterpreter', 'latex');
  set(gca,'TickDir','both');
  grid on;
  box on;
  xlabel([xlab]);
  ylabel([zlab]);
  title(titlefig);

  if(not(strcmp(wavepropdirection,'none')))
    % plot propagation direction arrow
    ax=gca();
    yl=ax.YLim;
    xl=ax.XLim;
    percent=5;
    per = (percent/100)*max(diff(yl),diff(xl));
    zarr = (max(yl) - per) * [1, 1];
    xarr = [min(xl), max(xl)] + per*[1,-1];
    P = [xarr;zarr]; % propagating towards right by default
    if(strcmp(wavepropdirection,'left'))
      P = fliplr(P);
    end
    dp=P(:,2)-P(:,1);
    hold on;
    quiver(P(1,1),P(2,1),dp(1),dp(2),0,'color','w','linewidth',2*lw,'maxheadsize',0.5);
    fs=ax.FontSize;
    text(mean(xl), P(2,1) + per*0.6, 'wave propagation direction', 'color', 'w', 'fontsize', 0.9*fs, 'HorizontalAlignment', 'center');
  end
  
  % try to get global rotation trend
  [th, ~] = cart2pol(sig_x,sig_z);
  vecn=vecnorm([sig_x;sig_z]);
  weight=vecn/max(vecn); % weight by the signal's norm
  WGrad=weight.*gradient(th,time); % weighted gradient
%   figure();
%   subplot(3,1,1);plot(time, vecn);
%   subplot(3,1,2);plot(time, th);
%   subplot(3,1,3);plot(time, WGrad);
  gradtheta = mean(WGrad); % derivative of angle wrt time
  disp(['[',mfilename,'] Weighted rotation derivative found to be ',num2str(gradtheta),' rad/s.']);
  
  ax=gca();
  yl=ax.YLim;
  xl=ax.XLim;
  d=0.9*min(diff(xl),diff(yl));
  hold on;
  rectangle('Position',[mean(xl)-d/2, mean(yl)-d/2, d, d],'Curvature',[1,1],'edgecolor','w','linewidth',lw, 'linestyle', '--');
  l=0.075*diff(yl);
  quiver(mean(xl)+d/2, mean(yl)-sign(gradtheta)*l/2, 0, sign(gradtheta)*l, 0,'color','w','linewidth',2*lw,'maxheadsize',50);
  quiver(mean(xl)-d/2, mean(yl)+sign(gradtheta)*l/2, 0, -sign(gradtheta)*l, 0,'color','w','linewidth',2*lw,'maxheadsize',50);
  
end