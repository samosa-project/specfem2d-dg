% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   plot_polarisation(time, sig_x, sig_y, lab_x, lab_y, titlefig, wavepropdirection)
% with:
%   TODO.
% yields:
%   TODO.

% One-liner to produce test data:
%   noi=0.75; rota=pi/6; t=0:.1:8*pi;v=[cos(t);2*sin(t)]'+noi*rand(numel(t),2);R=[cos(-rota),sin(-rota);-sin(-rota),cos(-rota)]; for i=1:numel(t);v(i,:)=(R*v(i,:)')';end;
%   plot_polarisation(t, detrend(v(:,1)), detrend(v(:,2)));

function plot_polarisation(time, sig_x, sig_y, lab_x, lab_y, titlefig, wavepropdirection)
  if(nargin<3)
    error(['[',mfilename,', ERROR] Must provide at least (time, sig_x, sig_y). Run ''help ',mfilename,'''.']);
  end
  if(not(exist('titlefig')))
    titlefig='';
  end
  if(not(exist('lab_x')))
    lab_x='$v_x$ [m/s]';
  end
  if(not(exist('lab_y')))
    lab_y='$v_z$ [m/s]';
  end
  if(not(exist('wavepropdirection')))
    wavepropdirection='none';
  end
  if(not(strcmp(wavepropdirection,'none') || strcmp(wavepropdirection,'left') || strcmp(wavepropdirection,'right')))
    error('kek');
  end
  time  = reshape(time,  numel(time),  1);
  sig_x = reshape(sig_x, numel(sig_x), 1);
  sig_y = reshape(sig_y, numel(sig_y), 1);
  if(not(numel(time)==numel(sig_x) & numel(sig_x)==numel(sig_y)))
    error(['[',mfilename,', ERROR] (time, sig_x, sig_y) must be the same length. Run ''help ',mfilename,'''.']);
  end
  
  doWindowSignal = -1;
  while (not((numel(doWindowSignal)==1 && doWindowSignal==0) || (numel(doWindowSignal)==2)))
    doWindowSignal = input(['[', mfilename, '] Window signal ([t1 t2] vector for yes, 0 for no)? > ']);
  end
  if(any(doWindowSignal))
    sel = (time>=min(doWindowSignal) & time<=max(doWindowSignal));
    time = time(sel);
    sig_x = sig_x(sel);
    sig_y = sig_y(sel);
  end
  
  overdrawcolor='k';
  
  % Parameters for plotting and displaying.
  lw=1;
  frmt_eigang_cnsl = '%6.1f';
  frmt_eigang_plot = '%.1f';
  frmt_eigval_cnsl = '%6.1e';
  frmt_eigval_plot = '%.1e';
  
  figure('units','normalized','outerposition',[0 0 0.5 1]);
%   tight_subplot(1, 1, 0.005, 0.05, [0.05, 0.01]); % not mandatory, but prettier % incompatible with daspect
  
  % Get and plot major/minor axes (under everything else).
  [eigval, eigvec] = simple2DPCA([sig_x,sig_y]);
  [th_eig,~] = cart2pol(eigvec(1,:),eigvec(2,:)); % send x-values and y-values, this is why it's eigvec(1,:) and not actually the first eigenvector eigvec(:,1) etc.
  
  [~,idsorted]=sort(eigval);
  ilmax=idsorted(end);
  ilmin=idsorted(1);
  
%   line(1e9*[-evec1(1),evec1(1)],1e9*[-evec1(2),evec1(2)])
  factor=2*max(max(abs(sig_x)),max(abs(sig_y)));
  i=1; eigVec1LineP=factor*eigvec(:,i)*[-1,1]; line(eigVec1LineP(1,:),eigVec1LineP(2,:),'color',overdrawcolor,'linewidth',0.1*lw, 'linestyle', '-'); hold on;
  i=2; eigVec2LineP=factor*eigvec(:,i)*[-1,1]; line(eigVec2LineP(1,:),eigVec2LineP(2,:),'color',overdrawcolor,'linewidth',0.1*lw, 'linestyle', '-'); hold on;
               disp(['[',mfilename,'] Ellipticity: great axis (l =',sprintf(frmt_eigval_cnsl,eigval(ilmax)),') is ',sprintf(frmt_eigang_cnsl,th_eig(ilmax)*180/pi),'° and']);
  disp([blanks(length(mfilename)+2),'              small axis (l =',sprintf(frmt_eigval_cnsl,eigval(ilmin)),') is ',sprintf(frmt_eigang_cnsl,th_eig(ilmin)*180/pi),'° counterclockwise from positive-pointing abscissas.']);
  
  
  % Actual plot of polarisation.
  plot(sig_x, sig_y, 'white', 'linewidth', 0.1*lw, 'linestyle',':'); hold on
  scatter(sig_x, sig_y, 15, time, 'filled'); hold on;
  daspect([1 1 1]); % incompatible with tight_subplot
  
  % Add a colorbar for time.
  colormap jet;
%   colormap hsv;
  h=colorbar;
  set(h, 'TickLabelInterpreter', 'latex');
  ticckks=unique(sort([h.Ticks, min(time), max(time)]));
  ticckks_lab=split(sprintf('%.3g ',ticckks),' '); ticckks_lab=ticckks_lab(~cellfun('isempty',ticckks_lab));
  set(h, 'Ticks', ticckks); set(h,'ticklabels',ticckks_lab);
  ylabel(h,'time $t$ [s]','interpreter','latex');
  
  % Customised axis cosmetics.
  minmax=[min(sig_x), max(sig_x), ...
          min(sig_y), max(sig_y)];
  tenPow_x=10.^floor(log10(abs(minmax))); % find power of ten
%   tenPow_x=max(tenPow_x)+tenPow_x*0; % save only highest power of ten on each axis
%   tenPow_x(1:2)=max(tenPow_x(1:2))+tenPow_x(1:2)*0; tenPow_x(3:4)=max(tenPow_x(3:4))+tenPow_x(3:4)*0; % save only highest power of ten on each axis
  minmaxrenormed = minmax./tenPow_x; % remove power of ten
  minmaxrenormed_roundedup = ceil(abs(minmaxrenormed)); % round up
%   minmaxrenormed_roundedup-minmaxrenormed
  selnotenoughmargin=(minmaxrenormed_roundedup-abs(minmaxrenormed)<0.5); % find where change was not enough
  minmaxrenormed_roundedup(selnotenoughmargin) = minmaxrenormed_roundedup(selnotenoughmargin)+1; % up again if change was not enough
  minmax_nextint=sign(minmax).*minmaxrenormed_roundedup.*tenPow_x;
  axis(minmax_nextint);
  if(overdrawcolor=='w')
    set(gca,'Color','k');
    set(gca,'GridColor','white');
  end
  set(gca, 'TickLabelInterpreter', 'latex');
  set(gca,'TickDir','both');
  grid on;
  box on;
  xlabel([lab_x]);
  ylabel([lab_y]);
  title(titlefig);
  
  % Get and plot global rotation trend.
  [th, ~] = cart2pol(sig_x,sig_y);
  vecn=vecnorm([sig_x;sig_y]);
  weight=vecn/max(vecn); % weight by the signal's norm
%   mean(gradient(th,time))
  WGrad=weight.*gradient(th,time); % weighted gradient
%   figure();
%   subplot(3,1,1);plot(time, vecn);
%   subplot(3,1,2);plot(time, th);
%   subplot(3,1,3);plot(time, WGrad);
  gradtheta = mean(WGrad); % derivative of angle wrt time
  disp(['[',mfilename,'] Weighted rotation derivative found to be ',num2str(gradtheta),' rad/s.']);
  ax=gca();
  yl = ax.YLim;
  xl = ax.XLim;
  d = 0.9*min(diff(xl),diff(yl));
  hold on;
  rectangle('Position',[mean(xl)-d/2, mean(yl)-d/2, d, d],'Curvature',[1,1],'edgecolor',overdrawcolor,'linewidth',lw, 'linestyle', '--');
  l = 0.075*diff(yl);
  quiver(mean(xl)+d/2, mean(yl)-sign(gradtheta)*l/2, 0, sign(gradtheta)*l, 0,'color',overdrawcolor,'linewidth',2*lw,'maxheadsize',50);
  quiver(mean(xl)-d/2, mean(yl)+sign(gradtheta)*l/2, 0, -sign(gradtheta)*l, 0,'color',overdrawcolor,'linewidth',2*lw,'maxheadsize',50);

  % Add eigvector angles.
  ax=gca();
  yl = ax.YLim;
  xl = ax.XLim;
  shift=0.02*min(diff(xl),diff(yl));
  fs=ax.FontSize;
  txteigvec=['EigAng (',sprintf(frmt_eigang_plot,th_eig(ilmax)*180/pi),', ',sprintf(frmt_eigang_plot,th_eig(ilmin)*180/pi),')$^\circ$'];
  txteigval=['EigVal (',sprintf(frmt_eigval_plot,eigval(ilmax)),', ',sprintf(frmt_eigval_plot,eigval(ilmin)),')'];
  text(min(xl)+shift, min(yl)+shift, [txteigvec, ', ', txteigval], 'color', overdrawcolor, 'fontsize', 0.5*fs, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
  
  % Eventually plot propagation direction arrow.
  if(not(strcmp(wavepropdirection,'none')))
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
    quiver(P(1,1),P(2,1),dp(1),dp(2),0,'color',overdrawcolor,'linewidth',2*lw,'maxheadsize',0.5);
    fs=ax.FontSize;
    text(mean(xl), P(2,1) + per*0.6, 'wave propagation direction', 'color', overdrawcolor, 'fontsize', 0.9*fs, 'HorizontalAlignment', 'center');
  end
end