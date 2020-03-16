function [fh] = plotDumpsWrapper(OFD, IT, verbose, nx, ny)
  tags = {'rho', 'vel', 'pre'};
  
  Ntags = numel(tags);
  
  [X,Y,V] = readDumpsUnique(OFD, IT, verbose);
  [Xi, Yi, Vi] = interpDumps(X, Y, V, nx, ny);
  
  fh = figure('units','normalized','outerposition',[0,0,1,1]);
  tightAxes = tight_subplot(1, Ntags, [0, 0.05], [.1, .08], [0.06, 0.04]); % gaph gapw marghlow marghupp margwlef margwrig
  
  for t=1:Ntags
    axes(tightAxes(t));
    if(not(isempty(V.(tags{t}))))
      pcolor(Xi, Yi, Vi.(tags{t}));
    end
    title(tags{t});
    shading interp;
    colormaps_fromPython('seismic', 1);
%   caxis([-1,1]*max(abs(Vi(:))));
    h = colorbar;
    xlabel(['$x$ [m]']);
    if(t==1)
      ylabel(['$z$ [m]']);
    else
      yticklabels({});
    end
  end
  linkaxes(tightAxes,'x');
  linkaxes(tightAxes,'y');
  linkprop(tightAxes,'xtick'); xticks('auto');
  linkprop(tightAxes,'ytick'); yticks('auto');
end