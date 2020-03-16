function [fh,X,Y,V,Xi, Yi, Vi] = plotDumpsWrapper(OFD, IT, verbose, nx, ny)
  removeMeanPressure = 1;
  
  tags = {'rho', 'pre', 'velx', 'velz'};
  units = {'kg/m$^3$', 'Pa', 'm/s', 'm/s'};
  
  Ntags = numel(tags);
  
  [X,Y,V] = readDumpsUnique(OFD, IT, verbose);
  [Xi, Yi, Vi] = interpDumps(X, Y, V, nx, ny);
  
  fh = figure('units','normalized','outerposition',[0,0,1,0.5],'name',[shorten_string(OFD,90,6), ' @ITERATION ',num2str(IT)]);
  tightAxes = tight_subplot(1, Ntags, [0.05, 0.05], [.15, .11], [0.06, 0.04]); % gaph gapw marghlow marghupp margwlef margwrig
  
  prefixes = {};
  factores = {};
  for t=1:Ntags
    axes(tightAxes(t));
    curVi = Vi.(tags{t});
    if(not(isempty(curVi)))
      if(removeMeanPressure && strcmp(tags{t}, 'pre'))
        meanPre = mean(curVi,'all');
        curVi = curVi - meanPre;
        tags{t} = [tags{t}, '$-',scientific_latex_notation(meanPre, 1),'$'];
      end
      [prefixes{t}, factores{t}] = prefix_factor_values({curVi});
      pcolor(Xi, Yi, curVi*factores{t});
    end
    title([tags{t}, ' [', prefixes{t},units{t},']']);
    shading interp;
    colormaps_fromPython('seismic', 1);
%   caxis([-1,1]*max(abs(Vi(:))));
    h{t} = colorbar;
%     ytl = split(sprintf('%.4f|',h{t}.Ticks),'|'); ytl(end)=[]; set(h{t},'ticklabels',ytl);
    xlabel(['$x$ [m]']);
    if(t==1)
      ylabel(['$z$ [m]']);
    else
      yticklabels({});
    end
    daspect([1,1,1]);
  end
  linkaxes(tightAxes,'x');
  linkaxes(tightAxes,'y');
  linkprop(tightAxes,'xtick'); xticks('auto');
  linkprop(tightAxes,'ytick'); yticks('auto');
end