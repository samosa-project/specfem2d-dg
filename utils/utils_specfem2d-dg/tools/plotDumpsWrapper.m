function [fh,X,Y,V,Xi, Yi, Vi] = plotDumpsWrapper(OFD, IT, verbose, forceDGMesh, nx, ny)
  removeMeanPressure = 1;
  do_factors = 0;
  
  tags = {'rho', 'pre', 'velx', 'velz'};
  units = {'kg/m$^3$', 'Pa', 'm/s', 'm/s'};
  
  Ntags = numel(tags);
  nrows = 2; ncols = ceil(Ntags/nrows);
  
  [X,Y,V] = readDumpsUnique(OFD, IT, verbose);
  [Xi, Yi, Vi] = interpDumps(X, Y, V, nx, ny, forceDGMesh);
  
  fh = figure('units','normalized','outerposition',[0,0,1,1],'name',[shorten_string(OFD,90,6), ' @ITERATION ',num2str(IT)]);
  tightAxes = tight_subplot(nrows, ncols, [0.15, 0.05], [.12, .06], [0.04, 0.02]); % gaph gapw marghlow marghupp margwlef margwrig
  
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
      if(do_factors)
        [prefixes{t}, factores{t}] = prefix_factor_values({curVi});
      else
        prefixes{t} = '';
        factores{t} = 1;
      end
      pcolor(Xi, Yi, curVi*factores{t});
      title([tags{t}, ' [', prefixes{t},units{t},']']);
      h{t} = colorbar;
      caxis([-1,1]*max(abs(curVi(:)*factores{t})));
    end
%     shading interp;
    shading flat;
    colormaps_fromPython('seismic', 1);
%     ytl = split(sprintf('%.4f|',h{t}.Ticks),'|'); ytl(end)=[]; set(h{t},'ticklabels',ytl);
    xlabel(['$x$ [m]']);
%     if(t==1)
      ylabel(['$z$ [m]']);
%     else
%       yticklabels({});
%     end
%     daspect([1,1,1]);
  end
  linkaxes(tightAxes,'x');
  linkaxes(tightAxes,'y');
  linkprop(tightAxes,'xtick'); xticks('auto');
  linkprop(tightAxes,'ytick'); yticks('auto');
end