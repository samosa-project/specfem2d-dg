function [] = plot_bg_model(DATAFILE)
  bgm = load_bg_model(DATAFILE);
  [order, tag, tex, unit] = order_bg_model();
  nb_qty = size(order, 1);
  
  dbgmx = diff(bgm.xx);
  dbgmz = diff(bgm.zz);
  if(not(any(diff(diff(find(dbgmx>0))))))
    % bgm.xx is a series of repeating xx values
    NZ = find(dbgmx>0,1,'first');
    NX = numel(bgm.xx)/NZ;
    for iqty = 1:nb_qty
      bgm.(order(iqty,:)) = reshape(bgm.(order(iqty,:)), NZ, NX);
    end
    isMatrix = 1;
  else
    isMatrix = 0;
  end
  
  toKm = 1;
  if(toKm)
    bgm.xx = bgm.xx/1e3;
    bgm.zz = bgm.zz/1e3;
    unit{1}=['k',unit{1}];
    unit{2}=['k',unit{2}];
  end
  
  thefig = figure('units','normalized','outerposition',[0,0,1,1]);
  nlin=2; ncol = 4;
  tightAxes = tight_subplot(nlin, ncol, [0.08,0.03], [0.12,0.12], [0.07,0.03]);
  XLAB = ['$',tex{1},'$ [',unit{1},']'];
  YLAB = ['$',tex{2},'$ [',unit{2},']'];

  for iqty = 3:nb_qty
    axes(tightAxes(iqty-2));
    toPlot = bgm.(order(iqty,:));
    if(strcmp(tag{iqty},'rho0[kg.m^{-3}]') || strcmp(tag{iqty},'p0[Pa]'))
      toPlot = log(toPlot);
      tex{iqty} = ['\log_{10}\left(',tex{iqty},'\right)'];
      doingLog = 1;
    else
      doingLog = 0;
    end
    if(isMatrix)
      pcolor(bgm.xx, bgm.zz, toPlot);
%       if(all(size(unique(toPlot))==[NZ, 1])==1)
%         % only varying vertically
%         plot(toPlot, bgm.zz);
%         onlyVaryingVertically = 1;
%       else
%         pcolor(bgm.xx, bgm.zz, toPlot);
%         onlyVaryingVertically = 0;
%       end
    else
      error('not implemented');
    end
    if(mod(iqty-2, ncol)==1)
      ylabel(YLAB);
    else
      yticklabels({});
    end
    if(iqty-2>(nlin-1)*ncol)
      xlabel(XLAB);
    else
      xticklabels({});
    end
    if(numel(unique(toPlot(:)))==1)
      % constant value
      cst = unique(toPlot(:));
      colormaps_fromPython('seismic', 1);
      TIT = ['$',tex{iqty},'=',num2str(cst),'$'];
    else
      hcb = colorbar;
      TIT = ['$',tex{iqty},'$'];
      if(sign(min(toPlot(:)))*sign(max(toPlot(:)))==-1 & not(doingLog))
        % opposed sign min and max, use seismic cb
        colormaps_fromPython('seismic', 1);
        caxis([-1,1]*max(abs(toPlot(:))));
      else
        % keep inferno
      end
    end
    title(TIT);
  end
  
  customSaveFig(thefig,regexprep(DATAFILE,'.bin',''),{'png','eps'}, 9999);
end

