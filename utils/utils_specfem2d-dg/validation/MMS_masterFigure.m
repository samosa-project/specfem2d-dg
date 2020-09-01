function MMS_masterFigure(testCase, errQtity, err_l2, err_rel, globalSave, GS_ID_NX, GS_ID_DT, GS_ID_IT, GS_ID_EPS, GS_ID_RELERR, savefigpath, plotProgrWRTTime)
  if(size(globalSave,1)<=1)
    error('nothing to do');
  end
  
  figErr_saveFigName_base = ['MES_Error__'];
  figErr_extsToSave = {'fig', 'png'};
  
  LS = ':'; LW = 2;
  
  TIT_L2 = [err_l2.name, ' ', err_l2.prefix,'\left(N\right)$ [',err_l2.unit,']'];
  TIT_RELERR = [err_rel.name, ' ',err_rel.prefix,'\left(N\right)$ [\%]'];
  XLAB = 'number of elements on one side $N$';
  % errName = [errorPrefix,'}\left(',num2str(NX),'\right)'];
%   YLAB_L2 = ['$L^2$ error ', errorPrefix,'}(N)$'];
%   YLAB_RELERR = ['relative error [\%]'];
  YLAB_L2 = ['']; YLAB_RELERR = [''];
  plotErr_xlim = [min(min(globalSave(:,:,GS_ID_NX))),max(max(globalSave(:,:,GS_ID_NX)))];

  % Plot last time as standalone.
  id_IT_to_plot = size(globalSave,2);
  N_EPS_RELERR = squeeze(globalSave(:,id_IT_to_plot,[GS_ID_NX,GS_ID_EPS,GS_ID_RELERR]));
  N = N_EPS_RELERR(:, 1);
  EPS = N_EPS_RELERR(:, 2);
  RELERR = N_EPS_RELERR(:, 3) * 100;
  [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, colourPlot, markerPlot, markerSize] = MMS_constants(testCase); % Get plot parameters as function of test case.
  
  fig_error = figure('units','pixels','position',[0,1e4,1900,800]);
  tightAxes = tight_subplot(1, 2, [-1,0.07], [0.13,0.08], [0.06, 0.02]);
  
  axes(tightAxes(1));
  loglog(N, EPS,'displayname',testCase,'linestyle', LS, 'linewidth', LW, 'marker',markerPlot,'markersize',markerSize,'color',colourPlot,'markerfacecolor',colourPlot);
  legend('location','northeast');
  xlabel(XLAB); ylabel(YLAB_L2); title(TIT_L2);
  yticks(logspace(log10(min(ylim)), log10(max(ylim)), log10(max(ylim))-log10(min(ylim))+1)); % make sre every tick is printed

  axes(tightAxes(2));
  loglog(N, RELERR,'displayname',testCase,'linestyle', LS, 'linewidth', LW, 'marker',markerPlot,'markersize',markerSize,'color',colourPlot,'markerfacecolor',colourPlot);
  xlabel(XLAB); ylabel(YLAB_RELERR); title(TIT_RELERR);
  YL=ylim; ylim(sort(unique([YL(YL<=100),100]))); % ensure upper ylim is 100
  yticks(logspace(log10(min(ylim)), log10(max(ylim)), log10(max(ylim))-log10(min(ylim))+1)); % make sre every tick is printed
  YTL = get(gca, 'yticklabels'); YTL{end}='100';YTL{end-1}='10';YTL{end-2}='1'; set(gca, 'yticklabels', YTL);
  legend('location','northeast');
  
  linkaxes(tightAxes, 'x');
  xlim(plotErr_xlim);
%   tickLabels_general_reskin_y(fig_error);
  prettyAxes(fig_error);
  figErr_saveFigFullPath = [savefigpath, figErr_saveFigName_base, datestr(now,'YYmmDD_HHMMss')];
  customSaveFig(fig_error, figErr_saveFigFullPath, figErr_extsToSave);
  
  
  if(plotProgrWRTTime)
    % Plot progression w.r.t. time.
    ids_IT_to_plot = 1:size(globalSave,2);
    colours = winter(size(globalSave,2));
    fig_error_progress = figure('units','pixels','position',FIGPOS);
    tight_subplot(1, 1, -1, [0.2, 0.11], [0.1, 0.04]);
    for id_IT_to_plot = ids_IT_to_plot
%       IT_N=00250dx = globalSave(1,id_IT_to_plot,2);
      DT_IT = squeeze(globalSave(:,id_IT_to_plot,[GS_ID_DT,GS_ID_IT]));
      if(not(numel(unique(prod(DT_IT,2))==1)))
        error(['[, ERROR] t is inconsistent between saved errors (',num2str(prod(DT_IT,2)),')']);
      end
      t = unique(prod(DT_IT,2));
      N_EPS_RELERR = squeeze(globalSave(:,id_IT_to_plot,[GS_ID_NX,GS_ID_EPS]));
      N = N_EPS_RELERR(:, 1);
      EPS = N_EPS_RELERR(:, 2);
      semilogy(N, EPS,'displayname',['@$t=',num2str(t),'$~s'],'color',colours(id_IT_to_plot,:));
      hold on;
    end
    legend();
    xlabel(XLAB); ylabel(YLAB_L2); title([TIT_L2,' and Simulation Time $t$']);
    xlim(plotErr_xlim);
    ylim([10^floor(log10(min(min(globalSave(:,:,GS_ID_EPS))))),10^ceil(log10(max(max(globalSave(:,:,GS_ID_EPS)))))]);
    prettyAxes(fig_error_progress);
    figErr_saveFigFullPath = [savefigpath, figErr_saveFigName_base, 'progress__', datestr(now,'YYmmDD_HHMMss')];
    customSaveFig(fig_error_progress, figErr_saveFigFullPath, figErr_extsToSave);
  end
end