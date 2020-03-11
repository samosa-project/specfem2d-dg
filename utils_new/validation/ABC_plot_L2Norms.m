function [fh] = ABC_plot_L2Norms(L2SqrdErr_FAF_grp, L2SqrdErr_BUF_grp, DSPLNM_strct, colours_runs, allCases, outputfigpath)
  
  FAFCol = colours_runs{2};
  BUFCol = colours_runs{3};
  
  fh = figure('units','normalized','outerposition',[0,0,0.4,1]);
  tightAxes = tight_subplot(1, 1, [1,1]*0, [.1, .11], [0.2, 0.02]);
  axes(tightAxes);
  
  MS = 70;
  LW = 4.5;
  lip = 0.2;
  shift_l = 0.03;
  shift_r = -shift_l;
  
  for i=1:3
    plot(i+shift_r, mean(L2SqrdErr_BUF_grp{i}), '.', 'color', BUFCol, 'markersize', MS, 'handlevisibility', 'off'); hold on;
  %   plot([1,1]*i, [min(sqrd_L2_Norm_of_Err_BUF), max(sqrd_L2_Norm_of_Err_BUF)], '-', 'color', BUFCol); hold on;
    linehandles = drawBrackets(gca, i+shift_r, [min(L2SqrdErr_BUF_grp{i}), max(L2SqrdErr_BUF_grp{i})], 2, lip); set(linehandles, 'linewidth', LW, 'color', BUFCol);
    
    plot(i+shift_l, mean(L2SqrdErr_FAF_grp{i}), '.', 'color', FAFCol, 'markersize', MS, 'handlevisibility', 'off'); hold on;
  %   plot([1,1]*i, [min(sqrd_L2_Norm_of_Err_FAF), max(sqrd_L2_Norm_of_Err_FAF)], '-', 'color', FAFCol); hold on;
    linehandles = drawBrackets(gca, i+shift_l, [min(L2SqrdErr_FAF_grp{i}), max(L2SqrdErr_FAF_grp{i})], 2, lip); set(linehandles, 'linewidth', LW, 'color', FAFCol);
  end
  % artificial plots for nicer legend
  plot([1,1]*1e1, [1,1]*1e1, '-', 'displayname', DSPLNM_strct.FAF, 'marker', '.', 'linewidth', LW, 'markersize', MS*0.75, 'color', BUFCol, 'handlevisibility', 'on'); hold on;
  plot([1,1]*1e1, [1,1]*1e1, '-', 'displayname', DSPLNM_strct.BUF, 'marker', '.', 'linewidth', LW, 'markersize', MS*0.75, 'color', FAFCol, 'handlevisibility', 'on'); hold on;
  
  set(tightAxes, 'yscale', 'log');
  
  xlim([0,4]); xticks(1:3); xticklabels(allCases);
  xlabel(['test case']);
  
  ylim([1e-11,1e-7])
  ylabel(['squared L$^2$ norm of error [Pa$\cdot$s]']);
  
  title({['L$^2$ Norms of Time Series'],['Errors for Each Test Case']});
  
  legend('location', 'northeast');
  
  customSaveFig(fh,outputfigpath,{'fig', 'png', 'eps', 'tex'});
end

