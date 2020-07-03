function ABC_plot_L2Norms(aimedAxe, L2SqrdErr_FAF_grp, L2SqrdErr_BUF_grp, DSPLNM_strct, colours_runs, allCases)
  
  DES_BARRES = 0;
  MS = 300;
  markfaf = 's';
  markbuf = 'h';
  
  MSdot = 70;
  LW = 4.5;
  lip = 0.2;
  shift_l = -0.08;
  shift_r = -shift_l;
  
  FAFCol = colours_runs{2};
  BUFCol = colours_runs{3};
  
%   fh = figure('units','normalized','outerposition',[0,0,0.4,1]);
%   tightAxes = tight_subplot(1, 1, [1,1]*0, [.1, .11], [0.2, 0.02]);
  axes(aimedAxe);
  
  for i=1:3
    if(DES_BARRES)
      h2=plot(i+shift_l, mean(L2SqrdErr_FAF_grp{i}), '.', 'color', FAFCol, 'markersize', MSdot, 'handlevisibility', 'off'); hold on;
      linehandles = drawBrackets(gca, i+shift_l, [min(L2SqrdErr_FAF_grp{i}), max(L2SqrdErr_FAF_grp{i})], 2, lip); set(linehandles, 'linewidth', LW, 'color', FAFCol);
      h1=plot(i+shift_r, mean(L2SqrdErr_BUF_grp{i}), '.', 'color', BUFCol, 'markersize', MSdot, 'handlevisibility', 'off'); hold on;
      linehandles = drawBrackets(gca, i+shift_r, [min(L2SqrdErr_BUF_grp{i}), max(L2SqrdErr_BUF_grp{i})], 2, lip); set(linehandles, 'linewidth', LW, 'color', BUFCol);
    else
      n=numel(L2SqrdErr_BUF_grp{i});
      plot(zeros(n,1)+i+shift_l, L2SqrdErr_FAF_grp{i}, '-', 'color', FAFCol, 'displayname', DSPLNM_strct.FAF, 'handlevisibility', 'off'); hold on;
      plot(zeros(n,1)+i+shift_r, L2SqrdErr_BUF_grp{i}, '-', 'color', BUFCol, 'displayname', DSPLNM_strct.BUF, 'handlevisibility', 'off'); hold on;
      h2 = scatter(zeros(n,1)+i+shift_l, L2SqrdErr_FAF_grp{i}, MS*0.8, markfaf, 'filled', 'markerfacecolor', FAFCol, 'markeredgecolor', FAFCol*0.5, 'linewidth', 2, 'displayname', DSPLNM_strct.FAF, 'handlevisibility', 'off'); hold on;
      h1 = scatter(zeros(n,1)+i+shift_r, L2SqrdErr_BUF_grp{i}, MS, markbuf, 'filled', 'markerfacecolor', BUFCol, 'markeredgecolor', BUFCol*0.5, 'linewidth', 2, 'displayname', DSPLNM_strct.BUF, 'handlevisibility', 'off'); hold on;
    end
    if(i==1)
      hleg = [h2,h1];
    end
  end
  % artificial plots for nicer legend
%   plot([1,1]*1e1, [1,1]*1e1, '-', 'displayname', DSPLNM_strct.FAF, 'marker', '.', 'linewidth', LW, 'markersize', MS*0.75, 'color', FAFCol, 'handlevisibility', 'on'); hold on;
%   plot([1,1]*1e1, [1,1]*1e1, '-', 'displayname', DSPLNM_strct.BUF, 'marker', '.', 'linewidth', LW, 'markersize', MS*0.75, 'color', BUFCol, 'handlevisibility', 'on'); hold on;
  
  set(aimedAxe, 'yscale', 'log');
  
  xlim([0,4]);
  xticks(1:3);
  xticklabels(allCases);
  xlabel(['test case']);
  
  ylim([1e-12,0.3e-7]);
  ylabel(['squared $L^2$ norm of error [Pa$\cdot$s]']);
  
%   title({['L$^2$ Norms of Time Series'],['Errors for Each Test Case']});
  
  legend(hleg, 'location', 'southeast');
  
%   customSaveFig(fh,outputfigpath,{'fig', 'png', 'eps', 'tex'});
end

