L2SqrdErr_FAF_grp = {};
L2SqrdErr_BUF_grp = {};

% Load everything.
for i_case = 1:numel(allCases)
  curCase = allCases{i_case};

  % build OUTPUT_FILES directories' paths for the current case
  curCase_OF_LAR = OFDIRS_strct.(curCase).LAR;
  curCase_OF_FAF = OFDIRS_strct.(curCase).FAF;
  curCase_OF_BUF = OFDIRS_strct.(curCase).BUF;
  % retrieve parameters for the current case
  [curCase_INFOS] = ABC_load_parameters(curCase_OF_LAR, curCase_OF_FAF, curCase_OF_BUF, curCase, BUF_params);
  % print out parameters to text file
  ABC_print_params(curCase_INFOS, outputDir);
  % prepare some parameters
  FIGTITLE = curCase;
  outputFigDir_wprefix=[outputDir,curCase,'_'];
  % set distance and factor by which multiply the absolute error between methods for the plot
  switch curCase
    case 'PW'
      distType_x0z1d2=1;
      errorFactor = 1e5;
    case 'PS'
      distType_x0z1d2=0;
      errorFactor = 1e2;
    case 'WPS'
      distType_x0z1d2=2;
      errorFactor = 1e2;
    otherwise
      error('sdfgjsmgkmsfkgl');
  end

  % check stations agree between runs for each case
  [~,diffSTATBufferLarge] = system(['diff ',curCase_OF_BUF,'STATIONS ',curCase_OF_LAR,'STATIONS']);
  [~,diffSTATFFLarge] = system(['diff ',curCase_OF_FAF,'STATIONS ',curCase_OF_LAR,'STATIONS']);
  if(numel(diffSTATBufferLarge)>0)
    error('STATIONS were not the same between runs buffer and large');
  end
  if(numel(diffSTATFFLarge)>0)
    error('STATIONS were not the same between runs farfield and large');
  end

  % load stations
  [x_stat, z_stat, y_stat, stat_file] = loadStations(curCase_OF_BUF);
  nbstats = numel(x_stat);

  % load actual time series, and compute and store errors
  [time, Zamp, NMs, COLs, LSs, ord, L2SqrdErr_FAF_grp{i_case}, L2SqrdErr_BUF_grp{i_case}] = ABC_load_TS_and_compute_error({curCase_OF_LAR, curCase_OF_FAF, curCase_OF_BUF}, nbstats, subsample, subsample_dt, errorFactor, colours_runs, DSPLNM_strct);
  
  if(strcmp(caseToPlot, curCase))
    x_stat_toplot = x_stat;
    z_stat_toplot = z_stat;
    time_toplot = time;
    Zamp_toplot = Zamp;
    NMs_toplot = NMs;
    COLs_toplot = COLs;
    LSs_toplot = LSs;
    errorFactor_toplot = errorFactor;
    distType_x0z1d2_toplot = distType_x0z1d2;
  end
end

% figure

% plot the one we care to plot
% Time series plot
nRepetitions = size(Zamp_toplot, 1)/nbstats;
switch(distType_x0z1d2_toplot)
  case 0
    distttt = x_stat_toplot; distSymbol='$x$';
  case 1
    distttt = z_stat_toplot; distSymbol='$z$';
  case 2
    % ASSUMING SOURCE IS AT (0, 0)
    % FLEMME TO DO IT WITH XSOURCEZSOURCE
    distttt = (x_stat_toplot.^2+z_stat_toplot.^2).^0.5; distSymbol='$d$';
  otherwise
    error('distance switch not implemented');
end
distanceee = repmat(distttt,[nRepetitions,1]);

% figure
fh_TS = figure('units','normalized','outerposition',[0,0,1,0.8]);
tightAxes = tight_subplot(1, 2, [0.,0.06], [0.13,0.07], [0.08, 0.01]);

ABC_plot_L2Norms(tightAxes(1), L2SqrdErr_FAF_grp, L2SqrdErr_BUF_grp, DSPLNM_strct, colours_runs, allCases);
title(['All $L^2$ Errors']);

axes(tightAxes(2));
distoverptp = 150;
hleg = [];
for i=1:size(Zamp_toplot, 1)
  h = plot(time_toplot, distanceee(i) + distoverptp*Zamp_toplot(i,:), 'displayname', NMs_toplot{i}, 'color', COLs_toplot{i}, 'linestyle', LSs_toplot{i}); hold on
  if(ismember(i,[1:numel(distttt):size(Zamp_toplot, 1)]))
    hleg = [hleg, h];
  end
end
xlabel(['time [s]']);
ylabel([distSymbol,' [m]']);
yticks(sort(distttt));
xlim([0, max(time)]);
ll = legend(hleg, 'units','normalized');
ll.Position([1,2]) = [0.528,0.69];
title(['Time-Distance Plot for the ',  caseToPlot, ' Case']);

hshift = 0.2;
set(tightAxes(1), 'position', tightAxes(1).Position - [0,0,hshift, 0]);
set(tightAxes(2), 'position', tightAxes(2).Position - [hshift,0,-hshift, 0]);

ll = add_labels_subplots(gcf, 0.9, 6);
