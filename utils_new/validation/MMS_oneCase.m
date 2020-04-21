function MMS_oneCase(testCase, OFDs, iterationsToPlot_forEachOFD, plotFields_do, verbose, savefigpath)

  % N for interpolation grid.
  Nexact = 2000; % Choose Nexact >> the N's we want to test. 10000 chugs HARD, 1000 runs ok. Up to N=100 points on the simulation's mesh, Nexact=1000 and Nexact=5000 yield the exact same L2 error.
  errorPrefix_base = ['$\varepsilon_{'];
  
  % Field plots?
  plotFields_do = 1; % do or do not plot fields themselves
  plotFields_saveFigName_prefix = ['comparedDump__'];
  plotFields_extsToSave={'jpg'};
  plotFields_xlab = '$x/L$';
  plotFields_ylab = '$z/L$';
  plotFields_NColorbarTicks = 9;
  plotFields_addErrorPanel = 1; % Set to 1 for debug, set to 0 for printing nice plots
  plotFields_CMAP_thresh = 0.01; % relative to max
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Treatment.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(not(numel(OFDs)==numel(iterationsToPlot_forEachOFD)))
    error(['[, ERROR] must provide same number of OFD as IDz']);
  else
    N_OFD = numel(OFDs);
  end

  globalSave = []; % prepare saving of errors w.r.t. iterations for each OFD
  % Storage IDs.
  GS_ID_NX = 1;
  GS_ID_DT = 2;
  GS_ID_IT = 3;
  GS_ID_EPS = 4;
  for ofdi = 1:N_OFD
    OFD = OFDs{ofdi};
    IDz = iterationsToPlot_forEachOFD{ofdi};
    if(not(OFD(end)=='/')); OFD = [OFD,'/']; end; % Safeguard.
    if(verbose)
      disp(['[',mfilename,'] Looking at ''',OFD,'''.']);
    end

    if(plotFields_do)
      % Prepare Figure saving w.r.t. OFD.
      spl=split(OFD,'/'); plotFields_saveFigName_base = [plotFields_saveFigName_prefix, regexprep(spl{end-2},[prefix,'_'],''),'__',regexprep(spl{end-1},'OUTPUT_FILES_',''),'__'];% if sub OUTPUT_FILES folders per case
    end

    %%%%%%%%%%%%%%%%
    % Load parameters from simulation.
    %%%%%%%%%%%%%%%%
    parfile = [OFD,'input_parfile'];
    MODEL = readExampleFiles_extractParam(parfile,'MODEL','string');
    USE_ISOTHERMAL_MODEL = readExampleFiles_extractParam(parfile,'USE_ISOTHERMAL_MODEL','bool');
    if(not(strcmp(MODEL,'default') & not(USE_ISOTHERMAL_MODEL)))
      error(['fluid model is wrong, set MODEL=default and USE_ISOTHERMAL_MODEL=.false. ']);
    end
    NX = readExampleFiles_extractParam(parfile,'nx','int');
    RHO0 = readExampleFiles_extractParam(parfile,'surface_density','float');
    sound_velocity = readExampleFiles_extractParam(parfile,'sound_velocity','float');
    V0_x = readExampleFiles_extractParam(parfile,'wind','float');
    GAM = readExampleFiles_extractParam(parfile,'constant_p','float') / readExampleFiles_extractParam(parfile,'constant_v','float');
    P0 = sound_velocity^2*RHO0/GAM; % deduce p0 from isobaric hypothesis
    E0 = P0/(GAM-1) + 0.5*RHO0*V0_x^2; % deduce E0 from isobaric hypothesis
    if(readExampleFiles_extractParam(parfile,'USE_LNS','bool'))
      LNSorFNS = 'LNS';
    else
      LNSorFNS = 'FNS';
    end
    synthname = [LNSorFNS,' Synthetic'];
    KAPPA = readExampleFiles_extractParam(parfile,'thermal_conductivity','float');
    MU = readExampleFiles_extractParam(parfile,'dynamic_viscosity','float');
    imagetype_wavefield_dumps = readExampleFiles_extractParam(parfile, 'imagetype_wavefield_dumps', 'int');
    [VISCOUS] = MMS_check(testCase, MU, KAPPA);
    %%%%%%%%%%%%%%%%

    % Save N.
    globalSave(ofdi,:,GS_ID_NX) = ones(1,numel(IDz))* NX; % save N in last dimension, for all iterations
    % Save DT.
    DT = readExampleFiles_extractParam(parfile,'DT','float');
    globalSave(ofdi,:,GS_ID_DT) = DT;

    %%%%%%%%%%%%%%%%
    % Loop dumps and plot agreement.
    %%%%%%%%%%%%%%%%
    for IT_id = 1:numel(IDz)
      IT = IDz(IT_id);
      % Load dump.
      [X,Y,V] = readDumpsUnique(OFD, IT, 0);
      if(verbose)
        disp(['[',mfilename,'] Read dump at iteration ',num2str(IT),' (t = ',num2str(IT*DT),' s).']);
      end

      % Create independent mesh onto which interpolate.
      [Xe, Ye, zexp] = interpDumps(X, Y, V, Nexact, Nexact);
      x_exact = unique(Xe); y_exact = unique(Ye);

      % Build analytic solution (cf. LNS_manufactured_solutions.mw).
      [drho_th, dvx_th, dvz_th, dE_th, dp_th] = MMS_analytic(testCase, Xe, Ye, RHO0, V0_x, E0, P0, GAM);

      switch(imagetype_wavefield_dumps)
        case 4
          zexp = zexp.pre;
          zth = dp_th;
          qtity = 'p';
          if(verbose)
            disp(['[',mfilename,'] Dumps are pressure dumps, computing and plotting error w.r.t. pressure.']);
          end
        case 2
          zexp = zexp.vel;
          zth = dvx_th;
          qtity = 'v_x';
          if(verbose)
            disp(['[',mfilename,'] Dumps are horizontal velocity dumps, computing and plotting error w.r.t. horizontal velocity.']);
          end
        otherwise
          error(['[] imagetype_wavefield_dumps not implemented']);
      end
      caxxxxiiss_exp = [min(min(zexp)), max(max(zexp))];
      caxxxxiiss_th = max(max(abs(zexp)))*[-1,1];
      errorPrefix = [errorPrefix_base, qtity];
      % Compute error.
      zerr = zexp-zth;
      errName = [errorPrefix,'}\left(',num2str(NX),'\right)'];
      caxxxxiiss_err = [min(min(zerr)), max(max(zerr))];
      L2NormOfError = (trapz(x_exact,trapz(y_exact,zerr.^2)))^0.5;
      if(verbose)
        disp(['[',mfilename,'] Integral of zth = ',num2str(trapz(x_exact,trapz(y_exact,zth)))]);
        disp(['[',mfilename,'] L_{inf} norm of error = ',num2str(max(max(zerr)))]);
      end

      globalSave(ofdi,IT_id,GS_ID_IT) = IT; % save iteration in last dimension
      globalSave(ofdi,IT_id,GS_ID_EPS) = L2NormOfError; % save error in last dimension

      % Plot error on field.
      if(plotFields_do)
        savefigname = [plotFields_saveFigName_base, 'IT',sprintf('%012.0f',IT)];
        savefigfullpath = [savefigpath, regexprep(savefigname,'\.','')]; % regexprep because latex crashed when filenames have dots in it
        MMS_plotFields;
      end
    end
  end
  %%%%%%%%%%%%%%%%

  squeeze(globalSave)

  %%%%%%%%%%%%%%%%
  % Error plots.
  %%%%%%%%%%%%%%%%

  % Error plots.
  plotProgrWRTTime = 0; % plot progression wrt time?
  figErr_saveFigName_base = ['MES_Error__'];
  figErr_extsToSave = {'jpg'};
  plotErr_title = '$L^2$ Error as Function of Number of Elements $N$';
  plotErr_xlab = 'number of elements on one side $N$';
  plotErr_figPos = [0,1e4,1600,600];

  if(size(globalSave,1)>1)
    % If more than one OFD.
    plotErr_ylab = ['$L^2$ error ', errorPrefix,'}(N)$'];
    plotErr_xlim = [min(min(globalSave(:,:,GS_ID_NX))),max(max(globalSave(:,:,GS_ID_NX)))];

    if(plotProgrWRTTime)
      % Plot progression w.r.t. time.
      ids_IT_to_plot = 1:size(globalSave,2);
      colours = winter(size(globalSave,2));
      fig_error_progress = figure('units','pixels','position',plotErr_figPos);
      tight_subplot(1, 1, -1, [0.2, 0.11], [0.1, 0.04]);
      for id_IT_to_plot = ids_IT_to_plot
  %       IT_N=00250dx = globalSave(1,id_IT_to_plot,2);
        DT_IT = squeeze(globalSave(:,id_IT_to_plot,[GS_ID_DT,GS_ID_IT]));
        if(not(numel(unique(prod(DT_IT,2))==1)))
          error(['[, ERROR] t is inconsistent between saved errors (',num2str(prod(DT_IT,2)),')']);
        end
        t = unique(prod(DT_IT,2));
        N_EPS = squeeze(globalSave(:,id_IT_to_plot,[GS_ID_NX,GS_ID_EPS]));
        N = N_EPS(:, 1);
        EPS = N_EPS(:, 2);
        semilogy(N, EPS,'displayname',['@$t=',num2str(t),'$~s'],'color',colours(id_IT_to_plot,:));
        hold on;
      end
      legend();
      xlabel(plotErr_xlab); ylabel(plotErr_ylab); title([plotErr_title,' and Simulation Time $t$']);
      xlim(plotErr_xlim);
      ylim([10^floor(log10(min(min(globalSave(:,:,GS_ID_EPS))))),10^ceil(log10(max(max(globalSave(:,:,GS_ID_EPS)))))]);
      prettyAxes(fig_error_progress);
      figErr_saveFigFullPath = [savefigpath, figErr_saveFigName_base, 'progress__', datestr(now,'YYmmDD_HHMMss')];
      customSaveFig(fig_error_progress, figErr_saveFigFullPath, figErr_extsToSave);
    end

    % Plot last time as standalone.
    id_IT_to_plot = size(globalSave,2);
    N_EPS = squeeze(globalSave(:,id_IT_to_plot,[GS_ID_NX,GS_ID_EPS]));
    N = N_EPS(:, 1);
    EPS = N_EPS(:, 2);
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, colourPlot, markerPlot, markerSize] = MMS_constants(testCase); % Get plot parameters as function of test case.
    fig_error = figure('units','pixels','position',plotErr_figPos);
    tight_subplot(1, 1, -1, [0.17,0.11], [0.1, 0.02]);
    loglog(N, EPS,'displayname',testCase,'marker',markerPlot,'markersize',markerSize,'color',colourPlot,'markerfacecolor',colourPlot);
    legend('location','northeast');
    xlabel(plotErr_xlab); ylabel(plotErr_ylab); title(plotErr_title);
    xlim(plotErr_xlim);
  %   tickLabels_general_reskin_y(fig_error);
    prettyAxes(fig_error);
    figErr_saveFigFullPath = [savefigpath, figErr_saveFigName_base, datestr(now,'YYmmDD_HHMMss')];
    customSaveFig(fig_error, figErr_saveFigFullPath, figErr_extsToSave);
  end
  %%%%%%%%%%%%%%%%
end