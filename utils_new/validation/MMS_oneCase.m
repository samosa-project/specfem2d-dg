function [errQtity, err_l2, err_rel, globalSave, GS_ID_NX, GS_ID_DT, GS_ID_IT, GS_ID_EPS, GS_ID_RELERR] = MMS_oneCase(prefix, testCase, OFDs, itToPlot_forEachOFD, plotFields_do, err_l2, err_rel, verbose, savefigpath)

  % N for interpolation grid.
  Nexact = 2000; % Choose Nexact >> the N's we want to test. 10000 chugs HARD, 1000 runs ok. Up to N=100 points on the simulation's mesh, Nexact=1000 and Nexact=5000 yield the exact same L2 error.
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Treatment.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(not(numel(OFDs)==numel(itToPlot_forEachOFD)))
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
  GS_ID_RELERR = 5;
  for ofdi = 1:N_OFD
    OFD = OFDs{ofdi};
    IDz = itToPlot_forEachOFD{ofdi};
    if(not(OFD(end)=='/')); OFD = [OFD,'/']; end; % Safeguard.
    if(verbose)
      disp(['[',mfilename,'] Looking at ''',OFD,'''.']);
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
    plotFields_synthname = [LNSorFNS,' Synthetic'];
%     plotFields_synthname = ['Synthetic'];
    KAPPA = readExampleFiles_extractParam(parfile,'thermal_conductivity','float');
    MU = readExampleFiles_extractParam(parfile,'dynamic_viscosity','float');
    imagetype_wavefield_dumps = readExampleFiles_extractParam(parfile, 'imagetype_wavefield_dumps', 'int');
    [~] = MMS_check(testCase, MU, KAPPA);
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

      % Build analytic solution (cf. LNS_manufactured_solutions.mw) over grid.
      [drho_th, dvx_th, dvz_th, dE_th, dp_th] = MMS_analytic(testCase, Xe, Ye, RHO0, V0_x, E0, P0, GAM);
      
      % Grab synthetic value.
      [Xquery, Yquery, syntheticValue] = MMS_grabSyntheticValue(OFD);
      % Build analytic solution at query point.
      [~, ~, ~, ~, dp_th_query] = MMS_analytic(testCase, Xquery, Yquery, RHO0, V0_x, E0, P0, GAM);

      switch(imagetype_wavefield_dumps)
        case 4
          zexp = zexp.pre;
          zth = dp_th;
          errQtity = 'p';
          QtityUnit = 'Pa';
          if(verbose)
            disp(['[',mfilename,'] Dumps are pressure dumps, computing and plotting error w.r.t. pressure.']);
          end
        case 2
          zexp = zexp.vel;
          zth = dvx_th;
          errQtity = 'v_x';
          QtityUnit = 'm/s';
          if(verbose)
            disp(['[',mfilename,'] Dumps are horizontal velocity dumps, computing and plotting error w.r.t. horizontal velocity.']);
          end
        otherwise
          error(['[] imagetype_wavefield_dumps not implemented']);
      end
      % update error names
      err_l2.unit = [QtityUnit,'$\cdot$m'];
      plotFields_YLABCB = ['${\Delta}',errQtity,'$ [',QtityUnit,']'];
      err_l2.prefix = ['$\',err_l2.symbol, '_{',errQtity,'}'];
      err_rel.prefix = ['$\',err_rel.symbol, '_{',errQtity,'}'];
      
      % Compute error.
      zerr = zexp-zth;
      err_l2_value = (trapz(x_exact,trapz(y_exact,zerr.^2)))^0.5;
      err_rel_value = abs((syntheticValue-dp_th_query)/dp_th_query);
      
      if(verbose)
        disp(['[',mfilename,'] Integral of zth = ',num2str(trapz(x_exact,trapz(y_exact,zth)))]);
        disp(['[',mfilename,'] L_{inf} norm of error = ',num2str(max(max(zerr)))]);
      end

      globalSave(ofdi,IT_id,GS_ID_IT) = IT; % save iteration in last dimension
      globalSave(ofdi,IT_id,GS_ID_EPS) = err_l2_value; % save error in last dimension
      globalSave(ofdi,IT_id,GS_ID_RELERR) = err_rel_value;

      % Plot error on field.
      if(plotFields_do)
        plotFields_saveFigName_prefix = ['comparedDump__'];
        plotFields_extsToSave={'png', 'eps'};
        spl=split(OFD,'/');
        plotFields_saveFigName_base = [plotFields_saveFigName_prefix, regexprep(spl{end-2},[prefix,'_'],''),'__',regexprep(spl{end-1},'OUTPUT_FILES_',''),'__'];% if sub OUTPUT_FILES folders per case
        savefigname = [plotFields_saveFigName_base, 'IT',sprintf('%012.0f',IT)];
        savefigfullpath = [savefigpath, regexprep(savefigname,'\.','')]; % regexprep because latex crashed when filenames have dots in it
        feg = MMS_plotFields(err_l2, err_l2_value, Xe, Ye, zth, zexp, zerr, IT, DT, NX, plotFields_synthname,plotFields_YLABCB);
        customSaveFig(feg, savefigfullpath, plotFields_extsToSave, 9999);
      end
      
    end % Endfor on IT.
  end % Endfor on OFD.
  %%%%%%%%%%%%%%%%

  squeeze(globalSave)

  %%%%%%%%%%%%%%%%
end