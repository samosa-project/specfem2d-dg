% Author:        Léo Martire.
% Description:   Validates 3 simulations between each other to validate
%                absorbing boundary conditions.
% Notes:         TODO.
%
% Usage:
%   1) Run the simulations (see './EXAMPLES/LNSABC_info').
%   2) Set the "Parameters section as you want.
%   3) Run this script.
% with:
%   TODO.
% yields:
%   TODO.

clear all;
% close all;
clc;
format compact;
addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new')); % plot_time_v_dist, plot_total_energy, truncToShortest, readAndSubsampleSynth
SPCFMEXloc = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputDir = [SPCFMEXloc, 'LNSABC_info/']; % don't forget ending '/'
allCases = {'PW', 'PS', 'WPS'};

% prepare colours once and for all
defaultColOrd = get(0,'defaultAxesColorOrder');
%   colour_LAR = 'k';
%   colour_FAF = 'c';
%   colour_BUF = 'm';
colour_LAR = defaultColOrd(1,:);
colour_FAF = defaultColOrd(2,:);
colour_BUF = defaultColOrd(3,:);
colours_runs = {colour_LAR, colour_FAF, colour_BUF};

subsample = 1;
subsample_dt = 1e-4;
LarRNam='large'; FafRNam='far-field'; BufRNam='buffers';

compareEnergies = 0; % compare energies for all cases at once
% compareStations = 0;
compareStations = 1; caseToPlot = ''; % plot simply errors
% compareStations = 1; caseToPlot = 'PW'; % compare time series for one case
% compareStations = 1; caseToPlot = 'PS'; % compare time series for one case
% compareStations = 1; caseToPlot = 'WPS'; % compare time series for one case

% Relative locations of OUTPUT_FILES directories. This should be the same across all test cases.
LAR_ext = '_large/OUTPUT_FILES_baseline/';
FAF_ext = '_FF/OUTPUT_FILES_FF/';
% BUF_ext   = ['_buffers/OUTPUT_FILES_eps0p25/']; BUF_varepsilon='$\varepsilon=0.25$';
% BUF_ext   = ['_buffers/OUTPUT_FILES_p3p25_q6_e1em4/']; BUF_varepsilon='$\varepsilon=10^{-4}$';
% BUF_ext   = ['_buffers/OUTPUT_FILES_p3p25_q6_e0p250/']; BUF_varepsilon='$\varepsilon=0.25$';
% BUF_ext   = ['_buffers/OUTPUT_FILES_p3p25_q6_e0p160_woE/']; BUF_varepsilon='$\varepsilon=0.16$';
% BUF_ext   = ['_buffers/OUTPUT_FILES_p3p25_q6_e0p001/']; BUF_varepsilon='$\varepsilon=10^{-3}$';
% BUF_ext   = ['_buffers/OUTPUT_FILES_p3p25_q6_e0p500/']; BUF_varepsilon='$\varepsilon=0.5$';
BUF_ext   = ['_buffers/OUTPUT_FILES_p3p25q6e0p001/']; BUF_params=[1e-3,3.25,6];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build database of OUTPUT_FILES locations.
OFDIRS_strct = {}; DSPLNM_strct = {};
OFDIRS_strct.PW = {}; bsnm = 'LNSABC_PW';
OFDIRS_strct.PW.LAR = [SPCFMEXloc,bsnm,LAR_ext]; DSPLNM_strct.LAR = LarRNam; LS_strct.LAR = '-';
OFDIRS_strct.PW.FAF = [SPCFMEXloc,bsnm,FAF_ext]; DSPLNM_strct.FAF = FafRNam; LS_strct.FAF = '--';
OFDIRS_strct.PW.BUF = [SPCFMEXloc,bsnm,BUF_ext]; DSPLNM_strct.BUF = BufRNam; LS_strct.BUF = ':';
OFDIRS_strct.PS = {}; bsnm = 'LNSABC_PS';
OFDIRS_strct.PS.LAR = [SPCFMEXloc,bsnm,LAR_ext];
OFDIRS_strct.PS.FAF = [SPCFMEXloc,bsnm,FAF_ext];
OFDIRS_strct.PS.BUF = [SPCFMEXloc,bsnm,BUF_ext];
OFDIRS_strct.WPS = {}; bsnm = 'LNSABC_WPS';
OFDIRS_strct.WPS.LAR = [SPCFMEXloc,bsnm,LAR_ext];
OFDIRS_strct.WPS.FAF = [SPCFMEXloc,bsnm,FAF_ext];
OFDIRS_strct.WPS.BUF = [SPCFMEXloc,bsnm,BUF_ext];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot energy
if(compareEnergies)
%   outputFigPath=[outputFigDir_wprefix,'KE'];
%   fKE = plot_total_energy({largeRunOF,farFieldRunOF,bufferRunOF},0,1,tit_all);
%   customSaveFig(outputFigPath, {'fig', 'jpg', 'eps'});
%   outputFigPath=[outputFigDir_wprefix,'PE'];
%   fPE = plot_total_energy({largeRunOF,farFieldRunOF,bufferRunOF},1,0,tit_all);
%   customSaveFig(outputFigPath, {'fig', 'jpg', 'eps'});
  
%   EEE_outputFigPath=[outputFigDir_wprefix,'E'];
%   fE = plot_total_energy({largeRunOF,farFieldRunOF,bufferRunOF},4,1,FIGTITLE,[],{'-','--',':'}, displaynames_large_ff_buf);
%   customSaveFig(fE, EEE_outputFigPath, {'fig', 'jpg', 'eps'});
%   disp(["figure(fE.Number); customSaveFig(EEE_outputFigPath, {'fig', 'jpg', 'eps'});"]);
  
  energy_outpath = [outputDir,'energy_summary'];
  [fE, axxx] = ABC_plot_NRJs(OFDIRS_strct,DSPLNM_strct,LS_strct,colours_runs,energy_outpath);
  %TODO: export descriptions (INFO_all) to a text file
end

% compare time series
if(compareStations)
  L2SqrdErr_FAF_grp = {};
  L2SqrdErr_BUF_grp = {};
  
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
    
    % plot the one we care to plot
    if(strcmp(caseToPlot, curCase))
      % Time series plot
      TS_outputFigPath=[outputFigDir_wprefix,'time_series'];
      nRepetitions = size(Zamp, 1)/nbstats;
      switch(distType_x0z1d2)
        case 0
          distttt = x_stat; distSymbol='$x$';
        case 1
          distttt = z_stat; distSymbol='$z$';
        case 2
          % ASSUMING SOURCE IS AT (0, 0)
          % FLEMME TO DO IT WITH XSOURCEZSOURCE
          distttt = (x_stat.^2+z_stat.^2).^0.5; distSymbol='$d$';
        otherwise
          error('distance switch not implemented');
      end
      fTvD = plot_time_v_dist(repmat(time,[nbstats*nRepetitions,1]), Zamp, repmat(distttt,[nRepetitions,1]), 0, FIGTITLE, distSymbol, 1, NMs, COLs, LSs);
      
      % cosmetics
      figure(fTvD);
      cax=gca();
      set(cax.Children(1:nRepetitions*(nbstats-1)),'handlevisibility','off');
      legend('location','northwest');
      prettyAxes(fTvD);
      customSaveFig(fTvD, TS_outputFigPath, {'fig', 'jpg', 'eps'});
      disp(["customSaveFig(fTvD, TS_outputFigPath, {'fig', 'jpg', 'eps'});"]);
    end
  end
  
  % plot L2 errors
  FAFCol = COLs{ord.Ff*nbstats+1};
  BUFCol = COLs{ord.Bu*nbstats+1};
  l2err_outpath = [outputDir, 'TS_L2_summary'];
  [fh] = ABC_plot_L2Norms(L2SqrdErr_FAF_grp, L2SqrdErr_BUF_grp, DSPLNM_strct, colours_runs, allCases, l2err_outpath);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One-liners for debug.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time series: plot only some
if(0)
%   or=orEFf;
  or=orEBu;
%   or=orLa;
%   or=orFf;
%   or=orBu;
  plot_time_v_dist(repmat(time,[nbstats,1]), Zamp(or*nbstats+(1:nbstats),:), repmat(z_stat,[1,1]), 0, '', '$z$', 1);
end
% time series: 1-by-1 plot
if(0)
  z=40;
  dataL = readAndSubsampleSynth(largeRunOF, find(z_stat==z), 'BXZ', 'semv', subsample, subsample_dt, 1);
  dataF = readAndSubsampleSynth(farFieldRunOF, find(z_stat==z), 'BXZ', 'semv', subsample, subsample_dt, 1);
  dataB = readAndSubsampleSynth(bufferRunOF, find(z_stat==z), 'BXZ', 'semv', subsample, subsample_dt, 1);
  figure();
  subplot(1,2,1); plot(dataL(:,2)); hold on; plot(dataF(:,2)); hold on; plot(dataB(:,2));
  subplot(1,2,2); plot(errorFactor*abs(dataL(:,2)-dataF(:,2))); hold on; plot(errorFactor*abs(dataL(:,2)-dataB(:,2)));
end

function autoSymLog(axxx, y)
  ydataamplitudes = cellfun(@max,get(ch(1).Children,'ydata'))-cellfun(@min,get(ch(1).Children,'ydata'));
  maxydataamplitudes = max(ydataamplitudes);
  symlog(ch(1),'y',ceil(maxydataamplitudes)-4);
end