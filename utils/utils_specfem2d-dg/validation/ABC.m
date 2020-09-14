% Author:        Léo Martire.
% Description:   Validates 3 simulations between each other to validate
%                absorbing boundary conditions.
% Notes:         TODO.
%
% Usage:
%   1) Run the simulations (see './EXAMPLES/LNSABC_info').
%   2) Set the "Parameters" section as you want.
%   3) Run this script.
% with:
%   TODO.
% yields:
%   TODO.

clear all;
% close all;
clc;

[SPCFMEXloc] = setup_overall();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputDir = [SPCFMEXloc, 'validation__lns_abc_setup',filesep]; % don't forget ending '/'
energy_outpath = [outputDir,'energy_summary'];
TS_outpath = [outputDir,'time_series_summary'];
allCases = {'PW', 'PS', 'WPS'};
LarRNam='large'; FafRNam='far-field'; BufRNam='buffers';

compareEnergies = 1; % compare energies for all cases at once
% compareStations = 0;
% compareStations = 1; caseToPlot = ''; % plot simply errors
% compareStations = 1; caseToPlot = 'PW'; % compare time series for one case
% compareStations = 1; caseToPlot = 'PS'; % compare time series for one case
compareStations = 1; caseToPlot = 'WPS'; % compare time series for one case

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

subsample = 1; subsample_dt = 1e-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare colours once and for all
defaultColOrd = get(0,'defaultAxesColorOrder');
%   colour_LAR = 'k';
%   colour_FAF = 'c';
%   colour_BUF = 'm';
colour_LAR = defaultColOrd(1,:);
colour_FAF = defaultColOrd(2,:);
colour_BUF = defaultColOrd(3,:);
colours_runs = {colour_LAR, colour_FAF, colour_BUF};

% Build database of OUTPUT_FILES locations.
OFDIRS_strct = {}; DSPLNM_strct = {};
OFDIRS_strct.PW = {}; bsnm = 'validation__lns_abc_PW';
OFDIRS_strct.PW.LAR = [SPCFMEXloc,bsnm,LAR_ext]; DSPLNM_strct.LAR = LarRNam; LS_strct.LAR = '-';
OFDIRS_strct.PW.FAF = [SPCFMEXloc,bsnm,FAF_ext]; DSPLNM_strct.FAF = FafRNam; LS_strct.FAF = '--';
OFDIRS_strct.PW.BUF = [SPCFMEXloc,bsnm,BUF_ext]; DSPLNM_strct.BUF = BufRNam; LS_strct.BUF = ':';
OFDIRS_strct.PS = {}; bsnm = 'validation__lns_abc_PS';
OFDIRS_strct.PS.LAR = [SPCFMEXloc,bsnm,LAR_ext];
OFDIRS_strct.PS.FAF = [SPCFMEXloc,bsnm,FAF_ext];
OFDIRS_strct.PS.BUF = [SPCFMEXloc,bsnm,BUF_ext];
OFDIRS_strct.WPS = {}; bsnm = 'validation__lns_abc_WPS';
OFDIRS_strct.WPS.LAR = [SPCFMEXloc,bsnm,LAR_ext];
OFDIRS_strct.WPS.FAF = [SPCFMEXloc,bsnm,FAF_ext];
OFDIRS_strct.WPS.BUF = [SPCFMEXloc,bsnm,BUF_ext];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot energy
if(compareEnergies)
  [fE, axxx] = ABC_plot_NRJs(OFDIRS_strct,DSPLNM_strct,LS_strct,colours_runs,energy_outpath);
end

% compare time series
if(compareStations)
  ABC_plot_TS();
  customSaveFig(fh_TS, TS_outpath, {'fig', 'eps'}, 9999);
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
