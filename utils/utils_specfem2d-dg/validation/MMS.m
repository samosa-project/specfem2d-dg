% Author:        Léo Martire.
% Description:   Verifies LNS simulations through the mean of manufactured
%                solutions.
% Notes:         TODO.
%
% Usage:
%   TODO.
% with:
%   N. A.
% yields:
%   N. A.

clear all;
% close all;
clc;

[SPCFMEXloc] = setup_overall();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 0;
prefix = 'validation__lns_manufactured';
savefigpath = [SPCFMEXloc,prefix,'_results/'];
err_l2.name = '$L^2$ Error';
err_l2.symbol = ['varepsilon'];
err_rel.name = 'Relative Error';
err_rel.symbol = ['epsilon'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % IDs of dumps to plot.
% IDz = 100;
% % OUTPUT_FILES directory to analyse.
% % OFDs = {[SPCFMEXloc,prefix,'/OUTPUT_FILES_dxN=00250_mu_dp__1p25em5']}; IDzs={[2]*10000}; plotFields_do = 0;
% OFDs = {[SPCFMEXloc,prefix,'_N=0025_iv/OUTPUT_FILES']}; IDzs={[100e3]}; plotFields_do = 1;
% % Test case to plot/compute.
% testCase = 'inviscid';
% % testCase = 'kappa';
% % testCase = 'mu';

% IDs of dumps to plot packed, for paper.
commonID2 = 4000;
generateCase = 4;
switch(generateCase)
  case 5
    % test case
    OFDs = {[SPCFMEXloc,prefix,'_N=0100_iv/OUTPUT_FILES']}; testCase = 'inviscid';
    iterationsToPlot_forEachOFD = {[commonID2*100]};
    plotFields_do = 1; % activate plotting of fields
  case 1
    OFDs = {[SPCFMEXloc,prefix,'_N=0001_iv/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0005_iv/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0025_iv/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0050_iv/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0100_iv/OUTPUT_FILES']}; testCase = 'inviscid';
    iterationsToPlot_forEachOFD = {[commonID2], [commonID2*5], [commonID2*25], [commonID2*50], [commonID2*100]};
    plotFields_do = 0; % deactivate plotting of fields
  case 2
    OFDs = {[SPCFMEXloc,prefix,'_N=0001_ka/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0005_ka/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0025_ka/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0050_ka/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0100_ka/OUTPUT_FILES']}; testCase = 'kappa';
    iterationsToPlot_forEachOFD = {[commonID2], [commonID2*5], [commonID2*25], [commonID2*50], [commonID2*100]};
    plotFields_do = 0; % deactivate plotting of fields
  case 3
    OFDs = {[SPCFMEXloc,prefix,'_N=0001_mu/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0002_mu/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0003_mu/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0004_mu/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0005_mu/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0025_mu/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0050_mu/OUTPUT_FILES'], ...
            [SPCFMEXloc,prefix,'_N=0100_mu/OUTPUT_FILES']}; testCase = 'mu';
    iterationsToPlot_forEachOFD = {[commonID2], [commonID2*2], [commonID2*3], [commonID2*4], [commonID2*5], [commonID2*25], [commonID2*50], [commonID2*100]};
    plotFields_do = 0; % deactivate plotting of fields
  case 4
    % illustrative plotfield
    OFDs = {[SPCFMEXloc,prefix,'_N=0025_iv/OUTPUT_FILES']}; testCase = 'inviscid';
    iterationsToPlot_forEachOFD = {[commonID2*25]};
    plotFields_do = 1; % activate plotting of fields
end

[errQtity, err_l2, err_rel, globalSave, IDNX, IDDT, IDIT, IDEPS, IDRELERR] = MMS_oneCase(prefix, testCase, OFDs, iterationsToPlot_forEachOFD, plotFields_do, err_l2, err_rel, verbose, savefigpath);

%%%%%%%%%%%%%%%%
% Error plots.
%%%%%%%%%%%%%%%%

% Error plots.
plotProgrWRTTime = 0; % plot progression wrt time?

MMS_masterFigure(testCase, errQtity, err_l2, err_rel, globalSave, IDNX, IDDT, IDIT, IDEPS, IDRELERR, savefigpath, plotProgrWRTTime);
