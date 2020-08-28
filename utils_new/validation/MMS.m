% Author:        Léo Martire.
% Description:   Verifies LNS simulations through the mean of manufactured
%                solutions.
% Notes:         TODO.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

clear all;
% close all;
clc;

verbose = 0;
addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new'));
SPCFM_EX_DIR = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefix = 'validation__lns_manufactured';
savefigpath = [SPCFM_EX_DIR,prefix,'_results/'];
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
% % OFDs = {[SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_dxN=00250_mu_dp__1p25em5']}; IDzs={[2]*10000}; plotFields_do = 0;
% OFDs = {[SPCFM_EX_DIR,prefix,'_N=0025_iv/OUTPUT_FILES']}; IDzs={[100e3]}; plotFields_do = 1;
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
    OFDs = {[SPCFM_EX_DIR,prefix,'_N=0100_iv/OUTPUT_FILES']}; testCase = 'inviscid';
    iterationsToPlot_forEachOFD = {[commonID2*100]};
    plotFields_do = 1; % activate plotting of fields
  case 1
    OFDs = {[SPCFM_EX_DIR,prefix,'_N=0001_iv/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0005_iv/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0025_iv/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0050_iv/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0100_iv/OUTPUT_FILES']}; testCase = 'inviscid';
    iterationsToPlot_forEachOFD = {[commonID2], [commonID2*5], [commonID2*25], [commonID2*50], [commonID2*100]};
    plotFields_do = 0; % deactivate plotting of fields
  case 2
    OFDs = {[SPCFM_EX_DIR,prefix,'_N=0001_ka/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0005_ka/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0025_ka/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0050_ka/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0100_ka/OUTPUT_FILES']}; testCase = 'kappa';
    iterationsToPlot_forEachOFD = {[commonID2], [commonID2*5], [commonID2*25], [commonID2*50], [commonID2*100]};
    plotFields_do = 0; % deactivate plotting of fields
  case 3
    OFDs = {[SPCFM_EX_DIR,prefix,'_N=0001_mu/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0002_mu/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0003_mu/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0004_mu/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0005_mu/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0025_mu/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0050_mu/OUTPUT_FILES'], ...
            [SPCFM_EX_DIR,prefix,'_N=0100_mu/OUTPUT_FILES']}; testCase = 'mu';
    iterationsToPlot_forEachOFD = {[commonID2], [commonID2*2], [commonID2*3], [commonID2*4], [commonID2*5], [commonID2*25], [commonID2*50], [commonID2*100]};
    plotFields_do = 0; % deactivate plotting of fields
  case 4
    % illustrative plotfield
    OFDs = {[SPCFM_EX_DIR,prefix,'_N=0025_iv/OUTPUT_FILES']}; testCase = 'inviscid';
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
