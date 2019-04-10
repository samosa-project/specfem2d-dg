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
set(0, 'DefaultLineLineWidth', 1.5); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'defaultTextFontSize', 12); set(0, 'defaultAxesFontSize', 12);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools'); % truncToShortest, readAndSubsampleSynth
addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new');  % plot_total_energy
addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/standalone'); % plot_time_v_dist
SPCFMEXloc = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputDir = [SPCFMEXloc, 'LNSABC_info/']; % don't forget ending '/'

subsample = 1;
subsample_dt = 1e-3;

compareEnergies = 1;
compareStations = 1;

basename = 'LNSABC_PW'; prefix='PW'; distType_x0z1d2=1; errorFactor=1e2; % factor by which multiply the absolute error between methods
% basename = 'LNSABC_PS'; prefix='PS'; distType_x0z1d2=0; errorFactor=1e2; % factor by which multiply the absolute error between methods
% basename = 'LNSABC_WPS'; prefix='WPS'; distType_x0z1d2=2; errorFactor=1e2; % factor by which multiply the absolute error between methods

% bufferRunOF   = [SPCFMEXloc,basename,'_buffers/OUTPUT_FILES_eps1em4']; suffixBu='$\varepsilon=10^{-4}$';
bufferRunOF   = [SPCFMEXloc,basename,'_buffers/OUTPUT_FILES_eps0p25']; suffixBu='$\varepsilon=0.25$';
% bufferRunOF   = [SPCFMEXloc,basename,'_buffers/OUTPUT_FILES_eps0p5']; suffixBu='$\varepsilon=0.5$';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build OUTPUT_FILES directories' paths
BufRNam='buffers';
largeRunOF    = [SPCFMEXloc,basename,'_large/OUTPUT_FILES_baseline']; LarRNam='large';
farFieldRunOF = [SPCFMEXloc,basename,'_FF/OUTPUT_FILES_baseline']; FafRNam='far-field';
% safeguard
if(numel(bufferRunOF)>0); if(not(strcmp(bufferRunOF(end),filesep))); bufferRunOF=[bufferRunOF,filesep]; end; end;
if(numel(largeRunOF)>0); if(not(strcmp(largeRunOF(end),filesep))); largeRunOF=[largeRunOF,filesep]; end; end;
if(numel(farFieldRunOF)>0); if(not(strcmp(farFieldRunOF(end),filesep))); farFieldRunOF=[farFieldRunOF,filesep]; end; end;
% retrieve parameters
OF=largeRunOF; [xmM_La, zmM_La, ~] = getRunParameters([OF,'input_parfile'], [OF,'SOURCE'], [OF,'input_interfaces']);
OF=bufferRunOF; [xmM_Bu, zmM_Bu, ~] = getRunParameters([OF,'input_parfile'], [OF,'SOURCE'], [OF,'input_interfaces']); LBUF=extractParamFromInputFile([OF,'input_parfile'], 'ABC_STRETCH_TOP_LBUF', 'float');
OF=farFieldRunOF; [xmM_Fa, zmM_Fa, ~] = getRunParameters([OF,'input_parfile'], [OF,'SOURCE'], [OF,'input_interfaces']);
% prepare information strings
tit_La = ['$[',num2str(min(xmM_La)),', ',num2str(max(xmM_La)),']\times[',num2str(min(zmM_La)),', ',num2str(max(zmM_La)),']$'];
tit_Bu = ['$[',num2str(min(xmM_Bu)),', ',num2str(max(xmM_Bu)),']\times[',num2str(min(zmM_Bu)),', ',num2str(max(zmM_Bu)),']$, Buffer=',num2str(LBUF),' [m], ',suffixBu];
tit_Fa = ['$[',num2str(min(xmM_Fa)),', ',num2str(max(zmM_Fa)),']\times[',num2str(min(zmM_Fa)),', ',num2str(max(zmM_Fa)),']$'];
% tit_La = ['$z\in[',num2str(min(zmM_La)),', ',num2str(max(zmM_La)),']$'];
% tit_Bu = ['$z\in[',num2str(min(zmM_Bu)),', ',num2str(max(zmM_Bu)),']$, Buffer=',num2str(LBUF),', ',suffixBu];
% tit_Fa = ['$z\in[',num2str(min(zmM_Fa)),', ',num2str(max(zmM_Fa)),']$'];
tit_all = {[LarRNam,'=[',tit_La,'], '],[BufRNam,'=[',tit_Bu,'], '],[FafRNam,'=[',tit_Fa,']']};

outputFigDir_wprefix=[outputDir,prefix,'_'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot energy
if(compareEnergies)
  outputFigPath=[outputFigDir_wprefix,'KE'];
  fKE = plot_total_energy({largeRunOF,farFieldRunOF,bufferRunOF},0,1,tit_all);
  customSaveFig(outputFigPath, {'fig', 'jpg', 'eps'});
  
  outputFigPath=[outputFigDir_wprefix,'PE'];
  fPE = plot_total_energy({largeRunOF,farFieldRunOF,bufferRunOF},1,0,tit_all);
  customSaveFig(outputFigPath, {'fig', 'jpg', 'eps'});
end

% plot time series
if(compareStations)
  [~,diffSTATBufferLarge]=system(['diff ',bufferRunOF,'STATIONS ',largeRunOF,'STATIONS']);
  [~,diffSTATFFLarge]=system(['diff ',farFieldRunOF,'STATIONS ',largeRunOF,'STATIONS']);
  if(numel(diffSTATBufferLarge)>0)
    error('STATIONS were not the same between runs buffer and large');
  end
  if(numel(diffSTATFFLarge)>0)
    error('STATIONS were not the same between runs farfield and large');
  end
  
%   STATIONS = [bufferRunOF,'STATIONS'];
  [x_stat, z_stat, y_stat, stat_file] = loadStations(bufferRunOF);
  nbstats=numel(x_stat);
  
  % Ordering for plot. Lowest is below (towards background), higher is above (towards foreground).
  orEFf=0;
  orEBu=1;
  orLa=2;
  orFf=3;
  orBu=4;
  
  % retrieve time series
  Zamp=[]; NMs={};
  for i=1:nbstats
    j=orLa*nbstats+i; iLa=j; % load LARGE run (assumed to be truth)
    [data, ~] = readAndSubsampleSynth(largeRunOF, i, 'BXZ', 'semv', subsample, subsample_dt, 1); Zamp(j,:)=data(:,2)'; if(i==1); time=data(:,1)'; end;
    NMs{j} = ['S',num2str(i),' ',LarRNam,'']; COLs{j} = 'k'; LSs{j} = '-';
    
    j=orBu*nbstats+i; iBUF=j; % load BUFFER run (or any other test ABC method)
    [data, ~] = readAndSubsampleSynth(bufferRunOF, i, 'BXZ', 'semv', subsample, subsample_dt, 1); Zamp(j,:)=data(:,2)';
    NMs{j} = ['S',num2str(i),' ',BufRNam,'']; COLs{j} = 'm'; LSs{j} = '--';
    j=orEBu*nbstats+i; % save error between the previous and LARGE
    Zamp(j,:) =  errorFactor* abs(Zamp(iLa,:) - Zamp(iBUF,:));
    NMs{j} = ['S',num2str(i),' $',num2str(errorFactor),'\times|$',LarRNam,'$-$',BufRNam,'$|$']; COLs{j} = 'm'; LSs{j} = ':';
    
    j=orFf*nbstats+i; iFAF=j; % load FARFIELD run (or any other test ABC method)
    [data, ~] = readAndSubsampleSynth(farFieldRunOF, i, 'BXZ', 'semv', subsample, subsample_dt, 1); Zamp(j,:)=data(:,2)';
    NMs{j} = ['S',num2str(i),' ',FafRNam,'']; COLs{j} = 'c'; LSs{j} = '--';
    j=orEFf*nbstats+i; % save error between the previous and LARGE
    Zamp(j,:) =  errorFactor* abs(Zamp(iLa,:) - Zamp(iFAF,:));
    NMs{j} = ['S',num2str(i),' $',num2str(errorFactor),'\times|$',LarRNam,'$-$',FafRNam,'$|$']; COLs{j} = 'c'; LSs{j} = ':';
  end
  
  % actual plot
  outputFigPath=[outputFigDir_wprefix,'time_series'];
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
      error('distance swtich not implemented');
  end
  fTvD = plot_time_v_dist(repmat(time,[nbstats*nRepetitions,1]), Zamp, repmat(distttt,[nRepetitions,1]), 0, tit_all, distSymbol, 1, NMs, COLs, LSs);
  % cosmetics
  figure(fTvD);
  cax=gca();
  set(cax.Children(1:nRepetitions*(nbstats-1)),'handlevisibility','off');
  legend('location','northwest');
  prettyAxes(fTvD);
  customSaveFig(outputFigPath, {'fig', 'jpg', 'eps'});
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