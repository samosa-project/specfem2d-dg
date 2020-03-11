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

subsample = 1;
subsample_dt = 1e-4;

compareEnergies = 1;
compareStations = 0;

% basename = 'LNSABC_PW'; prefix='PW'; distType_x0z1d2=1; errorFactor=1e5; % factor by which multiply the absolute error between methods
% basename = 'LNSABC_PS'; prefix='PS'; distType_x0z1d2=0; errorFactor=1e2; % factor by which multiply the absolute error between methods
basename = 'LNSABC_WPS'; prefix='WPS'; distType_x0z1d2=2; errorFactor=1e2; % factor by which multiply the absolute error between methods

% bufferRunOF   = [SPCFMEXloc,basename,'_buffers/OUTPUT_FILES_eps0p25']; suffixBu='$\varepsilon=0.25$';
% bufferRunOF   = [SPCFMEXloc,basename,'_buffers/OUTPUT_FILES_p3p25_q6_e1em4']; suffixBu='$\varepsilon=10^{-4}$';
% bufferRunOF   = [SPCFMEXloc,basename,'_buffers/OUTPUT_FILES_p3p25_q6_e0p250']; suffixBu='$\varepsilon=0.25$';
% bufferRunOF   = [SPCFMEXloc,basename,'_buffers/OUTPUT_FILES_p3p25_q6_e0p160_woE']; suffixBu='$\varepsilon=0.16$';
% bufferRunOF   = [SPCFMEXloc,basename,'_buffers/OUTPUT_FILES_p3p25_q6_e0p001']; suffixBu='$\varepsilon=10^{-3}$';
% bufferRunOF   = [SPCFMEXloc,basename,'_buffers/OUTPUT_FILES_p3p25_q6_e0p500']; suffixBu='$\varepsilon=0.5$';

bufferRunOF   = [SPCFMEXloc,basename,'_buffers/OUTPUT_FILES_p3p25q6e0p001']; suffixBu='$\varepsilon=10^{-3}$';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build OUTPUT_FILES directories' paths
BufRNam='buffers';
largeRunOF    = [SPCFMEXloc,basename,'_large/OUTPUT_FILES_baseline']; LarRNam='large';
farFieldRunOF = [SPCFMEXloc,basename,'_FF/OUTPUT_FILES_FF']; FafRNam='far-field';

% Locations of OUTPUT_FILES directories. Should be the same across all test cases.
LAR_ext = '_large/OUTPUT_FILES_baseline/';
FAF_ext = '_FF/OUTPUT_FILES_FF/';
BUF_ext = '_buffers/OUTPUT_FILES_p3p25q6e0p001/';

OFDIRS_strct = {}; DSPLNM_strct = {};
OFDIRS_strct.PW = {}; basename = 'LNSABC_PW';
OFDIRS_strct.PW.LAR = [SPCFMEXloc,basename,LAR_ext]; DSPLNM_strct.LAR = LarRNam; LS_strct.LAR = '-';
OFDIRS_strct.PW.FAF = [SPCFMEXloc,basename,FAF_ext]; DSPLNM_strct.FAF = FafRNam; LS_strct.FAF = '--';
OFDIRS_strct.PW.BUF = [SPCFMEXloc,basename,BUF_ext]; DSPLNM_strct.BUF = BufRNam; LS_strct.BUF = ':';
OFDIRS_strct.PS = {}; basename = 'LNSABC_PS';
OFDIRS_strct.PS.LAR = [SPCFMEXloc,basename,LAR_ext];
OFDIRS_strct.PS.FAF = [SPCFMEXloc,basename,FAF_ext];
OFDIRS_strct.PS.BUF = [SPCFMEXloc,basename,BUF_ext];
OFDIRS_strct.WPS = {}; basename = 'LNSABC_WPS';
OFDIRS_strct.WPS.LAR = [SPCFMEXloc,basename,LAR_ext];
OFDIRS_strct.WPS.FAF = [SPCFMEXloc,basename,FAF_ext];
OFDIRS_strct.WPS.BUF = [SPCFMEXloc,basename,BUF_ext];

% safeguard
if(numel(bufferRunOF)>0); if(not(strcmp(bufferRunOF(end),filesep))); bufferRunOF=[bufferRunOF,filesep]; end; end;
if(numel(largeRunOF)>0); if(not(strcmp(largeRunOF(end),filesep))); largeRunOF=[largeRunOF,filesep]; end; end;
if(numel(farFieldRunOF)>0); if(not(strcmp(farFieldRunOF(end),filesep))); farFieldRunOF=[farFieldRunOF,filesep]; end; end;
% retrieve parameters
OF=largeRunOF; [xmM_La, zmM_La, ~] = readExampleFiles([OF,'input_parfile'], [OF,'SOURCE'], [OF,'input_interfaces']);
OF=bufferRunOF; [xmM_Bu, zmM_Bu, ~] = readExampleFiles([OF,'input_parfile'], [OF,'SOURCE'], [OF,'input_interfaces']); LBUF=readExampleFiles_extractParam([OF,'input_parfile'], 'ABC_STRETCH_TOP_LBUF', 'float');
OF=farFieldRunOF; [xmM_Fa, zmM_Fa, ~] = readExampleFiles([OF,'input_parfile'], [OF,'SOURCE'], [OF,'input_interfaces']);
% prepare information strings
tit_La = ['$[',num2str(min(xmM_La)),', ',num2str(max(xmM_La)),']\times[',num2str(min(zmM_La)),', ',num2str(max(zmM_La)),']$'];
tit_Bu = ['$[',num2str(min(xmM_Bu)),', ',num2str(max(xmM_Bu)),']\times[',num2str(min(zmM_Bu)),', ',num2str(max(zmM_Bu)),']$, buffer=',num2str(LBUF),' [m], ',suffixBu];
tit_Fa = ['$[',num2str(min(xmM_Fa)),', ',num2str(max(xmM_Fa)),']\times[',num2str(min(zmM_Fa)),', ',num2str(max(zmM_Fa)),']$'];
% tit_La = ['$z\in[',num2str(min(zmM_La)),', ',num2str(max(zmM_La)),']$'];
% tit_Bu = ['$z\in[',num2str(min(zmM_Bu)),', ',num2str(max(zmM_Bu)),']$, Buffer=',num2str(LBUF),', ',suffixBu];
% tit_Fa = ['$z\in[',num2str(min(zmM_Fa)),', ',num2str(max(zmM_Fa)),']$'];
% get energyboxes
ENERGYBOX_La = [readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_XMIN', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_XMAX', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_ZMIN', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_ZMAX', 'float')];
ENERGYBOX_Bu = [readExampleFiles_extractParam([bufferRunOF,'input_parfile'], 'ENERGYBOX_XMIN', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_XMAX', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_ZMIN', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_ZMAX', 'float')];
ENERGYBOX_Fa = [readExampleFiles_extractParam([farFieldRunOF,'input_parfile'], 'ENERGYBOX_XMIN', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_XMAX', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_ZMIN', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_ZMAX', 'float')];
if(not(all(ENERGYBOX_La==ENERGYBOX_Bu & ENERGYBOX_Bu==ENERGYBOX_Fa)))
  error(['ENERGY BOXES WERE NOT THE SAME, ABORTING']);
else
  ENERGYBOX = ENERGYBOX_La;
  globxmin=min([min(xmM_La),min(xmM_Bu),min(xmM_Fa)]);
  globxmax=max([max(xmM_La),max(xmM_Bu),max(xmM_Fa)]);
  globzmin=min([min(zmM_La),min(zmM_Bu),min(zmM_Fa)]);
  globzmax=max([max(zmM_La),max(zmM_Bu),max(zmM_Fa)]);
  if(ENERGYBOX(1)<globxmin)
    ENERGYBOX(1)=globxmin;
  end
  if(ENERGYBOX(2)>globxmax)
    ENERGYBOX(2)=globxmax;
  end
  if(ENERGYBOX(3)<globzmin)
    ENERGYBOX(3)=globzmin;
  end
  if(ENERGYBOX(4)>globzmax)
    ENERGYBOX(4)=globzmax;
  end
end
INFO_all = {[LarRNam,'=[',tit_La,'], ',FafRNam,'=[',tit_Fa,']'],[BufRNam,'=[',tit_Bu,'], '],['energy box: ',['$[',num2str(min(ENERGYBOX(1:2))),', ',num2str(max(ENERGYBOX(1:2))),']\times[',num2str(min(ENERGYBOX(3:4))),', ',num2str(max(ENERGYBOX(3:4))),']$']]};
FIGTITLE = prefix;

displaynames_large_ff_buf = {[LarRNam], [FafRNam], [BufRNam]};

outputFigDir_wprefix=[outputDir,prefix,'_'];

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
  [fE, axxx] = ABC_plot_NRJs(OFDIRS_strct,DSPLNM_strct,LS_strct,);
  %TODO: export descriptions (INFO_all) to a text file
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
      error('distance swtich not implemented');
  end
  fTvD = plot_time_v_dist(repmat(time,[nbstats*nRepetitions,1]), Zamp, repmat(distttt,[nRepetitions,1]), 0, FIGTITLE, distSymbol, 1, NMs, COLs, LSs);
  % cosmetics
  figure(fTvD);
  cax=gca();
  set(cax.Children(1:nRepetitions*(nbstats-1)),'handlevisibility','off');
  legend('location','northwest');
  prettyAxes(fTvD);
  customSaveFig(fTvD, TS_outputFigPath, {'fig', 'jpg', 'eps'});
  disp(["figure(fTvD.Number); customSaveFig(TS_outputFigPath, {'fig', 'jpg', 'eps'});"]);
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