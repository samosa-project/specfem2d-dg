% Author:        Léo Martire.
% Description:   Estimate new run time based on previous runs' saved data.
% Last modified: See file metadata.
% Usage:         Fill section "Parameters" according to the calculation to
%                be performed.
%                Use the scripts in './utils_new/extract_run_information'
%                to fill data used for estimation (function 'load' below).
% Notes:         The method (N-D nearest point search) is approximate and
%                does not take into account the fact that some parameters
%                have more impact than others. All in all, it is a very
%                rough and very approximate method.

% clear all
% close all
clc
format longG;

addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools/');
[data, t, info, ~]=load_data(); % Load data (see function below).

plot_rate = 1;
plotrate_selfulldgonly = 1;
plotrate_snappercentthresh = 10; % in [%]
plotrateversion = 'FNS';
% plotrateversion='LNS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(not(plot_rate))
  % Choices. %%%%%%%%%%%%%%%%%%%%
  % nproc         = 8*48;
  % nbeltsElastic = 6600;
  nbeltsElastic = input(['[',mfilename,'] Number of viscoelastic elements? > ']);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Auto-read from parfile. %%%%%
  parfile='';
  while(not(exist(parfile)==2))
    parfile = input(['[',mfilename,'] Path to simulation parfile > '],'s');
  end
  nstations = sum(readExampleFiles_extractParam(parfile, 'nrec', 'float'));
  nstepseismo = readExampleFiles_extractParam(parfile, 'NSTEP_BETWEEN_OUTPUT_SEISMOS', 'float');
  extmesh = readExampleFiles_extractParam(parfile, 'read_external_mesh', 'bool');
  if(not(extmesh))
    command = ['cat ',parfile,' | grep -oP " *1 +[0-9]+ +[0-9]+ +[0-9]+ | tail -1"'];
    [~, r]=system(command);
    r = str2num(r);
    neltot    = r(end,end-2)*r(end,end);
  else
    disp(['[',mfilename,'] Need to read number of elements from the external mesh.']);
    spl=split(parfile,'/');
    spl{end}='';
    folderofparfile=strjoin(spl,filesep);
    tryMeshExtMeshLocalisation = [folderofparfile,'EXTMSH',filesep,'Mesh_extMesh'];
    if(not(exist(tryMeshExtMeshLocalisation)))
      error('ENTER THE NUMBER OF POINTS OMEGALUL');
    else
      [~,r]=system(['head ',tryMeshExtMeshLocalisation,' -n 1']);
      neltot    = str2num(r);
      disp(['[',mfilename,'] Total number of elements found to be ',num2str(neltot),'. Check your EXTMSH files.']);
      disp(['[',mfilename,'] Press any key to continue.']);
      pause();
    end
  end
  neldg       = neltot-nbeltsElastic;
  nstepsnap   = readExampleFiles_extractParam(parfile, 'NSTEP_BETWEEN_OUTPUT_IMAGES', 'float');
  nsteptot    = readExampleFiles_extractParam(parfile, 'NSTEP', 'float');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Define by hand. %%%%%%%%%%%%%
  % nstations   = 857;
  % nstepseismo = 100;
  % nstepsnap   = 5000;
  % nsteptot    = 70000;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Display data.               %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % for i=numel(t)-3:numel(t)
  for i=1:numel(t)
    ID=info{i}(1); ID=ID{1};
    disp(['[',mfilename,'] Run ID ',sprintf('%7d',ID),' (',sprintf('%5.1f',data(i,2)*100),'% DG):       ',sprintf('%.2e',t(i)),' s CPU, per el., per it.       (''',char(string(info{i}(2))),''').']);
  end
  disp(' ');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Computation of estimate.    %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  weigth=[1,100,25,1,50]; % Importance of each parameter.
  meandata=mean(data,1);

  nproc = -1;
  while(nproc~=0)
    nproc = -1;
    while(not(numel(nproc)==1 & nproc>=0))
      disp(' ');disp(' ');disp(' ');
      nproc = input(['[',mfilename,'] Number of procs (>=0, 0 to stop script)? > ']);
    end
    if(nproc==0)
      continue
    end
    point=[nstations,neldg/neltot,nstepsnap/nsteptot,nstepseismo/nsteptot,neltot/nproc]; % Format point as data format. Currently [% elements as DG, % timesteps as snapshots, elements per proc].

    disp(['[',mfilename,'] ',num2str(neltot),' elements including ',num2str(neldg),' DG elements. ',num2str(nsteptot),' time steps. ',num2str(nstations),' stations sampling every ',num2str(nstepseismo),' iterations. Snapshots taken every ',num2str(nstepsnap),' iterations. ',num2str(nproc),' CPUs.']);
    disp(" ");
    disp(['[',mfilename,']                     [        n_stations       percent_DG     percent_snap    percent_synth n_elems_per_proc]']);
    disp(strcat("[",mfilename,"] Current point:      [ ", sprintf("%17.3e", point), "]."));
    % idp=dsearchn(data,point);
    idp=dsearchn((data-meandata)./weigth,(point-meandata)./weigth); % Optimised search.
    disp(strcat("[",mfilename,"] Closest data point: [ ", sprintf("%17.3e", data(idp,:)), "] (",info{idp}{2},")."));
    disp(strcat("[",mfilename,"] Weights:            [ ", sprintf("%17.f", weigth), "]."));
    time=t(idp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display.                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(" ");

    cputime=time*neltot*nsteptot;
    % hms=fix(mod(cputime,[0,3600,60])./[3600,60,1]);
    % cputimestr=strcat(sprintf('%5.f',hms(1)),"h ",sprintf('%2.f',hms(2)),"m ",sprintf('%2.f',hms(3)),'s');
    cputimestr=formatSecondsToHHMMSS(cputime);

    realtime=time*neltot*nsteptot/nproc;
    % hms=fix(mod(realtime,[0,3600,60])./[3600,60,1]);
    % realtimestr=strcat(sprintf('%5.f',hms(1)),"h ",sprintf('%2.f',hms(2)),"m ",sprintf('%2.f',hms(3)),'s');
    realtimestr=formatSecondsToHHMMSS(realtime);

    disp(strcat("[",mfilename,"] Expected time per element, per iteration:          ",sprintf("%.3e", time), " s."));
    disp(strcat("[",mfilename,"] Expected run time:                        ",cputimestr, " (CPU), i.e.  ",sprintf('%15.0f',cputime)," s."));
    disp(strcat("[",mfilename,"]                                           ",realtimestr, " (real), i.e. ",sprintf('%15.0f',realtime)," s."));
    if(realtime>86400)
      disp(['[',mfilename,', WARNING] Expected run time is over one (1) day. Most HPC centres do not allow runs this long. To drop below 1 day, try out ',num2str(ceil(ceil(nproc*realtime/86400)/16)*16),' CPUs.']);
    end
    disp(['[',mfilename,', WARNING] Recall the method used for estimation is very rough and approximate. Do not take the estimation for granted.']);
  end
  
else % (if(plot_rate))
  disp(['[',mfilename,'] Plotting parallelisation rate.']);
  set(0, 'DefaultLineLineWidth', 2); % Default at 0.5.
  set(0, 'DefaultLineMarkerSize', 20); % Default at 6.
  set(0, 'defaultTextFontSize', 24);
  set(0, 'defaultAxesFontSize', 24); % Default at 10.
  set(0, 'DefaultTextInterpreter', 'latex');
  set(0, 'DefaultLegendInterpreter', 'latex');

  [data, t, info, dataraw] = load_data(); % Load data (see function below).
  % t = total CPU time [s]
  
  idX_SumNrec          = 1;
  idX_dgpercent        = 2;
  idX_snappercent      = 3;
  idX_synthpercent     = 4;
  idX_nbeltspproc      = 5;
  
  idRaw_nbelts         = 1;
  idRaw_nbeltsdg       = 2;
  idRaw_Cfl            = 3;
  idRaw_NSTEP          = 4;
  idRaw_NPROC          = 5;
  idRaw_NstepOutImages = 6;
  idRaw_SumNrec        = 7;
  idRaw_NstepOutSeismo = 8;
  idRaw_realtime       = 9;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Create selection.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sel = [];
  if(plotrate_selfulldgonly)
    sel = (data(:, idRaw_nbeltsdg)==1); % Select full DG only.
    disp(['[',mfilename,'] Selected only DG data points.']);
  else
    sel = logical(ones(size(t))); % Select all.
    disp(['[',mfilename,'] Selected all data points.']);
  end
  
  sel = sel & (data(:,idX_snappercent) < plotrate_snappercentthresh/100);
  disp(['[',mfilename,'] Restricted selection to only data points with snapshot frequency < ',num2str(plotrate_snappercentthresh),' %.']);
  
  % Remove some points from selection.
  if(strcmp(plotrateversion,'FNS'))
    % Remove all LNS data.
    for is=1:size(t); if(sel(is)); txt=info{is}(2); txt=txt{1}; if(not(isempty(regexp(txt,'LNS')))); sel(is)=0; end; end; end
    disp(['[',mfilename,'] Restricted selection to only FNS data points.']);
  elseif(strcmp(plotrateversion,'LNS'))
    % Remove all FNS data (select only LNS).
    for is=1:size(t); if(sel(is)); txt=info{is}(2); txt=txt{1}; if(isempty(regexp(txt,'LNS'))); sel(is)=0; end; end; end
    disp(['[',mfilename,'] Restricted selection to only LNS data points.']);
  else
    error('kek');
  end
  for is=1:size(t); if(sel(is)); txt=info{is}(2); txt=txt{1}; if(not(isempty(regexp(txt,'mb gmsh')))); txt, sel(is)=0; end; end; end
  for is=1:size(t); if(sel(is)); txt=info{is}(2); txt=txt{1}; if(not(isempty(regexp(txt,'mb huge')))); txt, sel(is)=0; end; end; end
  disp(['[',mfilename,'] Restricted selection to non-"mb gmsh" and non-"mb huge" runs.']);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Prepare plot.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  x = dataraw(sel, idRaw_NPROC); xlab=['number of CPUs $n$ [1]']; % Get number of CPUs.
  
%   y = t(sel); ylab=['CPU time, per element, per iteration [s]'];
%   theoretical_pow = 0; theoretical_lab=['theoretical best: ??'];
%   y = t(sel) ./ dataraw(sel, idRaw_NPROC); ylab={['mean CPU time per element per iteration [s]'],['(total CPU time / number of CPUS)']}; % Get real time (total CPU time p. el. p. it. / CPU).
%   theoretical_pow = -1; theoretical_lab=['theoretical best: $\propto n^{',num2str(theoretical_pow),'}$'];
%   y = t(sel) .* dataraw(sel, idRaw_nbelts) .* dataraw(sel,5); % Get CPU time for each simulation (CPU time p. el. p. it. * nbel.).
%   theoretical_pow = 0; theoretical_lab=['theoretical best: ??'];
  y = t(sel); ylab={['total CPU time per element per iteration [s]'], ['(real elapsed time $\times$ number of CPUS)']}; % Get real time (total CPU time p. el. p. it. / CPU).
  theoretical_pow = 0; theoretical_lab=['theoretical best: constant ($\propto n^{',num2str(theoretical_pow),'}$)'];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Choose color-code.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  c = log10(data(sel,3)); clab = ['log$_{10}$(snapshot frequency)']; % Color code with snapshot frequency.
  ux=sort(unique(x));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Build best w.r.t. color.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ubesty=[];
  for ix=1:size(ux)
    sely = y(x == ux(ix));
    selc = c(x == ux(ix));
    ubesty(ix, 1) = min(sely(selc == min(selc)));
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Build linear fits.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  beta_lfit_all = [ones(length(x),  1) log(x) ] \ log(y);
  beta_lfit_bst = [ones(length(ux), 1) log(ux)] \ log(ubesty);
  lfitlab_all = ['linear fit best: $\propto n^{', sprintf('%.2g', beta_lfit_all(2)), '}$'];
  lfitlab_bst = ['linear fit best: $\propto n^{', sprintf('%.2g', beta_lfit_bst(2)), '}$'];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Actual plot.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  f = figure('units','normalized','outerposition', [0 0 1 1]);
  scatter(x, y, 50, c, 'filled', 'displayname', ['data, color $\Leftrightarrow$ ', clab]); hold on;
  
  loglog(ux, exp(beta_lfit_all(1)) * ux.^beta_lfit_all(2), 'k-', 'displayname', ['linear fit all: $\propto n^{', sprintf('%.2g', beta_lfit_all(2)), '}$']);
  
  scatter(ux, ubesty, 'displayname', ['best data for each abscissa w.r.t. color']);
  
  loglog(ux, exp(beta_lfit_bst(1)) * ux.^beta_lfit_bst(2), 'k--','displayname', lfitlab_bst);
  
  loglog(ux, ux.^theoretical_pow * (mean(y(ux == min(ux)))/min(ux)^theoretical_pow), 'k:','displayname', theoretical_lab);
  
  set(gca, 'xscale', 'log');
  set(gca, 'yscale', 'log');
  xlim([0.9*min(ux), 1.1*max(ux)]);
  xlabel(xlab);
  ylabel(ylab);
  
  h = colorbar;
  set(h, 'ticklabelinterpreter','latex');
  legend('location', 'best');
  
  title(['Paralellisation Efficiency - ', plotrateversion]);
  
  prettyAxes(f);
  
%   % WEAK RATE
%   x = dataraw(sel,1); xlab=['number of elements $n$ [1]']; % Get nelems
% %   y = t(sel); ylab=['CPU time, per element, per iteration [s]']; 
%   y = t(sel).*dataraw(sel,1)./dataraw(sel,5); ylab=['real time, per element, per iteration [s]']; % Get real time (CPU time p. el. p. it. / CPU).
% %   y = t(sel).*dataraw(sel,1).*dataraw(sel,5); % Get CPU time for each simulation (CPU time p. el. p. it. * n. el.).
%   c = log(data(sel,3)); clab=['snapshot frequency']; % Color code with snapshot frequency.
%   ux=sort(unique(x));
%   % build best wrt c;
%   ubesty=[];
%   for ix=1:size(ux)
%     sely=y(x==ux(ix)); selc=c(x==ux(ix)); ubesty(ix,1)=min(sely(selc==min(selc))); % cheated a little bit here, chosen best time among times that are at best c
%   end
%   subplot(212);
% %   loglog(x,y,'.','markersize',20, 'displayname','all data');
%   scatter(x,y,50,c,'filled', 'displayname',['data, color $\Leftrightarrow$ log(',clab,')']); hold on;
%   beta=[ones(length(x),1) log(x)]\log(y);
%   loglog(ux,exp(beta(1))*ux.^beta(2), 'k-','displayname', ['linear fit all: $\propto n^{', sprintf('%.2g',beta(2)), '}$']);
%   scatter(ux,ubesty,'displayname',['best data for each abscissa w.r.t. color']);
%   beta=[ones(length(ux),1) log(ux)]\log(ubesty);
%   loglog(ux,exp(beta(1))*ux.^beta(2), 'k--','displayname', ['linear fit best: $\propto n^{', sprintf('%.2g',beta(2)), '}$']);
% %   loglog(ux,ux.^(-1)*mean(y(ux==min(ux)))/min(ux)^(-1), 'k:','displayname', ['theoretical best: $\propto n^{-1}$']);
%   
%   set(gca,'xscale','log'); set(gca,'yscale','log'); xlim([0.9*min(ux), 1.1*max(ux)]); xlabel(xlab); ylabel(ylab);
%   set(gca,'ticklabelinterpreter','latex'); h=colorbar; set(h,'ticklabelinterpreter','latex'); legend('location', 'best');
%   title(['Weak Rate']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function containing data.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, t, RUNINFO, RUN_RAWDATA]=load_data()
  col_nbelts   = 1;
  col_nbeltsdg = 2;
  col_Cfl      = 3;
  col_NSTEP    = 4;
  col_NPROC    = 5;
  col_NstepOutImages=6;
  col_SumNrec  = 7;
  col_NstepOutSeismo=8;
  col_realtime = 9;
  i=1;
  RUN_RAWDATA(i,:)=[   2500    2500 0.404   2000    4   50    0   0     82];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[   2500    2500 0.404   2000   16   50    0   0     24];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[   2500    2500 0.404   2000  256   50    0   0     11];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  10000   10000 0.404   4000    4   50    0   0    920];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  10000   10000 0.404   4000   16   50    0   0    257];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  10000   10000 0.404   4000  256   50    0   0     28];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  40000   40000 0.404   8000    4   50    0   0   8287];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  40000   40000 0.404   8000   16   50    0   0   2326];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  40000   40000 0.404   8000  256   50    0   0    157];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  40000   40000 0.404   8000  512   50    0   0    126];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  40000   40000 0.404   8000 1024   50    0   0    123];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 160000  160000 0.404  16000  256   50    0   0   1361];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 160000  160000 0.404  16000  512   50    0   0    828];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 160000  160000 0.404  16000 1024   50    0   0    576];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 640000  640000 0.404  32000  256   50    0   0  11963];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 640000  640000 0.404  32000  512   50    0   0   6945];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 640000  640000 0.404  32000 1024   50    0   0   4728];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[2560000 2560000 0.404  64000  256    0    0   0  85460];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[2560000 2560000 0.404  64000  512    0    0   0  41534];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[2560000 2560000 0.404  64000 1024    0    0   0  18854];  RUNINFO{i}={-1,'FNS test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[2016000 1830000 0.471 120000  256  200  116  50 124902];  RUNINFO{i}={583041, 'FNS OKQ45'}; i=i+1;
  RUN_RAWDATA(i,:)=[2016000 1830000 0.471 120000  256  200  116  50 123034];  RUNINFO{i}={586984, 'FNS OKQ0'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 227000  150000 0.443  30000  256  500    0   0   2920];  RUNINFO{i}={593959, 'FNS SH soft final'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 227000  150000 0.443  30000  256  500    0   0   2980];  RUNINFO{i}={593960, 'FNS SH hard final'}; i=i+1;
  RUN_RAWDATA(i,:)=[2892000 2832000 0.376  27000  512  200    0   0  23240];  RUNINFO{i}={594536, 'FNS StratoExplo66June1200 - not finished'}; i=i+1;
  RUN_RAWDATA(i,:)=[  13060   13060 0.156  29400   16  200    0   0   5025];  RUNINFO{i}={595104, 'FNS StratoExplo66June1200 - test'}; i=i+1;
  RUN_RAWDATA(i,:)=[2892000 2832000 0.376  40150  512  300    0   0  34470];  RUNINFO{i}={595500, 'FNS StratoExplo66June1200 - not finished'}; i=i+1;
  RUN_RAWDATA(i,:)=[3142000 3082000 0.376 140100  560  300    0   0 110028];  RUNINFO{i}={597316, 'FNS StratoExplo66June1200 - til n=140100'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 227000  150000 0.443  15750  256  500    2  50   1542];  RUNINFO{i}={610736, 'FNS SH soft final redone for data comparison'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 227000  150000 0.443  30000  256  500  102  50   3097];  RUNINFO{i}={616368, 'FNS SH soft final redone for data comparison'}; i=i+1;
  RUN_RAWDATA(i,:)=[  78000       0 0.443  30000  256  500  102  50    795];  RUNINFO{i}={623195, 'FNS SH soft final redone w/o fluid'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 748360  748360 0.345   8000   32  250  114  50  27994];  RUNINFO{i}={624650, 'FNS StratoBaro full with IBF'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 748360  748360 0.345   5750   64  250  114  50   9225];  RUNINFO{i}={637450, 'FNS StratoBaro re-run with EBF'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 748360  748360 0.345  16000   96  250  114  50  17354];  RUNINFO{i}={641616, 'FNS StratoBaro full with EBF'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 266162  266162 0.345  16000  128  250  114  50   5904];  RUNINFO{i}={656505, 'FNS microbaroms periodic with EBF'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 266162  266162 0.345  16000  256  250  114  50   4076];  RUNINFO{i}={655513, 'FNS microbaroms periodic with EBF'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 102150       0 0.443  30000   32  500   80  50    757];  RUNINFO{i}={660223, 'FNS SH soft axisym'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 102150       0 0.443  30000   32  500   80  50    744];  RUNINFO{i}={661601, 'FNS SH hard axisym'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 225000  150000 0.443  12250   32  500   82  50  10802];  RUNINFO{i}={668888, 'FNS SH soft first layer tweaked stopped 12kit'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 224000  150000 0.443  20000   64  500   82  50   7802];  RUNINFO{i}={668888, 'FNS SH soft first layer retweaked'}; i=i+1;
  RUN_RAWDATA(i,:)=[1866000 1680000 0.471 110000  256 5000  116  50 102012];  RUNINFO{i}={668833, 'FNS OKQ0 redone'}; i=i+1;
  RUN_RAWDATA(i,:)=[1866000 1680000 0.471 110000  256 5000  116  50 101666];  RUNINFO{i}={668844, 'FNS OKQ45 redone'}; i=i+1;
  RUN_RAWDATA(i,:)=[3546400 3546400 0.441  15500  480  500  452  25 250007];  RUNINFO{i}={642746, 'FNS mb huge 642746'}; i=i+1;
  RUN_RAWDATA(i,:)=[3627000 3627000 0.441  76500  480  500  452  25 250010];  RUNINFO{i}={672048, 'FNS mb huge 672048'}; i=i+1;
  RUN_RAWDATA(i,:)=[  29100   10500 0.471   2500   64  100    2  50     85];  RUNINFO{i}={71913,  'FNS OKQ45 small for test impedance'}; i=i+1;
  RUN_RAWDATA(i,:)=[  29100       0 0.471   2500   64  100    2  50     47];  RUNINFO{i}={71936,  'FNS OKQ45 small for test impedance but potential'}; i=i+1;
  RUN_RAWDATA(i,:)=[  57887   39553 0.441  28000   16  250   29  25  18010];  RUNINFO{i}={74752,  'FNS tir de mine heavy & incomplete'}; i=i+1;
  RUN_RAWDATA(i,:)=[  17220   17220 0.570  21400   32  500   52  25    1083]; RUNINFO{i}={74710,  'FNS mb gmsh'}; i=i+1;
  RUN_RAWDATA(i,:)=[  19425   13563 0.464  40000   16  250   29  25    6580]; RUNINFO{i}={75040,  'FNS tir de mine light & full'}; i=i+1;
  RUN_RAWDATA(i,:)=[  10000   10000 0.404  10000    4   50    1  25    1669]; RUNINFO{i}={830669, 'FNS visc'}; i=i+1;
  RUN_RAWDATA(i,:)=[  10000   10000 0.404  10000    4   50    1  25    1252]; RUNINFO{i}={830672, 'LNS visc'}; i=i+1; % First LNS. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  RUN_RAWDATA(i,:)=[  10000   10000 0.404  10000    4   50    1  25    1590]; RUNINFO{i}={830670, 'FNS novisc'}; i=i+1;
  RUN_RAWDATA(i,:)=[  10000   10000 0.404  10000    4   50    1  25     809]; RUNINFO{i}={830671, 'LNS novisc'}; i=i+1;
  RUN_RAWDATA(i,:)=[  29171   29171 0.421   4000   32  500   87  25     351]; RUNINFO{i}={120363, 'FNS mb gmsh'}; i=i+1;
  RUN_RAWDATA(i,:)=[  29171   29171 0.421   4000   32  500   87  25     325]; RUNINFO{i}={120414, 'LNS mb gmsh 200km'}; i=i+1;
  RUN_RAWDATA(i,:)=[  53719   53719 0.423  10000   48  500   87  25     931]; RUNINFO{i}={120557, 'LNS mb gmsh 400km (1)'}; i=i+1;
  RUN_RAWDATA(i,:)=[  53719   53719 0.423  10000   48  500   87  25     941]; RUNINFO{i}={120621, 'LNS mb gmsh 400km (2)'}; i=i+1;
  RUN_RAWDATA(i,:)=[7155000 7155000 0.310  40000  680  500    4  25   41137]; RUNINFO{i}={130337, 'FNS test DAG'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 195500  195500 0.490  24600  128  500   20  25    4080]; RUNINFO{i}={103088, 'FNS DAG'}; i=i+1;
  RUN_RAWDATA(i,:)=[  38257   38257 0.231  28600  288  500   87  25   12591]; RUNINFO{i}={1447857,'LNS mb gmsh longestyet'}; i=i+1;
  RUN_RAWDATA(i,:)=[  11178    6766 0.485 160000   96  250   10  25     972]; RUNINFO{i}={144802, 'FNS tir de mine'}; i=i+1;
  RUN_RAWDATA(i,:)=[  12944    6597 0.385 160000   96  250   18  25    1136]; RUNINFO{i}={147954, 'FNS tir de mine 400Hz'}; i=i+1;
  RUN_RAWDATA(i,:)=[1226000 1184000 0.428 113000  384  500  131  50   32583]; RUNINFO{i}={147921, 'FNS mars insight'}; i=i+1;
  RUN_RAWDATA(i,:)=[   5120    3440 0.428 113000   48  500   11  50     575]; RUNINFO{i}={1484867,'FNS mars insight cut'}; i=i+1;
  RUN_RAWDATA(i,:)=[   5120    3440 0.428   8000   48  500   11  50      41]; RUNINFO{i}={1485893,'FNS mars insight cut'}; i=i+1;
  RUN_RAWDATA(i,:)=[   5120    3440 0.428   8000   48  500   15  50      41]; RUNINFO{i}={1486004,'FNS mars insight cut'}; i=i+1;
  RUN_RAWDATA(i,:)=[   5120    3440 0.428   8000   48  500    8  50      41]; RUNINFO{i}={1486113,'FNS mars insight cut'}; i=i+1;
  RUN_RAWDATA(i,:)=[  10318    5451 0.496  64000   96 1000   18  25     401]; RUNINFO{i}={1508049,'FNS tir de mine 40Hz'}; i=i+1;
  RUN_RAWDATA(i,:)=[1328000 1284800 0.391 393800 1008 2000  164 100   74988]; RUNINFO{i}={1529789,'FNS mars insight 20h 3hz w/  perioBC'}; i=i+1;
  RUN_RAWDATA(i,:)=[1660000 1606000 0.391 600000 1680 4000  200 100   51256]; RUNINFO{i}={1538139,'FNS mars insight 22h 3hz w/o perioBC'}; i=i+1;
  RUN_RAWDATA(i,:)=[  80518   80518 0.254  20600 1152 5000   87  50   86388]; RUNINFO{i}={1560350,'LNS mb gmsh 400km (3) refined and atmmodel reg'}; i=i+1;
  RUN_RAWDATA(i,:)=[  80518   80518 0.254  37200 1152 5000   87  50   86386]; RUNINFO{i}={1560541,'FNS mb gmsh 400km (3) refined and atmmodel reg'}; i=i+1;
  RUN_RAWDATA(i,:)=[1660000 1606000 0.469 490000 1680 5000  204 100   83018]; RUNINFO{i}={1633618,'FNS mars insight 20h 3hz w/  perioBC'}; i=i+1;
  RUN_RAWDATA(i,:)=[  23520   20160 0.039  46200  240  100  400 100     598]; RUNINFO{i}={150400, 'FNS mars insight 20h 0.1hz'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 204980  198320 0.044  70000  384 5000  857 100    8541]; RUNINFO{i}={150395, 'FNS mars insight incidence 20h 3hz'}; i=i+1;
  RUN_RAWDATA(i,:)=[  23520   20160 0.039 236600  240  100  400 100    3036]; RUNINFO{i}={151120, 'FNS mars insight f=0.1Hz but crashed'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 415000  401500 0.469 490000  640 5000  669 100   81734]; RUNINFO{i}={151319, 'FNS mars insight incidence 20h 3hz larger'}; i=i+1;
  col_dgpercent=1;
  col_snappercent=2;
  col_synthpercent=3;
  col_nbeltspproc=4;
  col_totalcputime_h=5;
  col_totalcputimepelpit_s=6;
  RUN_LV1DATA(:,col_dgpercent)      = RUN_RAWDATA(:,col_nbeltsdg)./RUN_RAWDATA(:,col_nbelts);
  RUN_LV1DATA(:,col_snappercent)    = RUN_RAWDATA(:,col_NstepOutImages)./RUN_RAWDATA(:,col_NSTEP);
  RUN_LV1DATA(:,col_synthpercent)   = RUN_RAWDATA(:,col_NstepOutSeismo)./RUN_RAWDATA(:,col_NSTEP);
  RUN_LV1DATA(:,col_nbeltspproc)    = RUN_RAWDATA(:,col_nbelts)./RUN_RAWDATA(:,col_NPROC);
  RUN_LV1DATA(:,col_totalcputime_h) = RUN_RAWDATA(:,col_realtime) .* RUN_RAWDATA(:,col_NPROC)/3600; % Total CPU time = real elapsed time * number of CPUS. In [h].
  RUN_LV1DATA(:,col_totalcputimepelpit_s) = RUN_LV1DATA(:,col_totalcputime_h) * 3600 ./ (RUN_RAWDATA(:,col_nbelts) .* RUN_RAWDATA(:,col_NSTEP)); % Total CPU time, per element, per iteration. In [s].
  
  x = [RUN_RAWDATA(:,col_SumNrec), RUN_LV1DATA(:,col_dgpercent), RUN_LV1DATA(:,col_snappercent), RUN_LV1DATA(:,col_synthpercent), RUN_LV1DATA(:,col_nbeltspproc)];
  t = RUN_LV1DATA(:, col_totalcputimepelpit_s);
end

function timestr=formatSecondsToHHMMSS(s)
  hms=fix(mod(s,[0,3600,60])./[3600,60,1]);
  timestr=strcat(sprintf('%8.f',hms(1))," h ",sprintf('%2.f',hms(2))," m ",sprintf('%2.f',hms(3)),' s');
end
