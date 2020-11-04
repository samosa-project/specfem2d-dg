% Author:        Léo Martire.
% Description:   Estimate new run time based on previous runs' saved data.
% Last modified: See file metadata.
% Usage:         Fill section "Parameters" according to the calculation to
%                be performed.
%                Use the scripts in the 'extract_run_information' folder
%                to fill data used for estimation (function 'load' below).
% Notes:         The method (N-D nearest point search) is approximate and
%                does not take into account the fact that some parameters
%                have more impact than others. All in all, it is a very
%                rough and very approximate method.

% clear all
% close all
clc
format longG;

[~] = setup_overall();

[data, t, info, raw_data] = load_data(); % Load data (see function below).

plot_rate = 1;
plotrate_selfulldgonly = 0; % 1 = exclude simulations with elastic elements
plotrate_percentdgthresh = 90; % select only simulations with more than plotrate_percentdgthresh % dg elements
plotrate_snappercentthresh = 100; % in [%]
plotrateversion = 'FNS';
% plotrateversion = 'LNS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display data.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=numel(t)-3:numel(t)
for i=1:numel(t)
  ID=info{i}(1); ID=ID{1};
  disp(['[',mfilename,'] Run ID ',sprintf('%7d',ID),' (',sprintf('%5.1f',data(i,2)*100),'% DG):       ',sprintf('%.2e',t(i)),' s CPU, per el., per it.       (''',char(string(info{i}(2))),''').']);
end
disp(' ');
Nshow=10; [~,idsort_t]=sort(t);
disp(['[',mfilename,'] Best ',num2str(Nshow),' number of elements per proc in terms of total CPU time per element and per iteration:']);
format shortg; disp([data(idsort_t(1:Nshow),end),t(idsort_t(1:Nshow))]'); format longG;
disp(['[',mfilename,'] First best being ''',char(info{idsort_t(1)}(2)),'''. Run ''figure();loglog(t,data(:,end),''.''); hold on; CVH=convhull(t,data(:,end)); loglog(t(CVH),data(CVH,end)); xlabel(''t''); ylabel(''elements per CPU'');'' to plot a convex hull.']);
% figure();loglog(t,data(:,end),'.'); hold on; CVH=convhull(t,data(:,end)); loglog(t(CVH),data(CVH,end)); xlabel('t'); ylabel('elements per CPU');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(not(plot_rate))
  % Choices. %%%%%%%%%%%%%%%%%%%%
  % nproc         = 8*48;
  % nbeltsElastic = 6600;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Auto-read from parfile. %%%%%
  parfile='';
  while(not(exist(parfile)==2))
    parfile = input(['[',mfilename,'] Path to simulation parfile > '],'s');
  end
  
  
  
  
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
    if(not(exist(tryMeshExtMeshLocalisation, 'file')))
      error('CANNOT FIND Mesh_extMesh IN ORDER TO GUESS THE NUMBER OF POINTS');
    else
      [~,r]=system(['head ',tryMeshExtMeshLocalisation,' -n 1']);
      neltot    = str2num(r);
      disp(['[',mfilename,']   GUESS for number of elements: ',num2str(neltot),'.']);
      disp(['[',mfilename,']   Check your ''Mesh_extMesh'' file. Press any key to continue.']);
      pause();
    end
    tryMaterialExtMeshLocalisation = [folderofparfile,'EXTMSH',filesep,'Material_extMesh'];
    if(not(exist(tryMaterialExtMeshLocalisation, 'file')))
      error('CANNOT FIND Material_extMesh IN ORDER TO GUESS NUMBER OF VISCOELASTIC ELEMENTS');
    else
      kek = importdata(tryMaterialExtMeshLocalisation);
      nbeltsElastic_GUESS = numel(kek)-sum(kek==1);
      disp(['[',mfilename,']   GUESS (assuming material n°1 is the acoustic material) for number of viscoelastic elements: ',num2str(nbeltsElastic_GUESS),'.']);
    end
  end
  nbeltsElastic = input(['[',mfilename,'] Number of viscoelastic elements? > ']);
  
  nstations = sum(readExampleFiles_extractParam(parfile, 'nrec', 'float'));
  nstepseismo = readExampleFiles_extractParam(parfile, 'NSTEP_BETWEEN_OUTPUT_SEISMOS', 'float');
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
  % Computation of estimate.    %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  weigth=[1,100,25,1,50]; % Importance of each parameter.
  meandata=mean(data,1);

  nproc = -1;
  while(nproc~=0)
    nproc = -1;
    while(not(numel(nproc)==1 & nproc>=0))
      disp(' ');disp(' ');disp(' ');
      nproc = input(['[',mfilename,'] Number of procs (>=0, 0 to stop script, for now your parfile has NPROC=',num2str(readExampleFiles_extractParam(parfile, 'NPROC', 'float')),')? > ']);
    end
    if(nproc==0)
      continue
    end
    point=[nstations,neldg/neltot,nstepsnap/nsteptot,nstepseismo/nsteptot,neltot/nproc]; % Format point as data format. Currently [% elements as DG, % timesteps as snapshots, elements per proc].

    disp(['[',mfilename,'] ',num2str(neltot),' elements including ',num2str(neldg),' DG elements. ',num2str(nsteptot),' time steps. ',num2str(nstations),' stations sampling every ',num2str(nstepseismo),' iterations. Snapshots taken every ',num2str(nstepsnap),' iterations. ',num2str(nproc),' CPUs.']);
    disp(" ");
    disp(['[',mfilename,']                     [        n_stations       percent_DG     percent_snap    percent_synth n_elems_per_proc]']);
    disp(strcat("[",mfilename,"] Current point:      [ ", sprintf("%17.5f", point), "]."));
    % idp=dsearchn(data,point);
    idp=dsearchn((data-meandata)./weigth,(point-meandata)./weigth); % Optimised search.
    disp(strcat("[",mfilename,"] Closest data point: [ ", sprintf("%17.5f", data(idp,:)), "] (",info{idp}{2},")."));
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
    sel = (data(:, idRaw_nbeltsdg)>plotrate_percentdgthresh/100); % Select full DG only.
    disp(['[',mfilename,'] Selected only simulations containing more than ',num2str(plotrate_percentdgthresh),'% of DG elements.']);
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
  for is=1:size(t); if(sel(is)); txt=info{is}(2); txt=txt{1}; if(not(isempty(regexp(txt,'mb gmsh')))); txt; disp('      removed mb gmsh run'); sel(is)=0; end; end; end
  for is=1:size(t); if(sel(is)); txt=info{is}(2); txt=txt{1}; if(not(isempty(regexp(txt,'mb huge')))); txt; disp('      removed mb huge run'); sel(is)=0; end; end; end
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
  y = t(sel); ylab = {['total CPU time per element per iteration [s]'], ['(real elapsed time $\times$ number of CPUS)']}; % Get real time (total CPU time p. el. p. it. / CPU).
  theoretical_pow = 0; theoretical_lab=['theoretical best: constant ($\propto n^{',num2str(theoretical_pow),'}$)'];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Choose color-code.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  c = log10(data(sel,idX_snappercent)); clab = ['log$_{10}$(snapshot frequency)']; % Color code with snapshot frequency.
  ux=sort(unique(x));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Choose size-code.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sizz = 50;
  sizz_raw = data(sel, idX_nbeltspproc);
  sizz=sizz_raw;
  sizz = (sizz-min(sizz));
  sizz = sizz/range(sizz);
  sizz = sizz.^0.25;
  sizz = sizz*300 + 5;
  slab = ['(nb. elts. per CPU)$^{0.25}$'];
%   plot(sort(sizz))
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Build best w.r.t. color.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  u_best_y_wrt_col=[];
  for ix=1:size(ux)
    sely = y(x == ux(ix));
    selc = c(x == ux(ix));
    u_best_y_wrt_col(ix, 1) = min(sely(selc == min(selc)));
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Build linear fits.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  beta_lfit_all = [ones(length(x),  1) log(x) ] \ log(y);
  beta_lfit_bestcol = [ones(length(ux), 1) log(ux)] \ log(u_best_y_wrt_col);
  lfitlab_all = ['linear fit best: $\propto n^{', sprintf('%.2g', beta_lfit_all(2)), '}$'];
  lfitlab_bst = ['linear fit best: $\propto n^{', sprintf('%.2g', beta_lfit_bestcol(2)), '}$'];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Actual plot.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  f = figure('units','normalized','outerposition', [0 0 1 1]);
  dname = ['data, color $\Leftrightarrow$ ', clab, ', size $\Leftrightarrow$ ', slab];
  dname = ['data, size $\propto$ ', slab];
  scatter(x, y, sizz, c, 'filled', 'displayname', dname); hold on;
  
  loglog(ux, exp(beta_lfit_all(1)) * ux.^beta_lfit_all(2), 'k-', 'displayname', ['linear fit all: $\propto n^{', sprintf('%.2g', beta_lfit_all(2)), '}$']);
  
  scatter(ux, u_best_y_wrt_col, 'displayname', ['best data for each abscissa w.r.t. color']);
  
  loglog(ux, exp(beta_lfit_bestcol(1)) * ux.^beta_lfit_bestcol(2), 'k--','displayname', lfitlab_bst);
  
  loglog(ux, ux.^theoretical_pow * (mean(y(ux == min(ux)))/min(ux)^theoretical_pow), 'k:','displayname', theoretical_lab);
  
  set(gca, 'xscale', 'log');
  set(gca, 'yscale', 'log');
  xlim([0.9*min(ux), 1.1*max(ux)]);
  xlabel(xlab);
  ylabel(ylab);
  
  h = colorbar;
  set(h, 'ticklabelinterpreter','latex');
  ylabel(h,clab,'interpreter','latex', 'fontsize',26);
  legend('location', 'best');
  
  title(['Parallelisation Efficiency - ', plotrateversion]);
  
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

function [x, t, RUNINFO, RUN_RAWDATA] = load_data()
  col_nbelts   = 1;
  col_nbeltsdg = 2;
  col_Cfl      = 3;
  col_NSTEP    = 4;
  col_NPROC    = 5;
  col_NstepOutImages=6;
  col_SumNrec  = 7;
  col_NstepOutSeismo=8;
  col_realtime = 9;
  [RUN_RAWDATA, RUNINFO] = load_raw_data();
  
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

function [RUN_RAWDATA, RUNINFO] = load_raw_data()
  [~, rootFolder] = setup_overall();
  dataFile = [rootFolder, filesep, 'utils', filesep, 'utils_specfem2d-dg', filesep, 'extract_run_information', filesep, 'runs_information.dat'];
  
  fid = fopen(dataFile);
  RAWDATA = textscan(fid, '%d %d %f %d %d %d %d %d %d %d %s', 'Headerlines', 1, 'Delimiter', newline);
  fclose('all');
  
  RUN_RAWDATA = zeros(size(RAWDATA{1}, 1), numel(RAWDATA)-2);
  for c=1:size(RUN_RAWDATA,2)
    RUN_RAWDATA(:, c) = RAWDATA{c};
  end
  RUNINFO = {};
  for i=1:size(RAWDATA{1}, 1)
    RUNINFO{i} = {RAWDATA{end-1}(i), RAWDATA{end}{i}};
  end
end