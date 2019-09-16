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

format shortg;
set(0, 'DefaultLineLineWidth', 3); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 6); % Default at 6.
set(0, 'defaultTextFontSize', 24);
set(0, 'defaultAxesFontSize', 22); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools');
SPCFM_EX_DIR = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefix = 'validation_lns_manufactured';
savefigpath = [SPCFM_EX_DIR,prefix,'_info/'];
errorPrefix_base = ['$\varepsilon_{'];

% N for interpolation grid.
Nexact = 1000; % Choose Nexact >> the N's we want to test. 10000 chugs HARD, 1000 runs ok. Up to N=100 points on the simulation's mesh, Nexact=1000 and Nexact=5000 yield the exact same L2 error.

% Field plots?
plotFields_do = 1; % do or do not plot fields themselves
plotFields_saveFigName_prefix = ['comparedDump__'];
plotFields_extsToSave={'jpg'};
plotFields_xlab = '$x/L$';
plotFields_ylab = '$z/L$';
plotFields_NColorbarTicks = 9;
plotFields_addErrorPanel = 1; % Set to 1 for debug, set to 0 for printing nice plots
plotFields_CMAP_thresh = 0.01; % relative to max

% Error plots.
figErr_saveFigName_base = ['MES_Error__'];
figErr_extsToSave = {'jpg'};
plotErr_title = '$L^2$ Error as Function of Number of Elements $N$';
plotErr_xlab = '$N$';
plotErr_figPos = [0,1e4,1600,600];

% IDs of dumps to plot.
IDz = 100;

% OUTPUT_FILES directory to analyse.
% OFDs = {[SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_dx1p000_mu_dp__1p25em5']}; IDzs={[2]*10000}; plotFields_do = 0;
OFDs = {[SPCFM_EX_DIR,prefix,'_1p00_iv/OUTPUT_FILES']}; IDzs={[100e3]}; plotFields_do = 1;

% Test case to plot/compute.
testCase = 'inviscid';
% testCase = 'kappa';
% testCase = 'mu';

%%%%%%%%%%%%%%%%
% Packed, for paper.
%%%%%%%%%%%%%%%%
% OFDs = {[SPCFM_EX_DIR,prefix,'_5p00_iv/OUTPUT_FILES'], ...
%         [SPCFM_EX_DIR,prefix,'_1p00_iv/OUTPUT_FILES'], ...
%         [SPCFM_EX_DIR,prefix,'_0p50_iv/OUTPUT_FILES'], ...
%         [SPCFM_EX_DIR,prefix,'_0p25_iv/OUTPUT_FILES']}; testCase = 'inviscid';
% OFDs = {[SPCFM_EX_DIR,prefix,'_5p00_ka/OUTPUT_FILES'], ...
%         [SPCFM_EX_DIR,prefix,'_1p00_ka/OUTPUT_FILES'], ...
%         [SPCFM_EX_DIR,prefix,'_0p50_ka/OUTPUT_FILES'], ...
%         [SPCFM_EX_DIR,prefix,'_0p25_ka/OUTPUT_FILES']}; testCase = 'kappa';
OFDs = {[SPCFM_EX_DIR,prefix,'_5p00_mu/OUTPUT_FILES'], ...
        [SPCFM_EX_DIR,prefix,'_1p00_mu/OUTPUT_FILES'], ...
        [SPCFM_EX_DIR,prefix,'_0p50_mu/OUTPUT_FILES'], ...
        [SPCFM_EX_DIR,prefix,'_0p25_mu/OUTPUT_FILES']}; testCase = 'mu';
% commonID = [0.2,0.4,0.6,0.8,1]*1000;
commonID = [1]*20000;
IDzs = {[commonID], ...
        [commonID*5], ...
        [commonID*5*2], ...
        [commonID*5*2*2]};
plotFields_do = 0; % eventually deactivate plotting

% Parameters for analytic solution (cf. LNS_manufactured_solutions.mw).
% Should be the same as in source code.
switch(testCase)
  case 'inviscid'
    colourPlot = [1,0,0]; markerPlot = 'o';
    RHO_cst = 0.001;
    VX_cst = 0.;
    E_cst = 0.05;
  case 'kappa'
    colourPlot = [0,1,0]; markerPlot = 'diamond';
    RHO_cst = 0.;
    VX_cst = 0.;
    E_cst = 0.05;
  case 'mu'
    colourPlot = [0,0,1]; markerPlot = 'square';
    RHO_cst = 0.;
    VX_cst = 0.001;
    E_cst = 0.;
  otherwise
    error(['[, ERROR]']);
end
dRHO_x = 1;
dRHO_z = 2;
dVX_x = 2;
dVX_z = 0;
VZ_cst = 0.;
dVZ_x = 1;
dVZ_z = 2;
dE_x = 3;
dE_z = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(not(numel(OFDs)==numel(IDzs)))
  error(['[, ERROR] must provide same number of OFD as IDz']);
else
  N_OFD = numel(OFDs);
end

if(plotFields_do)
  % Prepare colorbar.
  colorOK = [1,1,1]; % colour for |value|<thresh
  colorZero = [255,228,181]/255; % moccasin
  % CMAP = 'jet'; % default colorbar
  CMAP = [0,0,0.25;
          0,0,0.5;
          0,0,1;
          colorOK; colorZero; colorOK;
          1,0,0;
          0.5,0,0;
          0.25,0,0];
  xcm=[-1, -0.85, -0.15, plotFields_CMAP_thresh]; xcm=[xcm, 0, -fliplr(xcm)];
  N=1000; % number of points for color values
  CMAP = interp1(xcm, CMAP, linspace(min(xcm),max(xcm),N)');
  colormap(CMAP);
end

globalSave = []; % prepare saving of errors w.r.t. iterations for each OFD
% Storage IDs.
GS_ID_NX = 1;
GS_ID_DT = 2;
GS_ID_IT = 3;
GS_ID_EPS = 4;
for ofdi = 1:N_OFD
  OFD = OFDs{ofdi};
  IDz = IDzs{ofdi};
  if(not(OFD(end)=='/')); OFD = [OFD,'/']; end; % Safeguard.
  disp(['[] Looking at ''',OFD,'''.']);
  
  if(plotFields_do)
    % Prepare Figure saving w.r.t. OFD.
    spl=split(OFD,'/'); plotFields_saveFigName_base = [plotFields_saveFigName_prefix, regexprep(spl{end-2},[prefix,'_'],''),'__',regexprep(spl{end-1},'OUTPUT_FILES_',''),'__'];% if sub OUTPUT_FILES folders per case
  end

  %%%%%%%%%%%%%%%%
  % Load parameters from simulation.
  %%%%%%%%%%%%%%%%
  MODEL = readExampleFiles_extractParam([OFD,'input_parfile'],'MODEL','string');
  USE_ISOTHERMAL_MODEL = readExampleFiles_extractParam([OFD,'input_parfile'],'USE_ISOTHERMAL_MODEL','bool');
  if(not(strcmp(MODEL,'default') & not(USE_ISOTHERMAL_MODEL)))
    error(['fluid model is wrong, set MODEL=default and USE_ISOTHERMAL_MODEL=.false. ']);
  end
  NX = readExampleFiles_extractParam([OFD,'input_parfile'],'nx','int');
  RHO0 = readExampleFiles_extractParam([OFD,'input_parfile'],'surface_density','float');
  sound_velocity = readExampleFiles_extractParam([OFD,'input_parfile'],'sound_velocity','float');
  V0_x = readExampleFiles_extractParam([OFD,'input_parfile'],'wind','float');
  GAM = readExampleFiles_extractParam([OFD,'input_parfile'],'constant_p','float') / readExampleFiles_extractParam([OFD,'input_parfile'],'constant_v','float');
  P0 = sound_velocity^2*RHO0/GAM; % deduce p0 from isobaric hypothesis
  E0 = P0/(GAM-1) + 0.5*RHO0*V0_x^2; % deduce E0 from isobaric hypothesis
  if(readExampleFiles_extractParam([OFD,'input_parfile'],'USE_LNS','bool'))
    LNSorFNS = 'LNS';
  else
    LNSorFNS = 'FNS';
  end
  synthname = [LNSorFNS,' Synthetic'];
  KAPPA = readExampleFiles_extractParam([OFD,'input_parfile'],'thermal_conductivity','float');
  MU = readExampleFiles_extractParam([OFD,'input_parfile'],'dynamic_viscosity','float');
  if(max(MU,KAPPA)>0)
    VISCOUS = 1;
    disp(['[] VISCOUS']);
    if(strcmp(testCase,'inviscid'))
      error(['[, ERROR] Cannot test inviscid with viscosity activated in this OUTPUT_FILES.'])
    end
    if(RHO_cst~=0)
      error(['[ERROR] for viscous tests, RHO_cst must be zero, and now is not']);
    end
    if(strcmp(testCase,'kappa') && not(VX_cst==0 & VZ_cst==0))
      error(['[ERROR] for viscous kappa test, VX_cst and VZ_cst must be zero, and now are not']);
    end
    if(strcmp(testCase,'kappa') && E_cst==0)
      error(['[ERROR] for viscous kappa test, E_cst must be non zero, and now is zero']);
    end
    if(strcmp(testCase,'mu') && not(E_cst==0))
      error(['[ERROR] for viscous mu test, E_cst must be zero, and now is not']);
    end
    if(strcmp(testCase,'mu') && VX_cst==0)
      error(['[ERROR] for viscous mu test, VX_cst must be non zero, and now is zero']);
    end
  else
    VISCOUS = 0;
    disp(['[] INVISCID']);
    if(strcmp(testCase,'kappa') || strcmp(testCase,'mu'))
%       error(['[, ERROR] Cannot test viscous with viscosity deactivated in this OUTPUT_FILES.'])
    end
    if(not(VX_cst==0 & VZ_cst==0))
%       error(['[ERROR] for inviscid test, VX_cst and VZ_cst must be zero, and now are not']);
    end
    if(RHO_cst==0)
      error(['[ERROR] for inviscid test, RHO_cst must be non zero, and now is zero']);
    end
  end
  %%%%%%%%%%%%%%%%
  
  % Save N.
  globalSave(ofdi,:,GS_ID_NX) = ones(1,numel(IDz))* NX; % save N in last dimension, for all iterations
  % Save DT.
  DT = readExampleFiles_extractParam([OFD,'input_parfile'],'DT','float');
  globalSave(ofdi,:,GS_ID_DT) = DT;
  
  %%%%%%%%%%%%%%%%
  % Loop dumps and plot agreement.
  %%%%%%%%%%%%%%%%
  for IT_id = 1:numel(IDz)
    IT = IDz(IT_id);
    % Load dump.
    [X,Y,V, imagetype_wavefield_dumps] = readDumps(OFD, IT, 0);
    disp(['[] Read dump at iteration ',num2str(IT),' (t = ',num2str(IT*DT),' s).']);
    if(ismember(imagetype_wavefield_dumps,[1,2,3]))
      disp(['[] Selecting x-axis values from dump.']);
      V = V(:,1);
    end
  %   if(ismember(imagetype_wavefield_dumps,[4]))
  %     disp(['[] Dump was dp/rho0, renormalise it to dp.']);
  %     V = V * RHO0;
  %   end
    
    if(plotFields_do)
      % Prepare output figure name based on dump ID.
      savefigname = [plotFields_saveFigName_base, 'IT',sprintf('%012.0f',IT)];
      savefigfullpath = [savefigpath, regexprep(savefigname,'\.','')]; % regexprep because latex crashed when filenames have dots in it
    end
    
    % Create independent mesh onto which interpolate.
    x_exact = linspace(min(X),max(X),Nexact);
    y_exact = linspace(min(Y),max(Y),Nexact);
    [Xe, Ye] = meshgrid(x_exact,y_exact);
    
    % Interpolate dump.
%     delauTri = delaunayTriangulation(X, Y);
%     tri = delauTri.ConnectivityList;
%     xi = delauTri.Points(:,1); yi = delauTri.Points(:,2);
%     xi = (xi-min(xi)) / peak2peak(xi); % bring back x to [0, 1] (safety)
%     yi = (yi-min(yi)) / peak2peak(yi); % bring back z to [0, 1] (safety)
    F = scatteredInterpolant(X,Y,V);
%     zexp = F(xi,yi);
    zexp = F(Xe, Ye);
    caxxxxiiss_exp = [min(min(zexp)), max(max(zexp))];

    % Build analytic solution (cf. LNS_manufactured_solutions.mw).
    % Constitutive variables.
    drho_th = RHO_cst * (sin(dRHO_x*pi*Xe) + sin(dRHO_z*pi*Ye));
    dvx_th  = VX_cst  * (sin(dVX_x*pi*Xe)  + sin(dVX_z*pi*Ye));
    dvz_th  = VZ_cst  * (sin(dVZ_x*pi*Xe)  + sin(dVZ_z*pi*Ye));
    dE_th   = E_cst   * (sin(dE_x*pi*Xe)   + sin(dE_z*pi*Ye));
    % Pressure.
    %dp_th = (GAM-1)*(-0.5*V0_x^2*drho_th + dE_th); % Compute directly.
    dp_th = (GAM-1)*(E0+dE_th - 0.5*(RHO0+drho_th).*((V0_x+dvx_th).^2 + dvz_th.^2)) - P0;
    
    switch(imagetype_wavefield_dumps)
      case 4
        disp(['[] Dumps are pressure dumps, computing and plotting error w.r.t. pressure.']);
        zth = dp_th;
        qtity = 'p';
      case 2
        disp(['[] Dumps are horizontal velocity dumps, computing and plotting error w.r.t. horizontal velocity.']);
        zth = dvx_th;
        qtity = 'v_x';
      otherwise
        error(['[] imagetype_wavefield_dumps not implemented']);
    end
    caxxxxiiss_th = max(max(abs(zexp)))*[-1,1];
    errorPrefix = [errorPrefix_base, qtity];
    disp(['[] Integral of zth = ',num2str(trapz(x_exact,trapz(y_exact,zth)))]);
    
    % Compute error.
    zerr = zexp-zth;
  %   errName = ['(', synthname, ' - Analytic)'];
    errName = [errorPrefix,'}\left(',num2str(NX),'\right)'];
    caxxxxiiss_err = [min(min(zerr)), max(max(zerr))];
    %L2NormOfError = (integrate2DDelaunayTriangulation(delauTri, zerr.^2))^0.5;
    L2NormOfError = (trapz(x_exact,trapz(y_exact,zerr.^2)))^0.5;
    disp(['[] L_{inf} norm of error = ',num2str(max(max(zerr)))]);
    
    globalSave(ofdi,IT_id,GS_ID_IT) = IT; % save iteration in last dimension
    globalSave(ofdi,IT_id,GS_ID_EPS) = L2NormOfError; % save error in last dimension
    
    %%%%%%%%%%%%%%%%
    % Plot fields.
    %%%%%%%%%%%%%%%%
    if(plotFields_do)
      % Colorbar scale
      all_caxx = [caxxxxiiss_th;caxxxxiiss_exp;caxxxxiiss_err];
      glob_caxx = caxxxxiiss_th;
      if(numel(unique(glob_caxx))==1)
        % If zth = cste, choose glob_caxx based on zexp
        glob_caxx = max(max(abs(zexp)))*[-1,1];
        % If again zexp = cste, choose glob_caxx arbitrarily
        if(numel(unique(glob_caxx))==1)
          if(glob_caxx(1)==0)
            % If cste=0, choose arbitrarily a maximum and minimum
            glob_caxx = [-1,1]*0.1;
          else
            glob_caxx = glob_caxx(1)*[0.9, 1.1];
          end
        end
      end
      % Figure. %%%%%%
      if(plotFields_addErrorPanel)
        nbSubplots = 3;
      else
        nbSubplots = 2;
      end
      feg = figure('units','normalized','outerposition',[0 0.4 1 0.6]);
      axxx = tight_subplot(1, nbSubplots, 0.012, [0.14,0.08], [0.05, 0.12]); % not mandatory, but prettier
      axes(axxx(1));
      % axxx(1)=subplot(1,3,1);
%       trisurf(tri,xi,yi,zth);
      surf(Xe,Ye,zth);
      shading interp; view([0,0,1]); colormap(CMAP); caxis(glob_caxx); %colorbar;
      xlabel(plotFields_xlab); ylabel(plotFields_ylab);
      title('Analytic');
      % Synthetic. %%%
      axes(axxx(2));
      % axxx(2)=subplot(1,3,2);
%       trisurf(tri,xi,yi,zexp);
      surf(Xe,Ye,zexp);
      shading interp; view([0,0,1]); colormap(CMAP); caxis(glob_caxx); %colorbar;
      yticklabels([]);
      xlabel(plotFields_xlab);
%       title([synthname,' ($n=',num2str(IT),'$)']);
      title([synthname,' ($t=',scientific_latex_notation(DT*IT,2),'$~s)']);
      if(plotFields_addErrorPanel)
      % Error. %%%%%%%
        axes(axxx(3));
        % axxx(3)=subplot(1,3,3);
  %       trisurf(tri,xi,yi,zerr);
        surf(Xe,Ye,zerr);
        shading interp; view([0,0,1]); colormap(CMAP); caxis(glob_caxx);
        xlabel(plotFields_xlab);
        yticklabels([]);
        title({[errName,'=',scientific_latex_notation(L2NormOfError,2),'$']});
      end
      % Color bar. %%%
      posAxxx3 = get(axxx(nbSubplots),'Position');
      h_cb = colorbar('Position', [posAxxx3(1)+posAxxx3(3)+0.03  posAxxx3(2)  0.03  posAxxx3(4)]);
      colorbarticks=linspace(-1,1,plotFields_NColorbarTicks); %colorbarticks=sort([colorbarticks, [-1,1]*max(thresh)]);
      if(max(max(abs(zth)))>0)
        colorbarticks = colorbarticks*max(max(abs(zth)));
      else
        colorbarticks = colorbarticks*max(max(abs(zexp)));
      end
      colorbarticks = colorbarticks;
      set(h_cb, 'ytick', colorbarticks);
      set(h_cb,'tickdir','both');
      % Link = linkprop(axxx,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
      Link = linkprop(axxx,{'CameraPosition', 'XLim', 'YLim'});
      setappdata(gcf, 'StoreTheLink', Link);
      prettyAxes(feg);
      customSaveFig(feg, savefigfullpath,plotFields_extsToSave);
    end
    %%%%%%%%%%%%%%%%
  end
  
  % Ended loop on iterations, save this OFD before going to the next.
%   globalSave(ofdi,:,:) = [save_IT;save_EPS];
end
%%%%%%%%%%%%%%%%

squeeze(globalSave)

%%%%%%%%%%%%%%%%
% Error plots.
%%%%%%%%%%%%%%%%
if(size(globalSave,1)>1)
  % If more than one OFD.
  plotErr_ylab = [errorPrefix,'}(N)$'];
  plotErr_xlim = [min(min(globalSave(:,:,GS_ID_NX))),max(max(globalSave(:,:,GS_ID_NX)))];

  % Plot progression w.r.t. time.
  ids_IT_to_plot = 1:size(globalSave,2);
  colours = winter(size(globalSave,2));
  fig_error_progress = figure('units','pixels','position',plotErr_figPos);
  tight_subplot(1, 1, -1, [0.2, 0.11], [0.1, 0.04]);
  for id_IT_to_plot = ids_IT_to_plot
    IT_1p000dx = globalSave(1,id_IT_to_plot,2);

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

  % Plot last time as standalone.
  fig_error = figure('units','pixels','position',plotErr_figPos);
  tight_subplot(1, 1, -1, [0.2,0.11], [0.14, 0.04]);
  loglog(N, EPS,'displayname',testCase,'linewidth',1.5,'marker',markerPlot,'markersize',10,'color',colourPlot,'markerfacecolor',colourPlot);
  legend('location','northeast');
  xlabel(plotErr_xlab); ylabel(plotErr_ylab); title(plotErr_title);
  xlim(plotErr_xlim);
  tickLabels_general_reskin_y(fig_error);
  prettyAxes(fig_error);
  figErr_saveFigFullPath = [savefigpath, figErr_saveFigName_base, datestr(now,'YYmmDD_HHMMss')];
  customSaveFig(fig_error, figErr_saveFigFullPath, figErr_extsToSave);
end
%%%%%%%%%%%%%%%%