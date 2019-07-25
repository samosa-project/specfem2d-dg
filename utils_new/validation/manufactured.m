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

format compact;
set(0, 'DefaultLineLineWidth', 3); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 6); % Default at 6.
set(0, 'defaultTextFontSize', 28);
set(0, 'defaultAxesFontSize', 24); % Default at 10.
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

% Field plots?
plotFields_saveFigName_base = ['comparedDump__'];
plotFields_extsToSave={'jpg'};
plotFields_do = 1; % do or do not plot fields themselves
plotFields_xlab = '$x/L$';
plotFields_ylab = '$z/L$';
plotFields_NColorbarTicks = 9;

% Error plots.
figErr_saveFigName_base = ['MES_Error__'];
figErr_extsToSave = {'jpg'};
plotErr_title = '$L^2$ Error as Function of Number of Elements $N$';
plotErr_xlab = '$N$';
plotErr_figPos = [0,1e4,1600,600];

% IDs of dumps to plot.
IDz = 100;

% OUTPUT_FILES directory to analyse.

OFDs = {[SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_dx1p000']}; IDzs={[5000]};
% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_dx0p500']; IDz=[10000];
% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_dx0p250']; IDz=[20000];

% OFDs = {[SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_dx1p000'], ...
%         [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_dx0p500'], ...
%         [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_dx0p250']};
% % commonID = 4*[1:5]*1000;
% commonID = 4*[5]*1000;
% IDzs = {[commonID/4], ...
%         [commonID/2], ...
%         [commonID]};

% Parameters for analytic solution (cf. LNS_manufactured_solutions.mw).
RHO_cst = 0.001; dRHO_x = 1; dRHO_z = 2;
E_cst = 50*RHO_cst; dE_x=3; dE_z=4;

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
  thresh = 0.05; % relative to max
  CMAP = [0,0,0.25;
          0,0,0.5;
          0,0,1;
          colorOK; colorZero; colorOK;
          1,0,0;
          0.5,0,0;
          0.25,0,0];
  xcm=[-1, -0.85, -0.15, thresh]; xcm=[xcm, 0, -fliplr(xcm)];
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
    spl=split(OFD,'/'); plotFields_saveFigName_base = [plotFields_saveFigName_base, regexprep(spl{end-1},'OUTPUT_FILES_',''),'__'];% if sub OUTPUT_FILES folders per case
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
  if(readExampleFiles_extractParam([OFD,'input_parfile'],'USE_LNS','bool'))
    LNSorFNS = 'LNS';
  else
    LNSorFNS = 'FNS';
  end
  synthname = [LNSorFNS,' Synthetic'];
  %%%%%%%%%%%%%%%%
  
  % Save N.
  globalSave(ofdi,:,GS_ID_NX) = ones(1,numel(IDz))* NX; % save N in last dimension, for all iterations
  % Save DT.
  globalSave(ofdi,:,GS_ID_DT) = readExampleFiles_extractParam([OFD,'input_parfile'],'DT','float');
  
  %%%%%%%%%%%%%%%%
  % Loop dumps and plot agreement.
  %%%%%%%%%%%%%%%%
  for IT_id = 1:numel(IDz)
    IT = IDz(IT_id);
    % Load dump.
    [X,Y,V, imagetype_wavefield_dumps] = readDumps(OFD, IT);
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

    % Interpolate dump.
    delauTri = delaunayTriangulation(X, Y);
    tri = delauTri.ConnectivityList ;
    xi = delauTri.Points(:,1) ; yi = delauTri.Points(:,2) ;
    xi = (xi-min(xi)) / peak2peak(xi); % bring back x to [0, 1] (safety)
    yi = (yi-min(yi)) / peak2peak(yi); % bring back z to [0, 1] (safety)
    F = scatteredInterpolant(X,Y,V);
    zexp = F(xi,yi);
    caxxxxiiss_exp = [min(zexp), max(zexp)];

    % Build analytic solution (cf. LNS_manufactured_solutions.mw).
    % Constitutive variables.
    drho_th = RHO_cst*(sin(dRHO_x*pi*xi) + sin(dRHO_z*pi*yi));
    dvx_th = 0.*xi;
    dvz_th = 0.*xi;
    dE_th = E_cst * (sin(dE_x*pi*xi) + sin(dE_z*pi*yi));
    % Pressure.
    dp_th = (GAM-1)*(-0.5*V0_x^2*drho_th + dE_th); % Compute directly.
%     p0 = sound_velocity^2 * RHO0/GAM; E0 = p0/(GAM-1) + 0.5*RHO0*V0_x^2; dp_th = (GAM-1)*(E0+dE_th - 0.5*(RHO0+drho_th).*((V0_x+dvx_th).^2 + dvz_th.^2)) - p0; % Deduced from (rho, v, E)
    
    switch(imagetype_wavefield_dumps)
      case 4
        zth = dp_th;
        qtity = 'p';
      case 2
        zth = dvx_th;
        qtity = 'v_x';
      otherwise
        error(['[] imagetype_wavefield_dumps not implemented']);
    end
    caxxxxiiss_th = max(abs(zexp))*[-1,1];
    errorPrefix = [errorPrefix_base, qtity];
    
    % Compute error.
    zerr = zexp-zth;
  %   errName = ['(', synthname, ' - Analytic)'];
    errName = [errorPrefix,',',num2str(NX),'}\left(',num2str(IT),'\right)'];
    caxxxxiiss_err = [min(zerr), max(zerr)];
    L2NormOfError = (integrate2DDelaunayTriangulation(delauTri, zerr.^2))^0.5;
    
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
        glob_caxx = max(abs(zexp))*[-1,1];
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
      feg = figure('units','normalized','outerposition',[0 0.4 1 0.6]);
      axxx = tight_subplot(1, 3, 0.012, [0.14,0.08], [0.05, 0.10]); % not mandatory, but prettier
      axes(axxx(1));
      % axxx(1)=subplot(1,3,1);
      % surf(Xmg,Ymg,Vmg); shading interp; view([0,0,1]); colormap jet;
      trisurf(tri,xi,yi,zth); shading interp; view([0,0,1]); colormap(CMAP); caxis(glob_caxx); %colorbar;
      xlabel(plotFields_xlab); ylabel(plotFields_ylab);
      title('Analytic');
      % Synthetic. %%%
      axes(axxx(2));
      % axxx(2)=subplot(1,3,2);
      trisurf(tri,xi,yi,zexp); shading interp; view([0,0,1]); colormap(CMAP); caxis(glob_caxx); %colorbar;
      yticklabels([]);
      xlabel(plotFields_xlab);
      title([synthname,' ($n=',num2str(IT),'$)']);
      % Error. %%%%%%%
      axes(axxx(3));
      % axxx(3)=subplot(1,3,3);
      trisurf(tri,xi,yi,zerr); shading interp; view([0,0,1]); colormap(CMAP); caxis(glob_caxx);
      xlabel(plotFields_xlab);
      yticklabels([]);
      title({[errName,'=',scientific_latex_notation(L2NormOfError,2),'$']});
      % Color bar. %%%
      posAxxx3 = get(axxx(3),'Position');
      h_cb = colorbar('Position', [posAxxx3(1)+posAxxx3(3)+0.03  posAxxx3(2)  0.03  posAxxx3(4)]);
      colorbarticks=linspace(-1,1,plotFields_NColorbarTicks); colorbarticks=sort([colorbarticks, [-1,1]*max(thresh)]);
      if(max(abs(zth))>0)
        colorbarticks = colorbarticks*max(abs(zth));
      else
        colorbarticks = colorbarticks*max(abs(zexp));
      end
      colorbarticks = colorbarticks;
      set(h_cb, 'ytick', colorbarticks);
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
  figErr_saveFigFullPath = [savefigpath, figErr_saveFigName_base, datestr(now,'YYmmDD_HHMMss')];
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

  % Plot last time as standalone.
  fig_error = figure('units','pixels','position',plotErr_figPos);
  tight_subplot(1, 1, -1, [0.2,0.11], [0.14, 0.04]);
  semilogy(N, EPS,'displayname',['@$t=',num2str(t),'$~s']);
  xlabel(plotErr_xlab); ylabel(plotErr_ylab); title(plotErr_title);
  xlim(plotErr_xlim);
  tickLabels_general_reskin_y(fig_error);
  prettyAxes(fig_error);
  customSaveFig(fig_error, figErr_saveFigFullPath, figErr_extsToSave);
end
%%%%%%%%%%%%%%%%