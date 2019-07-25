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

format compact;
set(0, 'DefaultLineLineWidth', 2); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 6); % Default at 6.
set(0, 'defaultTextFontSize', 18);
set(0, 'defaultAxesFontSize', 16); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools');
SPCFM_EX_DIR = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefix = 'validation_lns_manufactured';
savefigname_base = ['comparedDump__'];

% IDs of dumps to plot.
IDz = 100;

% OUTPUT_FILES directory to analyse.
% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_dp'];
% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_dp_190603'];
OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES']; IDz=[500,[1,2]*1000];
% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_exact']; IDz=[100];
% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_force_exact_flux'];
% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_jumpexactlambda'];
% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_test_dump_dvx'];
% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_test_dump_dvz'];

% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_baseline'];
% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_try_gmsh'];
% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_force_exact_flux_inside'];
% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_force_exact_flux_inside+source_at_center'];
% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_z-p5']; IDz=[1,2,3,4,5,6]*100;
% OFD = [SPCFM_EX_DIR,prefix,'/OUTPUT_FILES_x-p5']; IDz=[1,2,3]*1000;

% Parameters for analytic solution (cf. LNS_manufactured_solutions.mw).
RHO_cst = 0.001; dRHO_x = 1; dRHO_z = 2;
E_cst = 50*RHO_cst; dE_x=3; dE_z=4;

% Plots.
extsToSave={'jpg'};
% extsToSave={'jpg','eps'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(not(OFD(end)=='/')); OFD = [OFD,'/']; end; % Safeguard.
savefigpath = [SPCFM_EX_DIR,prefix,'_info/'];

spl=split(OFD,'/');
% savefigname = [savefigname, regexprep(spl{end-2},prefix,'')];% if one EXAMPLE per case
savefigname_base = [savefigname_base, regexprep(spl{end-1},'OUTPUT_FILES_',''),'__'];% if sub OUTPUT_FILES folders per case

%%%%%%%%%%%%%%%%
% Prepare colorbar.
%%%%%%%%%%%%%%%%
Ncolorbarticks = 9;
colorOK = [1,1,1]; % colour for |value|<thresh
colorZero = [255,228,181]/255; % moccasin
% CMAP = 'jet'; % default colorbar
thresh = 0.05; % relative to max
CMAP = [0,0,0.25;
        0,0,0.5;
        0,0,1;
        colorOK;
        colorZero;
        colorOK;
        1,0,0;
        0.5,0,0;
        0.25,0,0];
xcm=[-1, -0.85, -0.15, thresh]; xcm=[xcm, 0, -fliplr(xcm)];
N=1000; % number of points for color values
CMAP = interp1(xcm, CMAP, linspace(min(xcm),max(xcm),N)');
colormap(CMAP);
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
% Load parameters from simulation.
%%%%%%%%%%%%%%%%
MODEL = readExampleFiles_extractParam([OFD,'input_parfile'],'MODEL','string');
USE_ISOTHERMAL_MODEL = readExampleFiles_extractParam([OFD,'input_parfile'],'USE_ISOTHERMAL_MODEL','bool');
if(not(strcmp(MODEL,'default') & not(USE_ISOTHERMAL_MODEL)))
  error(['fluid model is wrong, set MODEL=default and USE_ISOTHERMAL_MODEL=.false. ']);
end
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

%%%%%%%%%%%%%%%%
% Loop dumps and plot agreement.
%%%%%%%%%%%%%%%%
for IT = IDz
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
  
  % Prepare output figure name based on dump ID.
  savefigname = [savefigname_base, 'IT',sprintf('%012.0f',IT)];
  savefigfullpath = [savefigpath, regexprep(savefigname,'\.','')]; % regexprep because latex crashed when filenames have dots in it

  % Interpolate dump.
  delauTri = delaunayTriangulation(X, Y);
  tri = delauTri.ConnectivityList ;
  xi = delauTri.Points(:,1) ; yi = delauTri.Points(:,2) ;
  xi = (xi-min(xi)) / peak2peak(xi); % bring back x to [0, 1] (safety)
  yi = (yi-min(yi)) / peak2peak(yi); % bring back z to [0, 1] (safety)
  F = scatteredInterpolant(X,Y,V);
  zexp = F(xi,yi);
  caxxxxiiss_exp = [min(zexp),max(zexp)];

  % Build analytic solution (cf. LNS_manufactured_solutions.mw).
  % Constitutive variables.
  drho_th = RHO_cst*(sin(dRHO_x*pi*xi) + sin(dRHO_z*pi*yi));
  dvx_th = 0.*xi;
  dvz_th = 0.*xi;
  dE_th = E_cst * (sin(dE_x*pi*xi) + sin(dE_z*pi*yi));
  % Pressure.
  dp_th = (GAM-1)*(-0.5*V0_x^2*drho_th + dE_th); % Compute directly.
%   p0 = sound_velocity^2 * RHO0/GAM; E0 = p0/(GAM-1) + 0.5*RHO0*V0_x^2; dp_th = (GAM-1)*(E0+dE_th - 0.5*(RHO0+drho_th).*((V0_x+dvx_th).^2 + dvz_th.^2)) - p0; % Deduced from (rho, v, E)
  
  switch(imagetype_wavefield_dumps)
    case 4
      zth = dp_th;
    case 2
      zth = dvx_th;
    otherwise
      error(['[] imagetype_wavefield_dumps not implemented']);
  end
  caxxxxiiss_th = max(abs(zexp))*[-1,1];

  % error
  zerr = zexp-zth; errName = ['(', synthname, ' - Analytic)'];
  caxxxxiiss_err = [min(zerr),max(zerr)];
  L2NormOfError = (integrate2DDelaunayTriangulation(delauTri, zerr.^2))^0.5;

  % colorbar scale
  all_caxx = [caxxxxiiss_th;caxxxxiiss_exp;caxxxxiiss_err];
%   glob_caxx = [min(all_caxx(:,1)),max(all_caxx(:,2))];
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

  %plot
  xlab = '$x/L$';
  ylab = '$z/L$';

  feg = figure('units','normalized','outerposition',[0 0.4 1 0.6]);
  axxx = tight_subplot(1, 3, 0.012, [0.14,0.08], [0.05, 0.10]); % not mandatory, but prettier
  
  axes(axxx(1));
  % axxx(1)=subplot(1,3,1);
  % surf(Xmg,Ymg,Vmg); shading interp; view([0,0,1]); colormap jet;
  trisurf(tri,xi,yi,zth); shading interp; view([0,0,1]); colormap(CMAP); caxis(glob_caxx); %colorbar;
  xlabel(xlab);
  ylabel(ylab);
  title('Analytic');

  axes(axxx(2));
  % axxx(2)=subplot(1,3,2);
  trisurf(tri,xi,yi,zexp); shading interp; view([0,0,1]); colormap(CMAP); caxis(glob_caxx); %colorbar;
  yticklabels([]);
  xlabel(xlab);
  title([synthname,' ($n=',num2str(IT),'$)']);

  axes(axxx(3));
  % axxx(3)=subplot(1,3,3);
  trisurf(tri,xi,yi,zerr); shading interp; view([0,0,1]); colormap(CMAP); caxis(glob_caxx);
  
  % Color bar.
  posAxxx3 = get(axxx(3),'Position');
  h_cb = colorbar('Position', [posAxxx3(1)+posAxxx3(3)+0.03  posAxxx3(2)  0.03  posAxxx3(4)]);
%   h_cb = colorbar;
  colorbarticks=linspace(-1,1,Ncolorbarticks); colorbarticks=sort([colorbarticks, [-1,1]*max(thresh)]);
  if(max(abs(zth))>0)
    colorbarticks = colorbarticks*max(abs(zth));
  else
    colorbarticks = colorbarticks*max(abs(zexp));
  end
  colorbarticks = colorbarticks;
  set(h_cb, 'ytick', colorbarticks);
  
  xlabel(xlab);
  yticklabels([]);
  title({[errName,', $L^2$ norm = $',scientific_latex_notation(L2NormOfError,2),'$']});
  
  % Link = linkprop(axxx,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
  Link = linkprop(axxx,{'CameraPosition', 'XLim', 'YLim'});
  setappdata(gcf, 'StoreTheLink', Link);

  prettyAxes(feg);
  customSaveFig(feg, savefigfullpath,extsToSave);
end
%%%%%%%%%%%%%%%%