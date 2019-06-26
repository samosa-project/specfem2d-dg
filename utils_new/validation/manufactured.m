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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefix = 'validation_lns_manufactured';
% output files
% OFD = ['/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/',prefix,'/OUTPUT_FILES_dp'];
% OFD = ['/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/',prefix,'/OUTPUT_FILES_dp_190603'];
% OFD = ['/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/',prefix,'/OUTPUT_FILES'];
% OFD = ['/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/',prefix,'/OUTPUT_FILES_force_exact_flux'];
OFD = ['/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/',prefix,'/OUTPUT_FILES_force_exact_flux_inside'];

if(not(OFD(end)=='/')); OFD = [OFD,'/']; end; % Safeguard.

savefigname = ['comparedDump__'];
savefigpath = ['/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/',prefix,'_info/'];

spl=split(OFD,'/');
% savefigname = [savefigname, regexprep(spl{end-2},prefix,'')];% if one EXAMPLE per case
savefigname = [savefigname, regexprep(spl{end-1},'OUTPUT_FILES_',''),'__'];% if sub OUTPUT_FILES folders per case

CMAP = 'jet';

CMAP = [0,0,1;
        0,153/255,153/255;
        1,1,1;
        1,1,1;
        1,0,0;
        188/255,0,1];
thresh=0.1;
xcm=[-2,-1,thresh]; xcm=[xcm,-fliplr(xcm)]; N=1000;
CMAP = interp1(xcm, CMAP, linspace(min(xcm),max(xcm),N)');

% for IT = (1300:100:2000)
% for IT = 5
for IT = 1000
% for IT = 2000
% for IT = 3000
% for IT = [4:6]*1000
% for IT = 5000
% for IT = 7000
% for IT = 10000
% for IT = 15000
% for IT = 20000
% for IT = [5,10,15,20]*1000
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Treatment.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % safety.
  OFD = char(OFD);
  if(not(strcmp(OFD(end),filesep)))
    OFD = [OFD,filesep];
  end

  % load dumps
  imagetype_wavefield_dumps = readExampleFiles_extractParam([OFD,'input_parfile'],'imagetype_wavefield_dumps','int');
  switch(imagetype_wavefield_dumps)
    case(4)
      disp(['reading pressure dumps']);
    otherwise
      error('no other type of dump can be read yet');
  end
  [X,Y,V] = readDumps(OFD, IT);

  % prepare output figure name.
%   savefigfullpath = [OFD,'comparedDump_IT',sprintf('%012.0f',IT)];
  savefigname = [savefigname, 'IT',sprintf('%012.0f',IT)];
  savefigfullpath = [savefigpath, regexprep(savefigname,'\.','')]; % regexprep because latex crashed when filenames have dots in it

  % interpolate dump
  delauTri = delaunayTriangulation(X,Y);
  tri = delauTri.ConnectivityList ;
  xi = delauTri.Points(:,1) ; yi = delauTri.Points(:,2) ;
  xi = (xi-min(xi)) / peak2peak(xi); % bring back x to [0, 1] (safety)
  yi = (yi-min(yi)) / peak2peak(yi); % bring back z to [0, 1] (safety)
  F = scatteredInterpolant(X,Y,V);
  zexp = F(xi,yi);
  caxxxxiiss_exp = [min(zexp),max(zexp)];

  % build analytic
  % cf. LNS_manufactured_solutions.mw
%   dEx=3; dEz=2;
  dEx=1.5; dEz=1.5;
  % dx=.01; x=0:dx:1; y=0:dx:1;
  % x = sort(unique(X));
  % y = sort(unique(Y));
  % [Xmg,Ymg]=meshgrid(x,y);
  % Vmg = sin(dEx*pi*Xmg)+sin(dEz*pi*Ymg); 
  zth = sin(dEx*pi*xi)+sin(dEz*pi*yi);
  caxxxxiiss_th = [min(zth),max(zth)];

  % error
  % zerr = zth-zexp; errName='Analytic - SPECFEM-DG-LNS';
  zerr = zexp-zth; errName='SPECFEM-DG-LNS - Analytic';
  caxxxxiiss_err = [min(zerr),max(zerr)];
  L2NormOfError = (integrate2DDelaunayTriangulation(delauTri, zerr.^2))^0.5;

  % colorbar scale
  all_caxx = [caxxxxiiss_th;caxxxxiiss_exp;caxxxxiiss_err];
%   glob_caxx = [min(all_caxx(:,1)),max(all_caxx(:,2))];
  glob_caxx = caxxxxiiss_th;

  %plot
  xlab = '$x/L$';
  ylab = '$z/L$';

  figure('units','normalized','outerposition',[0 0.4 1 0.6]);
  axxx = tight_subplot(1, 3, 0.012, [0.14,0.08], [0.05, 0.01]); % not mandatory, but prettier

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
  title(['SPECFEM-DG-LNS ($n=',num2str(IT),'$)']);

  axes(axxx(3));
  % axxx(3)=subplot(1,3,3);
  trisurf(tri,xi,yi,zerr); shading interp; view([0,0,1]); colormap(CMAP); caxis(glob_caxx);
  h_cb = colorbar;
  xlabel(xlab);
  yticklabels([]);
  title({[errName,', L$^2$-norm = ',num2str(L2NormOfError),'']});

  % Link = linkprop(axxx,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
  Link = linkprop(axxx,{'CameraPosition', 'XLim', 'YLim'});
  setappdata(gcf, 'StoreTheLink', Link);

  prettyAxes(gcf);
  customSaveFig(gcf, savefigfullpath,{'jpg','eps'});
end