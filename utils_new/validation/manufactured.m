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
close all;
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
addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/standalone');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output files
OFD = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/validation_lns_manufactured/OUTPUT_FILES_dp';

% IT = 5;
% IT = 5000;
IT = 10000;
% IT = 15000;
% IT = 20000;

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
figPath = [OFD,'comparedDump_IT',sprintf('%012.0f',IT)];

% interpolate dump
delauTri = delaunayTriangulation(X,Y);
tri = delauTri.ConnectivityList ;
xi = delauTri.Points(:,1) ; yi = delauTri.Points(:,2) ;
F = scatteredInterpolant(X,Y,V);
zexp = F(xi,yi);

% build analytic
% cf. LNS_manufactured_solutions.mw
dEx=3;
dEz=2;
% dx=.01; x=0:dx:1; y=0:dx:1;
% x = sort(unique(X));
% y = sort(unique(Y));
% [Xmg,Ymg]=meshgrid(x,y);
% Vmg = sin(dEx*pi*Xmg)+sin(dEz*pi*Ymg); 
zth = sin(dEx*pi*xi)+sin(dEz*pi*yi);

% error
% zerr = zth-zexp; errName='Analytic - SPECFEM-DG-LNS';
zerr = zexp-zth; errName='SPECFEM-DG-LNS - Analytic';
L2NormOfError = (integrate2DDelaunayTriangulation(delauTri, zerr.^2))^0.5;

%plot
figure('units','normalized','outerposition',[0 0.5 1 0.5]);

axxx(1)=subplot(1,3,1);
% surf(Xmg,Ymg,Vmg); shading interp; view([0,0,1]); colormap jet;
trisurf(tri,xi,yi,zth); shading interp; view([0,0,1]); colormap jet; colorbar;
title('Analytic');

axxx(2)=subplot(1,3,2);
trisurf(tri,xi,yi,zexp); shading interp; view([0,0,1]); colormap jet; colorbar;
yticklabels([]);
title('SPECFEM-DG-LNS');

axxx(3)=subplot(1,3,3);
trisurf(tri,xi,yi,zerr); shading interp; view([0,0,1]); colormap jet; colorbar;
yticklabels([]);
title({[errName,', L$^2$-norm = ',num2str(L2NormOfError),'']});

% Link = linkprop(axxx,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
Link = linkprop(axxx,{'CameraPosition', 'XLim', 'YLim'});
setappdata(gcf, 'StoreTheLink', Link);

prettyAxes(gcf);
customSaveFig(figPath,{'jpg','png'});