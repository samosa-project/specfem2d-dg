% Author:        Léo Martire.
% Description:   Given zonal and meridional winds, find the angle
% Notes:         Needs:
%                a) .m scripts and functions (if not alongside this
%                   script, recover via Léo):
%                  1) extract_atmos_model.m
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function [] = plot_wind(Z,WN,WE)
  [WTH,W]=cart2pol(WE,WN);
  WTH=WTH-pi/2; % Northward <-> theta=0. Eastward <-> theta=-pi/2. Southward <-> theta=pi. Westward <-> theta=pi/2.
  
  WTH(W==0)=0; % If wind is 0, zero the angle too.
  
  WTH=mod(WTH,2*pi);
  
  rad2deg=1; unit='[rad]';
  rad2deg=180/pi; unit='[deg]';
  
%   [WTH, W]

  WTHplot=rad2deg*WTH;
  figure('units','normalized','outerposition',[0 0 1 1]);
  subplot(121);
  plot(WTHplot,Z,'.');
  xlim(rad2deg*[0,2*pi]);
  ylabel(['$z$ [m]']);
  xlabel({['wind angle (from North, counter-clockwise) ',unit],['(',num2str(0*rad2deg),'=N, ',num2str(0.5*pi*rad2deg),'=W, ',num2str(1*pi*rad2deg),'=S, ',num2str(1.5*pi*rad2deg),'=E)']});
  xticks(rad2deg*2*pi*linspace(0,1,9));
  grid on; box on; set(gca,'ticklabelinterpreter','latex');
  subplot(122);
  plot(W,Z,'-');
  xlabel(['wind amplitude [m/s]']);
  grid on; box on; set(gca,'ticklabelinterpreter','latex');
end