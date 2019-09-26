% Author:        Léo Martire.
% Description:   Given zonal and meridional winds, find the angle
% Notes:         Needs:
%                a) .m scripts and functions (if not alongside this
%                   script, recover via Léo):
%                  1) extract_atmos_model.m
%
% Usage:
%   [fh] = plot_model_wind(Z, WN, WE)
%   [fh] = plot_model_wind(Z, WN, WE, existingfignumber)
%   [fh] = plot_model_wind(Z, WN, WE, existingfignumber, colour)
%   [fh] = plot_model_wind(Z, WN, WE, existingfignumber, colour, labl)
% with:
%   TODO.
% yields:
%   TODO.

function [fh] = plot_model_wind(Z, V, U, existingfignumber, colour, labl)%, customcolorbar)
  if(not(exist('existingfignumber')))
    existingfignumber=-1;
  end
  if(not(exist('colour')))
    colour='k';
  end
  if(not(exist('labl')))
    labl=-1;
  end
%   if(not(exist('customcolorbar')))
%     customcolorbar=0;
%   end
  [AZIMUTH, W] = cart2pol(V,U);
%   AZIMUTH = AZIMUTH-pi/2; % Northward <-> theta=0. Eastward <-> theta=pi/2. Southward <-> theta=pi. Westward <-> theta=-pi/2.
  
  AZIMUTH(W==0)=0; % If wind is 0, zero the angle too.
  
  AZIMUTH=mod(AZIMUTH,2*pi);
  
  rad2deg=1; unit='[rad]';
  rad2deg=180/pi; unit='[deg]';
  
%   [WTH, W]

  WTHplot=rad2deg*AZIMUTH;
  
  axx = [];
  if(existingfignumber>0)
    % if a figure is provided, plot on it
    fh=figure(existingfignumber);
  else
    % else, create new figure
    fh=figure('units','normalized','outerposition',[0 0 1 1]);
  end
  
%   if(any(customcolorbar))
%     npanel=8;
%     subplot(npanel,2,1:2:(npanel-2)*2);
%   else
    subplot(121);
%   end
  axx = [axx, gca()];
  plthndl=plot(WTHplot,Z,'.','color',colour);
  if(existingfignumber>0); hold on; end; % if a figure is provided, hold on
  xlim(rad2deg*[0,2*pi]);
  ylim([min(Z), max(Z)]);
  ylabel(['$z$ [m]']);
  xlabel({['wind azimuth ',unit],['(',num2str(0*rad2deg),'=N, ',num2str(0.5*pi*rad2deg),'=E, ',num2str(1*pi*rad2deg),'=S, ',num2str(1.5*pi*rad2deg),'=W)']});
  xticks(rad2deg*2*pi*linspace(0,1,9));
  grid on; box on; set(gca,'ticklabelinterpreter','latex');
  if(labl>-1)
    % if a label is provided, apply it and add legend
    set(plthndl,'displayname',labl);
    ll=legend('location', 'east','fontsize',12);
    posll=get(ll,'position');
    set(ll,'position',[0.5-0.26*posll(3),posll(2:4)]);
  end
  
%   if(any(customcolorbar))
%     subplot(npanel,2,2:2:(npanel-2)*2);
%   else
    subplot(122);
%   end
  axx = [axx, gca()];
  plot(W,Z,'-','color',colour);
  if(existingfignumber>0); hold on; end; % if a figure is provided, hold on
  xlabel(['wind amplitude [m/s]']);
  yticklabels([]);
  grid on; box on; set(gca,'ticklabelinterpreter','latex');
  
  linkaxes(axx, 'y');
  
%   if(any(customcolorbar))
%     subplot(npanel,2,npanel*2+[-1,0]);
%     yticks([]);
%     xticks(customcolorbar);
%     box on; set(gca,'ticklabelinterpreter','latex');
%   end
end