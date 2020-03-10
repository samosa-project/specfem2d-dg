% Author:        Léo Martire.
% Description:   From an altitude-temperature-soundspeed-wind tuple, try to
%                find the IDs of possible atmospheric refractions.
% Notes:         N/A.
%
% Usage:
%   TODO,
% with:
%   TODO,
% yields:
%   TODO.

function [idpacks_fwd, idpacks_bwd, IDs_refr_T] = atmospheric_refractions_wrapper(spcfm_file)
  set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8); set(0, 'defaultTextFontSize', 18); set(0, 'defaultAxesFontSize', 18);
  set(0, 'DefaultTextInterpreter', 'latex'); set(0, 'DefaultLegendInterpreter', 'latex');
  
  [Z, ~, T, C, ~, ~, ~, ~, ~, ~, ~, NORTHWIND, EASTWIND, PROJWIND, ~, ~, ~, ~, ~] = extract_atmos_model(spcfm_file, 3, 0, 0);
  
  rad2deg=180/pi;
  
  proj = [];
  while(not(numel(proj)==1 && ((proj>=0 && proj<=360) || proj==-1) ))
    proj = input(['[',mfilename,'] Azimuth, in [deg] (',num2str(0*rad2deg),'=N, ',num2str(0.5*pi*rad2deg),'=E, ',num2str(1*pi*rad2deg),'=S, ',num2str(1.5*pi*rad2deg),'=W, -1 for already projected wind) ? > ']);
  end
  
  zsource = [];
  while(not(numel(zsource)==1 && ((zsource>=0 && zsource<=max(Z)) || zsource==-1) ))
    zsource = input(['[',mfilename,'] Source altitude, in [m]? > ']);
  end
  
  if(proj==-1)
    W = PROJWIND;
  else
%     W = cos(proj)*NORTHWIND - sin(proj)*EASTWIND;
    W = uv2azimuth(proj, NORTHWIND, EASTWIND);
  end
  
%   if(zsource==-1)
%     [ID_refr_W_forward, ID_refr_W_backward, IDs_refr_T] = atmospheric_refractions(Z, T, C, W, 1, ['Slice projected onto ',num2str(proj),'$^\circ$ counter-clockwise from North'], 100);
%   else
  sourcetxt = ['source at $z=',num2str(zsource),'$ [m]'];
  if(proj==-1)
    [idpacks_fwd, idpacks_bwd, IDs_refr_T] = atmospheric_refractions(Z, [], C, W, zsource, 1, {sourcetxt});
  else
    [idpacks_fwd, idpacks_bwd, IDs_refr_T] = atmospheric_refractions(Z, [], C, W, zsource, 1, {sourcetxt,['Slice projected onto azimuth ',num2str(proj),'$^\circ$']});
  end
%   [ID_refr_W_forward, ID_refr_W_backward, IDs_refr_T] = atmospheric_refractions(Z, T, C, W, zsource, 1, ['Slice projected onto ',num2str(proj),'$^\circ$ counter-clockwise from North'], );
%   end
end