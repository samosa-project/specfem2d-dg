% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   modelconv_specfem2geoac(spcfm_file, geoac_file)
% with:
%   TODO.
% yields:
%   TODO.

function [] = modelconv_specfem2geoac(spcfm_file, geoac_file, maxz)

  if(nargin<2)
    error(['[',mfilename,', ERROR] Not enough input arguments. Needs ''spcfm_file, geoac_file''.']);
  end
  
  if(not(exist('maxz')))
    maxz=Inf;
  end
  
%   format compact;
%   set(0, 'DefaultLineLineWidth', 3); set(0, 'DefaultLineMarkerSize', 8);
%   set(0, 'defaultTextFontSize', 14); set(0, 'defaultAxesFontSize', 14);
%   set(0, 'DefaultTextInterpreter', 'latex');
%   set(0, 'DefaultLegendInterpreter', 'latex');
  
  [Z, RHO, T, ~, P, ~, ~, ~, ~, ~, ~, WN, WE, ~, ~, ~, ~] = extract_atmos_model(spcfm_file, 3, 0, 0);
  
  sel = find(Z<=min(maxz,max(Z)));
  
  Z   = Z/1000; % convert m to km
  RHO = RHO*1e-3; % convert kg/m^3 to g/cm^3
  P   = P * 1e-2; % convert Pa to mbar
  
%   geoac_file
  f_new = fopen(geoac_file, 'w');
  if(f_new==-1)
    error(strcat("[",mfilename,", ERROR] Cannot open new data file ", geoac_file,').'))
  end
  
  % z[km] T[K] w_M[m/s] w_Z[m/s] rho[kg/(m^3)] p[Pa]
  for ii=1:numel(sel)
    i=sel(ii);
    line_numbers=[Z(i), T(i), WN(i), WE(i), RHO(i), P(i)];
    fprintf(f_new,'%15.8e ',line_numbers);
    fprintf(f_new, "\n");
  end
  
  fclose(f_new);
  disp(strcat("[",mfilename,"] Model converted to: '", geoac_file, "'."));
end