% Author:        Léo Martire.
% Description:   Converts a 1D model from the SPECFEM2D-DG standard to the GeoAc standard.
% Notes:         TODO.
%
% Usage:
%   modelconv_specfem2geoac(spcfm_file, geoac_file)
% with:
%   spcfm_file path to an input atmospheric model in the SPECFEM2D-DG format,
%   geoac_file path for the output atmospheric model to be formatted under the GeoAc format,
%   maxz       [m] (optional) a maximum altitude at which truncate the input model (defaults to Inf),
% yields:
%   N. A.

function [] = modelconv_specfem2geoac(spcfm_file, geoac_file, maxz)
  if(nargin<2)
    error(['[',mfilename,', ERROR] Not enough input arguments. Needs ''spcfm_file, geoac_file''.']);
  end
  
  if(not(exist('maxz', 'var'))); maxz=Inf; end
  
  [Z, RHO, T, ~, P, ~, ~, ~, ~, ~, ~, WN, WE, ~, ~, ~, ~] = extract_atmos_model(spcfm_file, 3, 0, 0);
  
  sel = find(Z<=min(maxz,max(Z)));
  
  Z   = Z/1000; % convert m to km
  RHO = RHO*1e-3; % convert kg/m^3 to g/cm^3
  P   = P * 1e-2; % convert Pa to mbar
  
  f_new = fopen(geoac_file, 'w');
  if(f_new==-1)
    error(strcat("[",mfilename,", ERROR] Cannot open new data file ", geoac_file,').'))
  end
  
  % z[km] T[K] w_M[m/s] w_Z[m/s] rho[kg/(m^3)] p[Pa]
  for ii=1:numel(sel)
    i=sel(ii);
    line_numbers=[Z(i), T(i), WE(i), WN(i), RHO(i), P(i)];
    fprintf(f_new,'%15.8e ',line_numbers);
    fprintf(f_new, "\n");
  end
  
  fclose(f_new);
  disp(['[',mfilename,'] Model ''',spcfm_file,''' converted to: ''', geoac_file, '''.']);
end