% Author:        Léo Martire.
% Description:   Converts a 1D model from the SPECFEM2D-DG standard to the GeoAc standard.
% Notes:         TODO.
%
% Usage:
%   modelconv_specfem2geoac(spcfm_file, geoac_file)
% with:
%   spcfm_file    path to an input atmospheric model in the SPECFEM2D-DG format, 
%   ncpaprop_file path for the output atmospheric model to be formatted under the GeoAc format, 
%   maxz          [m] (optional) a maximum altitude at which truncate the input model (defaults to Inf), 
%   forceKeepC    (optional) a boolean (0 or 1) to ensure the speed of sound is conserved in this conversion at the cost of computing a sligthly different density (defaults to 0, meaning density is unchanged)
% yields:
%   N. A.

function [] = modelconv_specfem2ncpaprop(spcfm_file, ncpaprop_file, maxz, forceKeepC)
  if(nargin<2)
    error(['[', mfilename, ', ERROR] Not enough input arguments. Needs ''spcfm_file, geoac_file''.']);
  end
  if(not(exist('maxz', 'var')))
    maxz = Inf;
  end
  if(not(exist('forceKeepC', 'var')))
    forceKeepC = 0;
  end
  
  % Extract.
  if(forceKeepC)
    % Recompute a slightly off density, in order to keep speed of sound unchanged (see NPCAProp's AtmosphericProfile.cpp).
    % That is, ensure speed of sound in the input file is exactly the one understood by ncpaprop (see AtmosphericProfile.cpp), 
    %          at the cost of having a density ~4% (relative) wrong.
    [Z, ~, T, C, P, ~, ~, ~, ~, ~, ~, V, U, ~, ~, ~, ~] = extract_atmos_model(spcfm_file, 3, 0, 0);
    GAM_in_ncpaprop = 1.4; % this value of gamma is hard-coded in NCPAProp
    D = GAM_in_ncpaprop * P ./ (C.^2);
    % This ensures that when ncpaprop recomputes C as sqrt(GAM*P/D) (see AtmosphericProfile.cpp), C is the same as in the input file.
    disp(['[', mfilename, '] Ensuring C is conserved.']);
  else
    % Simply extract quantities and force feed them to the output file.
    [Z, D, T, ~, P, ~, ~, ~, ~, ~, ~, V, U, ~, ~, ~, ~] = extract_atmos_model(spcfm_file, 3, 0, 0);
    disp(['[', mfilename, '] C may not be conserved (because recomputed by ncpaprop using sqrt( gamma * P / rho )).']);
  end
  W = 0*U;
  
  % Convert units.
  maxz = maxz/1000;
  Z = Z / 1000; % convert m to km
  D = D / 1000; % convert kg/m^3 to g/cm^3
  P = P / 100; % convert Pa to hPa
  
  ZUVWTDP = [Z, U, V, W, T, D, P]';
  
	% output file
  ncpaprop_file = [ncpaprop_file, '_zuvwtdp.dat'];
  f_new = fopen(ncpaprop_file, 'w');
  if(f_new==-1)
    error(['[', mfilename, ', ERROR] Cannot open new data file ''', ncpaprop_file, '''.']);
  end
  
  % Selection.
  sel = find(Z<=min(maxz, max(Z)));
  fprintf(f_new, '%12.6e   %12.6e   %12.6e   %12.6e   %12.6e   %12.6e   %12.6e   \n', ZUVWTDP(:, sel));
  
  fclose(f_new);
  disp(strcat("[", mfilename, "] Model converted to: '", ncpaprop_file, "'."));
end