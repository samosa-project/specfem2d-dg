% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   modelconv_specfem2ncpaprop(spcfm_file, ncpaprop_file)
%   modelconv_specfem2ncpaprop(spcfm_file, ncpaprop_file, maxz)
% with:
%   ncpaprop_file the output file (without any 'zuvwtdp' suffix, or extension)
%   maxz          (optional) the maximum altitude to output, [m] (default Inf)
% yields:
%   TODO.

function [] = modelconv_specfem2ncpaprop(spcfm_file, ncpaprop_file, maxz, forceKeepC)
  if(nargin<2)
    error(['[',mfilename,', ERROR] Not enough input arguments. Needs ''spcfm_file, geoac_file''.']);
  end
  if(not(exist('maxz','var')))
    maxz = Inf;
  end
  if(not(exist('forceKeepC','var')))
    forceKeepC = 0;
  end
  
  % Extract.
  if(forceKeepC)
    % Recompute a slightly off density, in order to keep speed of sound unchanged (see AtmosphericProfile.cpp).
    % That is, ensure speed of sound in the input file is exactly the one understood by ncpaprop (see AtmosphericProfile.cpp),
    %          at the cost of having a density ~4% (relative) wrong.
    [Z, ~, T, C, P, ~, ~, ~, ~, ~, ~, V, U, ~, ~, ~, ~] = extract_atmos_model(spcfm_file, 3, 0, 0);
%     GAM_in_ncpaprop = 1.4;
    D = GAM_in_ncpaprop * P ./ (C.^2);
    % This ensures that when ncpaprop recomputes C as sqrt(GAM*P/D) (see AtmosphericProfile.cpp), C is the same as in the input file.
    disp(['[] Ensuring C is conserved.']);
  else
    % Simply extract quantities and force feed them to the output file.
    [Z, D, T, ~, P, ~, ~, ~, ~, ~, ~, V, U, ~, ~, ~, ~] = extract_atmos_model(spcfm_file, 3, 0, 0);
    disp(['[] C may not be conserved (because recomputed by ncpaprop using sqrt( gamma * P / rho )).']);
  end
  W = 0*U;
  
  % Convert units.
  maxz = maxz/1000;
  Z = Z / 1000; % convert m to km
  D = D / 1000; % convert kg/m^3 to g/cm^3
  P = P / 100; % convert Pa to hPa
  
  ZUVWTDP = [Z,U,V,W,T,D,P]';
  
	% output file
  ncpaprop_file = [ncpaprop_file,'_zuvwtdp.dat'];
  f_new = fopen(ncpaprop_file, 'w');
  if(f_new==-1)
    error(strcat("[",mfilename,", ERROR] Cannot open new data file ", ncpaprop_file,').'))
  end
  
  % Selection.
  sel = find(Z<=min(maxz,max(Z)));
  fprintf(f_new, '%12.6e   %12.6e   %12.6e   %12.6e   %12.6e   %12.6e   %12.6e   \n', ZUVWTDP(:,sel));
  
  fclose(f_new);
  disp(strcat("[",mfilename,"] Model converted to: '", ncpaprop_file, "'."));
end