% Author:        Léo Martire, Guerman Poler.
% Description:   Retrieves and computes useful quantities from an ERA5 file, previously downloaded from ECMWF's servers.
% Last modified: See file metadata.
% Usage:         1) Make sure an ERA5 file was previously downloaded.
%                2) Call the function with the path to that ERA5 file as argument.
% Notes:         If this function fails, ncdisp(ecmwf_era5_file) can help debugging.

function [points,z,T,p_half,p_full,g,w_M,w_Z]=retrieve_ECMWF(ecmwf_era5_file)
  % Default paths.
  ecmwf_coefs='ECMWF/ecmwf_coeffs.dat';

  % Find necessary components.
  % Get this function's folder.
  this_func_path_splitted=split(mfilename('fullpath'),'/'); this_func_path_splitted{end}=''; this_func_dir=join(this_func_path_splitted,'/'); clear('this_func_path_splitted');
  % Normally, this function is alongside this function, thus the full path is the following:
  ecmwf_coefs_fullpath = char(strcat(this_func_dir,ecmwf_coefs));
  
  % Constants
  Rd = 287.06; % ??
  g_0 = 9.80665; % Gravitional pull.
  GM_e = (6.67e-11)*(5.97219e24); % G * Earth mass.
  R_e = 6371010.0; % Earth radius.

  % Get model coefficients.
  ecmwf_coefs = load(ecmwf_coefs_fullpath);
  ecmwf_coef_a = ecmwf_coefs(:,2); % Coefficient A du modèle.
  ecmwf_coef_b = ecmwf_coefs(:,3); % Coefficient B du modèle.
  clear('ecmwf_coefs');

  % Read points and dimensions.
  long = ncread(ecmwf_era5_file, 'longitude');
  lat = ncread(ecmwf_era5_file, 'latitude');
  lvl = ncread(ecmwf_era5_file, 'level');
  time = ncread(ecmwf_era5_file, 'time');

  % Save dimensions.
  dimensions_len = [length(long), length(lat), length(lvl), length(time)]; % Get number of points for each dimension.
  dimension_of_levels=3; % Save dimension ID where levels are (for flipping later).

  % Read quantities.
  lnsp_full = ncread(ecmwf_era5_file, 'lnsp'); % Logarithme de pression surfacique
  z_full = ncread(ecmwf_era5_file, 'z'); % Géopotentiel au niveau du sol
  q = ncread(ecmwf_era5_file, 'q'); % Humidité spécifique
  T = ncread(ecmwf_era5_file, 't'); % Température
  w_Z = ncread(ecmwf_era5_file, 'u'); % Vent zonal
  w_M = ncread(ecmwf_era5_file, 'v'); % Vent méridional

  % Change of the model indexing to a more natural reading direction (ground to sky).
  % For original model, level 137 = surface.
  % For new model, level 1 = surface.
  q=flip(q,dimension_of_levels);
  T=flip(T,dimension_of_levels);
  w_Z=flip(w_Z,dimension_of_levels);
  w_M=flip(w_M,dimension_of_levels);
  
  % Computation of geopotential and pressure on model levels (based on model coefficients).
  disp(['  [',mfilename,'] Computation of geopotential and pressure on model levels (based on model coefficients).']);
  % Half pressure is at midpoint of level. Full pressure is at top of level.
  % Algorithm below is from 'compute_geopotential.py', found at
  % https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+geopotential+on+model+levels.
  p_half = zeros(dimensions_len); 
  p_full = zeros(dimensions_len);
  sp = exp(lnsp_full(:,:,1,:));
  z_h = z_full(:,:,1,:);
  for k=1:dimensions_len(3)
    t_moist = T(:,:,k,:).*(1 + 0.609133*q(:,:,k,:));
    p_half(:,:,k,:) = ecmwf_coef_a(end-(k-1)) + ecmwf_coef_b(end-(k-1))*sp;
    if (k == dimensions_len(3))
      % Last level.
      dlogP = log(p_half(:,:,k,:)/0.1);
      alpha = log(2);
      p_full(:,:,k,:) = 1/2*(p_half(:,:,k,:) + 0.1);
    else
      % Intermediate level.
      ph_lvlplusone = ecmwf_coef_a(end-k) + ecmwf_coef_b(end-k)*sp;
      dlogP = log(p_half(:,:,k,:)./ph_lvlplusone);
      dP = p_half(:,:,k,:)-ph_lvlplusone;
      alpha = 1 - (ph_lvlplusone./dP).*dlogP;
      p_full(:,:,k,:) = 1/2*(ph_lvlplusone + p_half(:,:,k,:));
    end
    TRd = t_moist*Rd;
    z_full(:,:,k,:) = z_h + TRd.*alpha;
    z_h = z_h + TRd.*dlogP;
  end

  % Conversion of time from ECMWF format (hours since 01/01/1900 00:00:00) to Matlab's datenum format (days since 01/01/0000 00:00:00).
  time = double(time)/24 + datenum('01/01/1900 00:00:00');
  
  % Computation of the effective altitude from geopotential altitude.
  z = z_full/g_0;
  % Computation of the effective gravity, based on geopotential height.
  g = GM_e./(R_e + z).^2;
  
  % Save array of points where model was computed.
  points={long,lat,lvl,time};
  disp(['  [',mfilename,'] Model loaded. Points={longitude,latitude,level,time}.']);
end