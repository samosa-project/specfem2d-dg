% Author:        Léo Martire.
% Description:   Builds a 1D atmospheric model file using MSISE and ECMWF models, and prints it in a format compatible with SPECFEM2D-DG.
%                The components are as follows:
%                  From ECMWF exclusively: T, 
%                                          p, 
%                                          g, 
%                                          meridional wind, 
%                                          zonal wind.
%                  From MSISE exclusively: rho, 
%                                          scale_height, 
%                                          kappa, 
%                                          mu, mu_vol, 
%                                          c_p, c_v, 
%                                          gamma.
%                  Computed from both:     soundspeed (using gamma MSISE, p ECMWF, and rho MSISE), 
%                                          N^2 (gamma MSISE, g ECMWF, rho MSISE, and p ECMWF).
% Notes:         An ERA5 file must have previously been downloaded. See the provided Matlab script 'prepare_call_ECMWF_API' in order to do so.
%                Anyhow, one cannot hope to make a single script performing both the API query and the model output, because the ECMWF queries tend to be relatively slow.
%
% Usage:         1) Download an ERA5 from ECMWF's API (use the function 'prepare_call_ECMWF_API').
%                2) Call this script.

clear all;
% close all;
clc;

[SPCFMEXloc] = setup_overall();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametrisation.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECMWF_DATAFILE = '/home/l.martire/Downloads/Atmospheric_model/Data/ERA5/era5.nc';
ECMWF_DATAFILE = input(['[', mfilename, '] Input ERA5 file to use > '], 's');
if(not(exist(ECMWF_DATAFILE, 'file')))
  error(['[', mfilename, ', ERROR] This ERA5 file does not exist.']);
end
threshold_ok_latlon = 0.5; % If point is < 0.5 ° away, consider it ok.
threshold_ok_time = 10; % If time is < 10 minutes away, consider it ok.
R = 1.380649e-23 * 6.02214076e23; % gas constant = k * N_A
M_dryair = 0.0288576070373303; % molar mass of dry air on Earth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin treatment,        %
% extraction, and output.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load ERA5.
disp(['[', mfilename, '] Selected ERA5 file is: ', ECMWF_DATAFILE]);
[pts, z_ecmwf, T_ecmwf, p_half, p_full_ecmwf, g_ecmwf, w_M_ecmwf, w_Z_ecmwf] = retrieve_ECMWF(ECMWF_DATAFILE);
ecmwf_lon = pts{1};
ecmwf_lat = pts{2};
ecmwf_time = pts{4};
[era5lonspan, era5latspan, era5timespan] = retrieve_ECMWF_span(pts);

% Ask user input.
aim_long = input(['[', mfilename, '] Input wanted longitude (must be in ', era5lonspan, ') > ']);
aim_lat = input(['[', mfilename, '] Input wanted latitude (must be in ', era5latspan, ') > ']);
time_ok = -1;
while(not(time_ok == 1))
  aim_time = datenum(input(['[', mfilename, '] Input wanted date and time (format ''YYYY/[M]M/[D]D [[H]H:[M]M:[S]S]'', must be in ', era5timespan, ') > '], 's'));
  time_ok = input(['[', mfilename, '] Wanted time is ', datestr(aim_time), ' UT, is that ok (1 for yes, 0 for no)? > ']);
end

aim_long = mod(aim_long, 360); % Put it in [0, 360].

% Check user choices are actually possible within the selected ERA5 file.
minimum_lon_gap = min(abs(ecmwf_lon-aim_long));
long_found = minimum_lon_gap<threshold_ok_latlon;
minimum_lat_gap = min(abs(ecmwf_lat-aim_lat));
lat_found = minimum_lat_gap<threshold_ok_latlon;
minimum_time_gap = min(abs(datetime(datevec(ecmwf_time))-datetime(datevec(aim_time))));
time_found = minutes(minimum_time_gap)<threshold_ok_time;

disp(['[', mfilename, '] Closest (lon, lat) point found is (', num2str(minimum_lon_gap), ', ', num2str(minimum_lat_gap), ') away from wanted (lon, lat) point.']);
disp(['[', mfilename, '] Closest timestamp found is ', char(minimum_time_gap), ' away from wanted timestamp.']);

if(not(lat_found && long_found && time_found))
  error(['[', mfilename, ', ERROR] Point not found in ECMWF data, or too far.']);
else
  disp(['[', mfilename, '] Point found in ECMWF data, or close enough.']);
end

% Choose the best data point from the ERA5 file.
id_lon = find(abs(aim_long-ecmwf_lon) == min(abs(aim_long-ecmwf_lon))); 
id_lat = find(abs(aim_lat-ecmwf_lat) == min(abs(aim_lat-ecmwf_lat)));
id_time = find(abs(aim_time-ecmwf_time) == min(abs(aim_time-ecmwf_time)));
final_lat = ecmwf_lat(id_lat);
final_lon = ecmwf_lon(id_lon);
final_time = ecmwf_time(id_time);
final_time_vec = datevec(final_time);
annee = final_time_vec(1)-2000;
jours = day(datetime(final_time_vec), 'dayofyear');
secondes = sum(final_time_vec(4:end).*[3600, 60, 1]);

% Remove useless points, keep only the one we're interested in.
z_ecmwf = squeeze(z_ecmwf(id_lon, id_lat, :, id_time));
T_ecmwf = squeeze(T_ecmwf(id_lon, id_lat, :, id_time));
p_half = squeeze(p_half(id_lon, id_lat, :, id_time));
p_full_ecmwf = squeeze(p_full_ecmwf(id_lon, id_lat, :, id_time));
g_ecmwf = squeeze(g_ecmwf(id_lon, id_lat, :, id_time));
w_M_ecmwf = squeeze(w_M_ecmwf(id_lon, id_lat, :, id_time));
w_Z_ecmwf = squeeze(w_Z_ecmwf(id_lon, id_lat, :, id_time));

p_ecmwf = p_full_ecmwf;

if(z_ecmwf(1)<0)
  error(['[', mfilename, ', ERROR] First ECMWF level has negative altitude. You have to check what that means, and correct this function accordingly.']);
end

dz_ecmwf = min(diff(z_ecmwf)); % Get finest step from ECMWF.
min_alt = input(['[', mfilename, '] Minimum altitude (m)? > ']);
max_alt = input(['[', mfilename, '] Maximum altitude (m, <', num2str(z_ecmwf(end)), ')? > ']);
nsteps = ceil((max_alt-min_alt)/dz_ecmwf);
nsteps_input = input(['[', mfilename, '] About to call MSISEHWM wrapper with ', num2str(nsteps), ' layers (dz = ', num2str((max_alt-min_alt)/nsteps), '). Continue? 1 for yes; new number of layers if wanted (max 5000)? > ']);
if(not(isempty(nsteps_input) | nsteps_input == 1))
  nsteps = nsteps_input;
end
clear('nsteps_input');

% wind_proj = input(['[', mfilename, '] Input projection angle for projected wind (0 is full zonal positive eastward, 90 is full meridional positive northward, 180 is full backward zonal positive southward, 270 is full backward meridional positive westward)? > ']);
wind_proj = input(['[', mfilename, '] Input projection angle for projected wind ([deg] from North counter-clockwise)? > ']);
% run plot_model_wind(z_ecmwf, w_M_ecmwf, w_Z_ecmwf).
flipwind = -1;
while(not(numel(flipwind) == 1 & ismember(flipwind, [0, 1])))
  disp(['[', mfilename, '] Flip wind? 0 for positive wind in West->East (/South->North) direction; 1 for positive in East->West (/North->South).']);
  flipwind = input(['[', mfilename, ']            This is still physical, it simply means we look on either side of the simulation plane. > ']);
end

o_file = model_namer(min_alt, max_alt, nsteps, final_lat, final_lon, annee, jours, secondes, wind_proj);

o_folder = '.';
wind_project_angle = 0;
build_MSISEHWM(min_alt, max_alt, nsteps, final_lat, final_lon, annee, jours, secondes, o_file, o_folder);
% Remark: wind_proj not used in call to MSISEHWM since we will use the winds from ECMWF.
[Z, ~, ~, ~, ~, H_m, ~, ~, ~, ~, MUVOL_m, ~, ~, ~, CP_m, CV_m, GAMMA_m] = extract_atmos_model(o_file, 3, 0, 0);

% Interpolate ECMWF on MSISE grid.
disp(['[', mfilename, '] Interpolating ECMWF quantities on MSISE (finer) grid with splines.']);
T_e = spline(z_ecmwf, T_ecmwf, Z);
p_e = spline(z_ecmwf, p_full_ecmwf, Z);
g_e = spline(z_ecmwf, g_ecmwf, Z);
w_North_e = spline(z_ecmwf, w_M_ecmwf, Z);
W_East_e = spline(z_ecmwf, w_Z_ecmwf, Z);

disp(['[', mfilename, '] Summary: loaded (Z, H, MUVOL, CP, CV, GAMMA) from MSISE, .']);
disp([blanks(length(mfilename)+2), '             and (T, P, G, WN, WE)            from ECMWF.']);
  
% Projecting wind.
disp(['[', mfilename, '] Built W by projection of (WN, WE) onto ', num2str(wind_proj), '° from North counter-clockwise.']); % Projection angle, from North, counter-clockwise, [rad].
% w_P_e = cos(wind_proj*pi/180.)*W_Z_e+sin(wind_proj*pi/180.)*w_M_e;
w_P_e = w_North_e*cos(wind_proj*pi/180) - W_East_e*sin(wind_proj*pi/180);
if(flipwind)
  w_P_e = -w_P_e;
end
smthpar = 1e-9; spline = fit(Z, w_P_e, 'smoothingspline', 'smoothingparam', smthpar); nW = spline(Z); disp(['[', mfilename, '] Spline W (par = ', num2str(smthpar), ').']);
% Apodise wind.
apodisewind = -1;
while(not(numel(apodisewind) == 1 & ismember(apodisewind, [0, 1, 2])))
  apodisewind = input(['[', mfilename, '] Apodise wind (0 for no, 1 for yes and w(z = 0) = 0, 2 for yes and w(z = 0) unchanged from model)? > ']);
end
if(apodisewind)
  width = -1;
  while(not(numel(width) == 1 & width>0 & width<max(Z)-min(Z)))
    width = input(['[', mfilename, '] Apodisation width (>0, <', num2str(max(Z)-min(Z)), ')? > ']);
  end
%   b = -1.82138636771844967304021031862099524348122888360095; % 1e-2 level.
  b = -2.75106390571206079614551316854267817109677902559646; % 1e-3 level.
%   b = -3.45891073727950002215092763595756951991566980804288674707621013; % 1e-6 level.
  a = -2*b/width; apoWind = erf(a*Z+b)/2+0.5; apoWind(Z == min(Z)) = 0;
  if(apodisewind == 1)
    nW = apoWind.*nW; disp(['[', mfilename, '] Apodised W up to ', num2str(width), ' m such that W(Z == ', num2str(min(Z)), ') = 0.']);
  elseif(apodisewind == 2)
    nW = w_P_e+apoWind.*(nW-w_P_e);
    dzmovingavg = 100; nW(Z<=width) = smooth(Z(Z<=width), nW(Z<=width), 'moving', floor(dzmovingavg/mean(diff(Z(Z<=width)))));
    disp(['[', mfilename, '] Apodise W up to ', num2str(width), ' m such that W(Z == ', num2str(min(Z)), ') is unchanged. Also applied a moving average on [', num2str(min(Z)), ', ', num2str(width), '] m with width ', num2str(dzmovingavg), ' m.']);
  else
    error('kek');
  end
end
w_P_treated = nW;

% Compute hydrostatic rho from ECMWF p.
P = p_e; G = g_e;
nRHO = - differentiation_matrix(Z, 0) * P ./ G; % $\rho = -\partial_z{P} / g_z$
disp(['[', mfilename, '] Bruteforce RHO from ECMWF''s P.']); % Regularise hydrostatic ratio by bruteforcing $\rho = -\partial_z{P} / g_z$.
smthpar = 1e-12; spline = fit(Z, log(nRHO), 'smoothingspline', 'smoothingparam', smthpar); nRHO = exp(spline(Z)); disp(['[', mfilename, '] Spline log(RHO) (par = ', num2str(smthpar), ').']);
% smthpar = 1e-15; spline = fit(Z, LRHO, 'smoothingspline', 'smoothingparam', smthpar);       figure();subplot(121);plot(LRHO, Z, spline(Z), Z);subplot(122);plot(D*LRHO, Z, D*spline(Z), Z);
rho_treated = nRHO;

% Compute dynamic viscosity from ECMWF quantities.
disp(['[', mfilename, '] Computing MU and KAPPA from ECMWF''s (RHO, T, P).']);
TP = thermodynamicalParameters;
mu_e = TP.mu.idealN2(rho_treated, T_e, p_e);
kappa_e = TP.kappa.ussa76(T_e);

% Compute missing quantities:
% - sound speed can be computed from T_ECMWF and gamma_MSISE (c = sqrt(gamma*R*T/M) = sqrt(gamma*P/rho)), 
% - N^2 can be computed from gamma_MSISE, g_ECMWF, and T_ECMWF (N^2 = (gamma-1)*g^2/(gamma*R*T)).
disp(['[', mfilename, '] Computing missing quantities (sound speed, Brunt-Väisälä frequency) from MSISE''s GAMMA, ECMWF''s (P, RHO, G).']);
% Version 1, assuming the molar mass of air does not change and is equal to dry air molar mass. Relying on 1 MSISE quantity.
% soundspeed_merged = sqrt(GAMMA.*R.*T_e2m/M_dryair);
% Nsquared_merged = (GAMMA-1).*g_e2m.^2.*M_dryair./(GAMMA.*R.*T_e2m);
% Version 2, more exact. Relying on 2 MSISE quantities.
% soundspeed_merged = sqrt(GAMMA.*p_e./DENSITY); % unused by specfem
soundspeed_merged = sqrt(GAMMA_m.*p_e./rho_treated); % unused by specfem
Nsquared_merged = (GAMMA_m-1).*g_e.^2.*rho_treated./(GAMMA_m.*p_e); % unused by specfem

killMU = -1;
while(not(numel(killMU) == 1 & ismember(killMU, [0, 1])))
  killMU = input(['[', mfilename, '] Set MU (viscosity) to zero (0 for no, 1 for yes)? > ']);
end
killKAPPA = -1;
while(not(numel(killKAPPA) == 1 & ismember(killKAPPA, [0, 1])))
  killKAPPA = input(['[', mfilename, '] Set KAPPA (thermal conductivity) to zero (0 for no, 1 for yes)? > ']);
end
if(killMU)
  factorMU = 0; suffixMU = '_noMU';
  disp(['[', mfilename, '] Killing MU (MU is now 0 everywhere).']);
else
  factorMU = 1; suffixMU = '';
end
if(killKAPPA)
  factorKAPPA = 0; suffixKAPPA = '_noKAPPA';
  disp(['[', mfilename, '] Killing KAPPA (KAPPA is now 0 everywhere).']);
else
  factorKAPPA = 1; suffixKAPPA = '';
end

plotricharddd = -1;
while(not(numel(plotricharddd) == 1 & ismember(plotricharddd, [0, 1])))
  plotricharddd = input(['[', mfilename, '] Plot Richardson number (0 for no, 1 for yes)? > ']);
end
if(plotricharddd)
  figure();
  semilogx(Richardson_number(w_P_e, differentiation_matrix(Z, 0), Nsquared_merged.^0.5), Z); hold on;
  semilogx(Richardson_number(w_P_treated, differentiation_matrix(Z, 0), Nsquared_merged.^0.5), Z); hold on;
  plot(0.25*[1, 1], [min(Z), max(Z)]); hold on;
  plot(1*[1, 1], [min(Z), max(Z)]);
  xlim([1e-1, 1e2]);
  ylim([min(Z), max(Z)]);
end

% Output newly-built file.
disp(['[', mfilename, '] Rewriting model to another file.']);
new_o_file = ['msiseecmwf_', o_file, suffixMU, suffixKAPPA];
kek = what(o_folder);
fullnewpath = [kek.path, filesep, new_o_file];
% rewrite_atmos_model(new_o_file, o_file, ALTITUDE, DENSITY, T_e2m, soundspeed_merged, p_e2m, LOCALPRESSURESCALE, g_e2m, Nsquared_merged, KAPPA, MU, MUVOL, w_M_e2m, W_Z_e2m, w_P, CP, CV, GAMMA);
rewrite_atmos_model(fullnewpath, o_file, Z, rho_treated, T_e, ...
                    soundspeed_merged, p_e, H_m, g_e, Nsquared_merged, ...
                    factorKAPPA*kappa_e, factorMU*mu_e, factorMU*MUVOL_m, ...
                    w_North_e, W_East_e, w_P_treated, CP_m, CV_m, GAMMA_m);

disp(['[', mfilename, '] Model stored in: ''', fullnewpath, '''.']);
plot_model(fullnewpath, '-', 'k', []);