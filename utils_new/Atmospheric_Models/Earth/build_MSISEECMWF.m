% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Builds an atmospheric model file (directly compatible
%                with SPECFEM) using MSISE and ECMWF models.
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
%                  Computed from both: soundspeed (gamma MSISE, p ECMWF,
%                                                  rho MSISE),
%                                      N^2 (gamma MSISE, g ECMWF,
%                                           rho MSISE, p ECMWF).
% Last modified: See file metadata.
% Usage:         1) Download an ERA5 from ECMWF's API (use the function 'prepare_call_ECMWF_API').
%                2) Call this script.
% Notes:         An ERA5 file must have previously been downloaded. See
%                the provided Matlab script 'prepare_call_ECMWF_API' in
%                order to do so.
%                See LaTeX project
%                https://v2.overleaf.com/project/5c4827468093415759c04743.

clear all;
close all;
clc;

tmp=evalc('which extract_atmos_model');tmp=split(tmp);tmp=tmp{1};tmp=split(tmp,'/');tmp{end}='';tmp=join(tmp,'/');tmp=tmp{1};addpath(tmp); clear('tmp'); % Addpath for the 'extract_atmos_model' function.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametrisation.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECMWF_DATAFILE='/home/l.martire/Downloads/Atmospheric_model/Data/ERA5/era5.nc';
ECMWF_DATAFILE=input(['[',mfilename,'] Input ERA5 file to use > '],'s');
threshold_ok_latlon=0.5; % If point is < 0.5 ° away, consider it ok.
threshold_ok_time=10; % If time is < 10 minutes away, consider it ok.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=8.3144598;
M_dryair=2.89645e-2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin treatment,            %
% extraction, and output.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load ERA5.
disp(['[',mfilename,'] Selected ERA5 file is: ',ECMWF_DATAFILE]);
[pts,z_ecmwf,T_ecmwf,p_half,p_full_ecmwf,g_ecmwf,w_M_ecmwf,w_Z_ecmwf]=retrieve_ECMWF(ECMWF_DATAFILE);
ecmwf_lon=pts{1};
ecmwf_lat=pts{2};
ecmwf_time=pts{4};
[era5lonspan,era5latspan,era5timespan] = retrieve_ECMWF_span(pts);

% Ask user input.
aim_long=input(['[',mfilename,'] Input wanted longitude (must be in ',era5lonspan,') > ']);
aim_lat=input(['[',mfilename,'] Input wanted latitude (must be in ',era5latspan,') > ']);
time_ok=-1;
while(not(time_ok==1))
  aim_time=datenum(input(['[',mfilename,'] Input wanted date and time (format ''YYYY/[M]M/[D]D [[H]H:[M]M:[S]S]'', must be in ',era5timespan,') > '],'s'));
  time_ok=input(['[',mfilename,'] Wanted time is ',datestr(aim_time),', is that ok (1 for yes, 0 for no)? > ']);
end

aim_long=mod(aim_long,360); % Put it in [0,360].

% Check user choices are actually possible within the selected ERA5 file.
minimum_lon_gap=min(abs(ecmwf_lon-aim_long));
long_found=minimum_lon_gap<threshold_ok_latlon;
minimum_lat_gap=min(abs(ecmwf_lat-aim_lat));
lat_found=minimum_lat_gap<threshold_ok_latlon;
minimum_time_gap=min(abs(datetime(datevec(ecmwf_time))-datetime(datevec(aim_time))));
time_found=minutes(minimum_time_gap)<threshold_ok_time;

disp(['[',mfilename,'] Closest (lon, lat) point found is (',num2str(minimum_lon_gap),', ',num2str(minimum_lat_gap),') away from wanted (lon, lat) point.']);
disp(['[',mfilename,'] Closest timestamp found is ',char(minimum_time_gap),' away from wanted timestamp.']);

if(not(lat_found && long_found && time_found))
  error(['[',mfilename,', ERROR] Point not found in ECMWF data, or too far.']);
else
  disp(['[',mfilename,'] Point found in ECMWF data, or close enough.']);
end

% Choose the best data point from the ERA5 file.
id_lon=find(abs(aim_long-ecmwf_lon)==min(abs(aim_long-ecmwf_lon))); 
id_lat=find(abs(aim_lat-ecmwf_lat)==min(abs(aim_lat-ecmwf_lat)));
id_time=find(abs(aim_time-ecmwf_time)==min(abs(aim_time-ecmwf_time)));
final_lat=ecmwf_lat(id_lat);
final_lon=ecmwf_lon(id_lon);
final_time=ecmwf_time(id_time);
final_time_vec=datevec(final_time);
annee=final_time_vec(1)-2000;
jours=day(datetime(final_time_vec), 'dayofyear');
secondes=sum(final_time_vec(4:end).*[3600,60,1]);

% Remove useless points, keep only the one we're interested in.
z_ecmwf=squeeze(z_ecmwf(id_lon,id_lat,:,id_time));
T_ecmwf=squeeze(T_ecmwf(id_lon,id_lat,:,id_time));
p_half=squeeze(p_half(id_lon,id_lat,:,id_time));
p_full_ecmwf=squeeze(p_full_ecmwf(id_lon,id_lat,:,id_time));
g_ecmwf=squeeze(g_ecmwf(id_lon,id_lat,:,id_time));
w_M_ecmwf=squeeze(w_M_ecmwf(id_lon,id_lat,:,id_time));
w_Z_ecmwf=squeeze(w_Z_ecmwf(id_lon,id_lat,:,id_time));

p_ecmwf=p_full_ecmwf;

if(z_ecmwf(1)<0)
  error(['[',mfilename,', ERROR] First ECMWF level has negative altitude. You have to check what that means, and correct this function accordingly.']);
end

dz_ecmwf=min(diff(z_ecmwf)); % Get finest step from ECMWF.
min_alt=input(['[',mfilename,'] Minimum altitude (m)? > ']);
max_alt=input(['[',mfilename,'] Maximum altitude (m, <',num2str(z_ecmwf(end)),')? > ']);
nsteps=ceil((max_alt-min_alt)/dz_ecmwf);
nsteps_input=input(['[',mfilename,'] About to call MSISEHWM wrapper with ',num2str(nsteps),' layers (dz=',num2str((max_alt-min_alt)/nsteps),'). Continue (1 for yes, new value for no)? > ']);
if(not(isempty(nsteps_input) | nsteps_input==1))
  nsteps=nsteps_input;
end
clear('nsteps_input');

wind_proj=input(['[',mfilename,'] Input projection angle for projected wind (0 is full zonal positive eastward, 90 is full meridional positive northward, 180 is full backward zonal positive southward, 270 is full backward meridional positive westward)? > ']);

o_file=model_namer(min_alt, max_alt, nsteps, final_lat, final_lon, annee, jours, secondes, wind_proj);

o_folder='.';
wind_project_angle=0;
build_MSISEHWM(min_alt, max_alt, nsteps,final_lat,final_lon,annee,jours,secondes,o_file,o_folder);
% Remark: wind_proj not used in call to MSISEHWM since we will use the winds from ECMWF.
[ALTITUDE, DENSITY, TEMPERATURE, SOUNDSPEED, PRESSURE, LOCALPRESSURESCALE, G, NBVSQ, KAPPA, MU, MUVOL, NORTHWIND, EASTWIND, ~, CP, CV, GAMMA] = extract_atmos_model(o_file,3,0,0);

% Interpolate ECMWF on MSISE grid.
disp(['[',mfilename,'] Interpolating ECMWF quantities on MSISE (finer) grid.']);
T_e2m=spline(z_ecmwf,T_ecmwf,ALTITUDE);
p_e2m=spline(z_ecmwf,p_full_ecmwf,ALTITUDE);
g_e2m=spline(z_ecmwf,g_ecmwf,ALTITUDE);
w_M_e2m=spline(z_ecmwf,w_M_ecmwf,ALTITUDE);
W_Z_e2m=spline(z_ecmwf,w_Z_ecmwf,ALTITUDE);

% Compute missing quantities:
% - sound speed can be computed from T_ECMWF and gamma_MSISE (c=sqrt(gamma*R*T/M)=sqrt(gamma*P/rho)),
% - N^2 can be computed from gamma_MSISE, g_ECMWF, and T_ECMWF (N^2=(gamma-1)*g^2/(gamma*R*T)).
disp(['[',mfilename,'] Computing missing quantities (sound speed, Brunt-Väisälä frequency).']);
% Version 1, assuming the molar mass of air does not change and is equal to dry air molar mass. Relying on 1 MSISE quantity.
% soundspeed_merged=sqrt(GAMMA.*R.*T_e2m/M_dryair);
% Nsquared_merged=(GAMMA-1).*g_e2m.^2.*M_dryair./(GAMMA.*R.*T_e2m);
% Version 2, more exact. Relying on 2 MSISE quantities.
soundspeed_merged=sqrt(GAMMA.*p_e2m./DENSITY);
Nsquared_merged=(GAMMA-1).*g_e2m.^2.*DENSITY./(GAMMA.*p_e2m);
w_P=cos(wind_proj*pi/180.)*W_Z_e2m+sin(wind_proj*pi/180.)*w_M_e2m;

% Output newly-built file.
disp(['[',mfilename,'] Rewriting model to another file.']);
new_o_file=['msiseecmwf_',o_file];
rewrite_atmos_model(new_o_file, o_file, ALTITUDE, DENSITY, T_e2m, soundspeed_merged, p_e2m, LOCALPRESSURESCALE, g_e2m, Nsquared_merged, KAPPA, MU, MUVOL, w_M_e2m, W_Z_e2m, w_P, CP, CV, GAMMA);

kek=what(o_folder);
fullnewpath=[kek.path,filesep,new_o_file];
disp(['[',mfilename,'] Model stored in: ''',fullnewpath,'''.']);
