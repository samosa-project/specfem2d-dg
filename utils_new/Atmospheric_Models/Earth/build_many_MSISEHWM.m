% Author:        Léo Martire.
% Description:   Builds many atmospheric models files (directly compatible with SPECFEM) using MSISE and HWM models.
% Notes:         N/A.
%
% Usage:         Configure request using the example request below.

clear all;
close all;
clc;

projection_azimuth = 0; % [°], [0, 360]. % Projection azimuth for 2D horizontal winds.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example request.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arctic polar circle: 66°33'47''N = 66.56306°.
% Tropic of Cancer: 23°26′13.0''N =  23.43694°.
% Tropic of Capricorn: 23°26′13.0''S =  -23.43694°.
% Antarctic polar circle: 66°33'47''S = -66.56306°.
% Nothern summer solstice (21th of June) of 2000 is the 173rd day of the year.
% Nothern winter solstice (21th of December) of 2000 is the 356th day of the year.
output_foldername = 'stratospheric';
altitude_min = 0;
altitude_max = 150e3;
nsteps = 1501;
lats  =        [66.56306, 45, 23.43694, 0, -23.43694, -66.56306]; % [°], [-90, 90]. Size must be paired with longitudes' vector.
lons  =        [ 0,        0,  0,       0,   0,         0      ]; % [°], [0, 360]. Size must be paired with latitudes' vector.
annee = 0;
days  =        [356, 356, 173, 173]; % [ ], [0, 366]. Size must be paired with seconds' vector.
secs  = 3600 * [0,    12,   0,  12]; SAT1_UT0=0; % Giving UT. [s], [0, 86400]. Size must be paired with days' vector.
% LTH = [0, 12, 0, 12]; SAT1_UT0=1; % Giving solar apparent time. [hrs], [0, 24]. Size must be paired with days' vector.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin treatment,            %
% extraction, and output.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produce subfolder.
output_folder = ['./msishwm_', output_foldername, '/'];
if (exist(output_folder, 'dir') ~= 7)
  mkdir(output_folder);
end

disp(['[',mfilename,'] MSISE-HWM models'' extraction started.']);
for i = 1:numel(lats)
  lat = lats(i);
  lon = lons(i);
  for j = 1:numel(days)
    jour = days(j);
    if(SAT1_UT0)
      [~, ~, ~, sec] = UT2SAT(LTH(j)*3600, lon, 1); % If solar apparent time is given, convert it to UT.
    else
      sec = secs(j); % If UT is given, use it.
    end
    output_file = model_namer(altitude_min, altitude_max, nsteps, lat, lon, annee, jour, sec, projection_azimuth);
    o_f_stored = build_MSISEHWM(altitude_min, altitude_max, nsteps, lat, lon, annee, jour, sec, output_file, output_folder, projection_azimuth);
  end
end

disp(' ');
disp(['[',mfilename,'] MSISE-HWM models'' extraction finished. See files in ''',output_folder,'''.']);