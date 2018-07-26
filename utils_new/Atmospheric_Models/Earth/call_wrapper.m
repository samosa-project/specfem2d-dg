% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         Make sure path to executable is correct. Configure
%                request using the example request below.
% Notes:         N/A.

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametrisation.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_msisehwm = "./wrapper/msisehwm"; % Path to wrapper executable.
projection_angle = 0; % [°], [0, 360]. % Projection of horizontal winds (useful for 2D simulations). Leave 0 if you do not know what you're doing.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example request.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arctic polar circle: 66°33'47''N = 66.56306°.
% Tropic of Cancer: 23°26′13.0''N =  23.43694°.
% Tropic of Capricorn: 23°26′13.0''S =  -23.43694°.
% Antarctic polar circle: 66°33'47''S = -66.56306°.
% Nothern summer solstice (21th of June) of 2000 is the 173rd day of the year.
% Nothern winter solstice (21th of December) of 2000 is the 356th day of the year.
foldername = 'stratospheric';
altitude_min = 0;
altitude_max = 150e3;
nsteps = 1501;
lats =        [66.56306, 45, 23.43694, 0, -23.43694, -66.56306]; % [°], [-90, 90]. Must be paired with longitudes' vector.
lons =        [ 0,        0,  0,       0,    0,        0      ]; % [°], [0, 360]. Must be paired with latitudes' vector.
year = 0;
days =        [356, 356, 173, 173]; % [ ], [0, 366]. Must be paired with seconds' vector.
secs = 3600 * [0,    12,   0,  12]; localtime=0; % Giving UT. [s], [0, 86400]. Must be paired with days' vector.
% LTH = [0, 12, 0, 12]; localtime=1; % Giving local time. [s], [0, 86400]. Must be paired with days' vector.
F107As = 106.7160 * ones(size(days)); % 81-day average of F10.7 radio flux (centered on the day of the event). Must be paired with days' vector.
F107s = 131 * ones(size(days)); % Daily F10.7 flux (for previous day). Must be paired with days' vector.
APs = 37 * ones(size(days)); % Daily magnetic index. Must be paired with days' vector.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin treatment,            %
% extraction, and output.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produce subfolder.
folder = ['./msishwm_', foldername, '/'];
if (exist(folder, 'dir') ~= 7)
  system(strcat("mkdir ", folder));
end

disp('MSISE-HWM models'' extraction started.');
for i = 1:numel(lats)
  lat = lats(i);
  lon = lons(i);

  for j = 1:numel(days)
    day = days(j);
    if (localtime == 0)
      % If UT is given, use it.
      sec = secs(j);
    else
      % If local time is given, convert it to UT.
      sec = mod(LTH(j) - lon / 15, 24) * 3600;
    end
    F107A = F107As(j);
    F107 = F107s(j);
    AP = APs(j);
 
    daystr = pad(num2str(day), 3, 'left', "0"); % For output file name formatting.
    secstr = pad(sprintf("%.1f",sec), 7, 'left', "0"); % For output file name formatting.
    latstr = strcat(pad(sprintf("%.5f", sign(lat)*lat), 8, 'left', "0")); % For output file name formatting.
    if (sign(lat) > 0)
      latstr = strcat("+", latstr);
    else
      latstr = strcat("-", latstr);
    end
    lonstr = pad(sprintf("%.5f", 270), 9, 'left', "0"); % For output file name formatting.
    output_file = strcat(num2str(2000 + year), "_", daystr, "_", secstr, "_", latstr, "_", lonstr, "_", num2str(altitude_min), "_", num2str(altitude_max), "_", num2str(nsteps), "_", sprintf("%.1f",projection_angle));
    command = strcat(path_to_msisehwm,     " ", num2str(altitude_min),  " ", num2str(altitude_max), " ", num2str(nsteps), " ", ...
                     sprintf("%.5f", lat), " ", sprintf("%.5f", lon),   " ", num2str(year),         " ", num2str(day), " ", ...
                     sprintf("%.5f", sec), " ", sprintf("%.5f", F107A), " ", sprintf("%.5f", F107), " ", ...
                     sprintf("%.5f", AP),  " ", output_file,            " ", sprintf("%.5f", projection_angle));
%     command % Debug print.
    ie = system(command);
    system(strcat("mv ", output_file, " ", folder, output_file));
    output_file = strcat(folder, output_file);
  end
end

disp(' ');
disp(['MSISE-HWM models'' extraction finished. See files in ''',folder,'''.']);