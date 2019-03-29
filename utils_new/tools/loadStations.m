% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   [x_stat, z_stat, y_stat, stat_file] = loadStations(OFdir)
% with:
%   TODO.
% yields:
%   TODO.

function [x_stat, z_stat, y_stat, stat_file] = loadStations(OFdir)
  try
    stat_file = importdata([OFdir, 'STATIONS']);
  catch
    disp(['[',mfilename,'] STATIONS file not found in OUTPUT_FILES directory.']);
    try
%       stat_file = importdata([rootDir,'/DATA/STATIONS']);
      stat_file = importdata([OFdir,'../DATA/STATIONS']);
      disp(['[',mfilename,'] STATIONS file found in root directory (OUTPUT_FILES*/../DATA/ folder).']);
    catch
      error(['[',mfilename,', ERROR] Cannot find STATIONS file.']);
    end
  end
  pos_stations = [stat_file.data(:, 1) stat_file.data(:, 3) stat_file.data(:, 2)];
  x_stat = pos_stations(:, 1);
  z_stat = pos_stations(:, 3);
  y_stat = pos_stations(:, 2);
end