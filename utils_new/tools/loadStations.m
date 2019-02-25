% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function [x_stat,z_stat,stat_file] = loadStations(rootDir, OFdir)
  try
    stat_file = importdata(strcat(OFdir, 'STATIONS'));
  catch
    disp(['[',mfilename,'] STATIONS file not found in OUTPUT_FILES directory.']);
    try
      stat_file = importdata(strcat(rootDir, '/DATA/STATIONS'));
      disp(['[',mfilename,'] STATIONS file found in root directory (OUTPUT_FILES*/../DATA/ folder).']);
    catch
      error(['[',mfilename,', ERROR] Cannot find STATIONS file.']);
    end
  end
  pos_stations = [stat_file.data(:, 1) stat_file.data(:, 2)];
  x_stat = pos_stations(:, 1);
  z_stat = pos_stations(:, 2);
end