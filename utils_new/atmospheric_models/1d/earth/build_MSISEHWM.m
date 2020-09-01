% Author:        Léo Martire.
% Description:   Builds a 1D atmospheric model file using MSISE and HWM models, and prints it in a format compatible with SPECFEM2D-DG.
% Notes:         Make sure the wrapper was compiled (see variable 'executable_fullpath' below).
%
% Usage:
%   [output_file] = build_MSISEHWM(altMin, altMax, nLayers, lat, lon, yyears, ddays, ssecs, outFile, outFold, wind_projection_azimuth, F107A, F107, AP)
% with:
%   altMin      [m] minimum altitude (should be >=0),
%   altMax      [m] maximum altitude,
%   nLayers     [1] number of layers,
%   lat         [deg] latitude,
%   lon         [deg] longitude,
%   yyears      [years] number of years since 2000,
%   ddays       [days] number of days since beginning of year,
%   ssecs       [s] number of seconds since beginning of day,
%   outFile     output file name,
%   outFold     output folder to which the output file will be moved,
%   projAz      [deg] projection azimuth for wind,
%   F107A       MSISE parameter (optional, defaults to 106.7160),
%   F107        MSISE parameter (optional, defaults to 131),
%   AP          MSISE parameter (optional, defaults to 37),
% yields:
%   output_file the path to the produced output file.

function [output_file] = build_MSISEHWM(altMin, altMax, nLayers, lat, lon, yyears, ddays, ssecs, outFile, outFold, projAz, F107A, F107, AP)
  % Path to executable.
  thisFolder = [regexprep(mfilename('fullpath'), mfilename, '')];
  executable_fullpath = [thisFolder, filesep, 'MSISE_HWM_wrapper/msisehwm'];
  
  % Input check.
  if(nargin<10)
    error(['  [', mfilename, 'ERROR] Not enough arguments.']);
  end
  % Set default values for optional arguments.
  if(~exist('F107A', 'var'))
    F107A=106.7160;
  end
  if(~exist('F107', 'var'))
    F107=131;
  end
  if(~exist('AP', 'var'))
    AP=37;
  end
  if(~exist('wind_projection_azimuth', 'var'))
    projAz = 0;
  end
  
  % Check existence of executable.
  if(~exist(executable_fullpath, 'file'))
    error(['[', mfilename, ' ERROR] Wrapper executable ''msisehwm'' not found. If not compiled, compile it beforehand. If compiled, make sure it is in ''./wrapper/msisehwm'' relative to this function.']);
  end
  
  % Build call.
  command = [executable_fullpath, ' ', num2str(altMin),    ' ', num2str(altMax),   ' ', num2str(nLayers), ' ', ...
             sprintf('%.5f', lat), ' ', sprintf('%.5f', lon), ' ', num2str(yyears),   ' ', num2str(ddays), ' ', ...
             sprintf('%.5f', ssecs), ' ', sprintf('%.5f', F107A), ' ', sprintf('%.5f', F107), ' ', ...
             sprintf('%.5f', AP), ' ', outFile,            ' ', sprintf('%.5f', projAz)];
  
  % Call.
  system(command);
  
  % Check existence of output folder, create it if not found.
  if(outFold(end)~='/')
    outFold = [outFold, '/'];
  end
  if (exist(outFold, 'dir') ~= 7)
    mkdir(outFold);
  end
  
  % Build stored filename.
  output_file = [outFold, filesep, outFile];
  
  % Move created file to output directory.
  system(['mv ', outFile, ' ', output_file]);
end

