% Author:        Léo Martire.
% Description:   Builds an atmospheric model file (directly compatible with SPECFEM) using MSISE and HWM models.
% Last modified: See file metadata.
% Usage:         1) Make sure the wrapper was compiled.
%                2) Make sure the path to the wrapper's executable is correctly set (hardcoded below).
%                3) Call the function with the wanted parameters.
% Notes:         N/A.
%
% altMin  minimum altitude (should be >=0)
% altMax  maximum altitude
% nLayers number of layers
% lat      latitude
% lon      longitude
% yyears   number of years since 2000
% ddays    number of days since beginning of year
% ssecs    number of seconds since beginning of day
% outFile  output file path
% outFold  output folder (to which the output file will be moved, somewhat
%          redundant with previous parameter for single calls but useful for multiple calls)
% wind_projection_azimuth projection azimuth for wind in [deg] (optional, useful for SPECFEM2D simulations)
% F107A    TODO (optional)
% F107     TODO (optional)
% AP       TODO (optional)

function output_file_stored = build_MSISEHWM(altMin, altMax, nLayers, lat, lon, yyears, ddays, ssecs, outFile, outFold, wind_projection_azimuth, F107A, F107, AP)
  % Default paths.
  default_executable_path="MSISE_HWM_wrapper/msisehwm";
  
  % Find necessary components.
  % Get this function's folder.
  this_func_path_splitted=split(mfilename('fullpath'),'/'); this_func_path_splitted{end}=''; this_func_dir=join(this_func_path_splitted,'/'); clear('this_func_path_splitted');
  % Normally, this function is alongside this function, thus the full path is the following:
  executable_fullpath = strcat(this_func_dir,default_executable_path);
  
  % Input check.
  if(nargin<10)
    error(['  [',mfilename,'ERROR] Not enough arguments.']);
  end
  % Set default values for optional arguments.
  if(~exist('F107A','var'))
    F107A=106.7160;
  end
  if(~exist('F107','var'))
    F107=131;
  end
  if(~exist('AP','var'))
    AP=37;
  end
  if(~exist('wind_projection_azimuth','var'))
    wind_projection_azimuth = 0;
  end  
%   [F107A,F107,AP,wind_project_angle]
  
  % Check existence of executable.
  if(~exist(executable_fullpath,'file'))
    error(strcat(['  [',mfilename,' ERROR] Wrapper executable ''msisehwm'' not found. If not compiled, compile it beforehand. If compiled, make sure it is in ''./wrapper/msisehwm'' relative to this function.']));
  end
  
  % Build call.
  command = strcat(executable_fullpath,    " ", num2str(altMin),        " ", num2str(altMax),       " ", num2str(nLayers), " ", ...
                   sprintf("%.5f", lat),   " ", sprintf("%.5f", lon),   " ", num2str(yyears),       " ", num2str(ddays), " ", ...
                   sprintf("%.5f", ssecs), " ", sprintf("%.5f", F107A), " ", sprintf("%.5f", F107), " ", ...
                   sprintf("%.5f", AP),    " ", outFile,                " ", sprintf("%.5f", wind_projection_azimuth));
  
  % Call.
  system(command);
  
  % Check existence of output folder, create it if not found.
  if(outFold(end)~='/')
    outFold = [outFold, '/'];
  end
  if (exist(outFold, 'dir') ~= 7)
    system(strcat("mkdir ", outFold));
  end
  
  % Build stored filename.
  output_file_stored = strcat(outFold, outFile);
  
  % Move created file to output directory.
  system(strcat("mv ", outFile, " ", output_file_stored));
end

