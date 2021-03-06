% Author:        Léo Martire.
% Description:   Calls ECMWF API with requested parameters.
% Notes:         Dates should be under format 2016-05-27.
%                Times should be under format HH:MM:SS[/HH:MM:SS[/HH:MM:SS[/...]]]
%                An API key is necessary to run the Python script calling the API. To obtain one, a demand must be forwarded to the ECMWF website.
%
% Usage:         1) Make sure the Python script is alongside this Matlab script.
%                2) Call this function with wanted parameters.
%                3) Execute provided command in a dedicated terminal.

function prepare_call_ECMWF_API(start_date, end_date, times, minmax_lon, minmax_lat, outputFileName)
  python_script = 'era5.py';
  apikeyfilename = '.ecmwfapirc';
  
  % API key data.
  % TODO: at some point, use a dedicated one instead of Guerman's.
  apiurl='https://api.ecmwf.int/v1';
  apikey=[];
  apimail=[];
  
  % Input treatment.
  times=char(times);
  outputFileName=char(outputFileName);
  
  % Save calling directory.
  calling_dir=pwd;
  
  % Get user HOME and check if API key is here.
  cd('~');
  if(not(exist(apikeyfilename,'file')))
    % If it does not exist, create it.
    disp(['  [',mfilename,'] API key file does not exist, creating it.']);
    apikey_fid = fopen(apikeyfilename, 'w');
    fprintf(apikey_fid, ['{\n    "url"   : "',apiurl,'",\n    "key"   : "',apikey,'",\n    "email" : "',apimail,'"\n}']);
    fclose('all');
  end
  
  % Get this function's folder and cd into it.
  this_func_path_splitted=split(mfilename('fullpath'),'/');
  this_func_path_splitted{end}='';
  this_func_dir=join(this_func_path_splitted,'/');
  this_func_dir=this_func_dir{1};
  clear('this_func_path_splitted');
  cd(this_func_dir);
  
  % Prepare call.
  area_str=[num2str(min(minmax_lat)),'/',num2str(min(minmax_lon)),'/',num2str(max(minmax_lat)),'/',num2str(max(minmax_lon))];
  args=[start_date, ' ', end_date, ' ', times, ' ', area_str, ' ', outputFileName];
  command=['python ',this_func_dir,python_script,' ',args];
  
  disp(['  [',mfilename,'] Open a dedicated terminal, and copy-paste this command to call the ECMWF API:']);
  disp(' ');
  disp(command);
  disp(' ');
  disp(['  [',mfilename,'] Then, wait until request is treated.']);

  % cd back to previous folder.
  cd(calling_dir);
end