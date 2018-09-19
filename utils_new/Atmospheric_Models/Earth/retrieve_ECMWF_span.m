% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Retrieves the area and time spanned by a previously retrieved ERA5 netcdf file.
% Last modified: See file metadata.
% Usage:         1) Retrieve the ERA5 file's points (using the script 'retrieve_ECMWF').
%                2) Call this function using the points as arguments.
% Notes:         N/A.

function [lonSpan,latSpan,timeSpan] = retrieve_ECMWF_span(ecmwfPts)
  ecmwf_lon=ecmwfPts{1};
  ecmwf_lat=ecmwfPts{2};
  ecmwf_time=ecmwfPts{4};
  lonSpan=['(',num2str(min(ecmwf_lon)),', ',num2str(max(ecmwf_lon)),')°E'];
  latSpan=['(',num2str(min(ecmwf_lat)),', ',num2str(max(ecmwf_lat)),')°N'];
  timeSpan=['(',datestr(min(ecmwf_time),'YYYY/mm/dd HH:MM:SS'),', ',datestr(max(ecmwf_time),'YYYY/mm/dd HH:MM:SS'),') UT'];
  disp(['  [',mfilename,'] ERA5 file spans: ',lonSpan,',']);
  disp(['  [',mfilename,']                  ',latSpan,',']);
  disp(['  [',mfilename,']                  ',timeSpan,'.']);
end