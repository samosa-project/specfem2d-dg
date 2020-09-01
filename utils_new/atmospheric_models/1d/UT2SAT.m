% Author:        Léo Martire.
% Description:   Computes local solar apparent time (SAT, SAT) from UT, or vice versa.
%                Outputs a date string by default.
%                Can output years since 2000, days since beginning of year, and seconds since beginning of day. Useful for calls to MSISE and HWM.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         Needs 'UTC2SolarApparentTime' (download from https://fr.mathworks.com/matlabcentral/fileexchange/32804-convert-utc-to-solar-apparent-time).
%
% t_in  input time (formatted as 'YYYY/MM/DD [[H]H:[M]M:[S]S]' or as '[H]H:[M]M:[S]S', or given as a double)
% lon   aimed longitude (in (0,360)°E)
% swi   switch (optional, 0 for UT to SAT, 1 for SAT to UT)
%
% t_out  output time, format 'yyyy/mm/dd HH:MM:SS'
% ys2000 years since 2000
% dsboy  days since beginning of year
% ssbod  seconds since beginning of day

function [t_out,ys2000,dsboy,ssbod]=UT2SAT(t_in,lon,swi)
  
  % Input check.
  if(nargin<2)
    error(['  [',mfilename,'ERROR] Not enough arguments.']);
  end
  
  if(not(exist('swi','var')))
    swi=-1;
    while(not(ismember(swi,[0,1])))
      swi=input('  UT to SAT (0) or SAT to UT (1)? > ');
    end
  end
  
  threshold_year=1900; % Threshold below which a date is considered not valid.
  
  if(ischar(t_in))
    if(length(t_in)<=8)
      t=datenum([date, ' ', t_in]);
      disp(['[',mfilename,'] Only time was given. Assuming date is today. Treating input as ',datestr(t),'.']);
    else
      t=t_in;
      disp(['[',mfilename,'] Treating input as ',datestr(t),'.']);
    end
  elseif(isnumeric(t_in))
    yyyy=year(datetime(datevec(t_in)))<1900;
    if(yyyy<threshold_year)
      t=datenum(date)+t_in/(24*3600);
      disp(['[',mfilename,'] Year (',num2str(yyyy),') is below threshold (',num2str(threshold_year),'). Assuming input is seconds since beginning of today. Treating input as ',datestr(t),'.']);
    else
      t=t_in;
      disp(['[',mfilename,'] Treating input as ',datestr(t),'.']);
    end
  end
  
  if(swi==0)
    % UT to SAT.
    disp(['[',mfilename,'] Converting UT to SAT at longitude ',num2str(lon),'.']);
    t_out = UTC2SolarApparentTime(datestr(t,'yyyy/mm/dd HH:MM:SS'),lon);
  else
    % SAT to UT.
    disp(['[',mfilename,'] Converting SAT at longitude ',num2str(lon),' to UT.']);
    t_out = UTC2SolarApparentTime(datestr(t,'yyyy/mm/dd HH:MM:SS'),-lon);
  end
  
  disp(['[',mfilename,'] Output assigned, format ''yyyy/mm/dd HH:MM:SS''.']);
  
  % Compute years since 2000, days since beginning of year, and seconds since beginning of day.
  ys2000=year(datetime(datevec(t_out)))-2000;
  dsboy=day(datetime(datevec(t_out)),'dayofyear');
  dd=split(datestr(datenum(t_out)));
  dd=dd{1};
  ssbod=(datenum(t_out)-datenum(dd))*24*3600;
end
