% Author:        Léo Martire.
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

%function [datestring, posstr, secondaryinfo, year, daysincenewyear, secondssincenewday, lat, lon, f107a, f107, ap] = extract_atmos_model_setup(DATAFILE)
function [dateString, positionString, secondaryInfoString, info] = extract_atmos_model_setup(DATAFILE)
  planettt = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'};
  
  fid = fopen(DATAFILE);
  if(fid==-1)
    error(strcat("Cannot open file ", DATAFILE,').'))
  end
  
  dateString='';
  positionString='';
  secondaryInfoString=[];
  year=-1;
  daysincenewyear=-1;
  secondssincenewday=-1;
  lat=-1;
  lon=-1;
  f107a=-1;
  f107=-1;
  ap=-1;
  info=struct();
  
  expr_float = '[0-9]+\.[0-9]+';
  expr_float_possibly_neg = '-?[0-9]+\.[0-9]+';
  expr_int = '[0-9]+';
  
  stop=0;
  count=0;
  while(stop==0)
    count=count+1;
    line = fgetl(fid);
    if(count==1)
      % try to find planet
      PLA_expr=['^ *PLANET *= *',expr_int];
      PLA = regexp(line, PLA_expr, 'match'); % match form
      if(isempty(PLA))
        PLA = 3; % earth by default
        disp(['[',mfilename,'] Planet string not found, assuming earth.']);
      else
        PLA = regexp(PLA{1}, expr_int, 'match'); % extract float from form
        PLA = str2num(PLA{1}); % convert to int
      end
      info.planet=PLA;
      
      % remove planet string and continue
      line = regexprep(line, PLA_expr, '');
      
      switch(PLA)
        case 3
          % Earth, format is "year 0 day 0 seconds 54000.0000 lat 4.502384 lon 135.623447"
          % try to find year
          year = regExpQuantity(line, 'year', expr_int); % convert to float
          if(isnan(year))
            year=2000; % safeguard if not found
          end
          info.year=year;
%           year = str2num(regexprep(regexprep(line, 'day.*', ''), 'year', ''));
%           if(isempty(year))
%             disp(['[',mfilename,'] Regexp on ''year'' has found nothing, assuming header is empty.']);
%             break;
%           end
          % try to find day since new year
          daysincenewyear = str2num(regexprep(regexprep(line, 'seconds.*', ''), 'year *[0-9]+ *day', ''));
          if(isempty(daysincenewyear))
            daysincenewyear=1; % safeguard if not found
          end
          info.doy=daysincenewyear;
        case 4
          % mars 
          LS = regExpQuantity(line, 'LS', expr_float); % try to find LS
          info.ls=LS;
          LT = regExpQuantity(line, 'LT', expr_float); % try to find LT
          info.lt=LT;
        otherwise
          error(['should not have reached this part of code'])
      end
      % always try to find lat and lon
%       lat = str2num(regexprep(regexprep(line, 'year *[0-9]+ *day *[0-9]+ *seconds *[0-9]+\.[0-9]+ *lat', ''), 'lon.*', ''))
%       lat = str2num(regexprep(regexprep(line, ['year *[0-9]+ *day *[0-9]+ *seconds *',expr_float,'+ *lat'], ''), 'lon.*', ''))
%       lat = regexp(line, ['lat ',expr_float], 'match'); % match form
%       lat = regexp(lat{1}, expr_float, 'match'); % extract float from form
      lat = regExpQuantity(line, 'lat', expr_float_possibly_neg); % convert to float
      if(isnan(lat))
        lat=0; % safeguard if not found
      end
      info.lat=lat;
%       lon = str2num(regexprep(line, 'year *[0-9]+ *day *[0-9]+ *seconds *[0-9]+\.[0-9]+ *lat *-?[0-9]+\.[0-9]+ *lon', ''))
%       lon = regexprep(line, ['year *[0-9]+ *day *[0-9]+ *seconds *[0-9]+\.[0-9]+ *lat *-?[0-9]+\.[0-9]+ *lon'], '');
%       lon = regexp(line, ['lon ',expr_float], 'match'); % match form
%       lon = regexp(lon{1}, expr_float, 'match'); % extract float from form
%       lon = str2num(lon{1}); % convert to float
      lon = regExpQuantity(line, 'lon', expr_float_possibly_neg); % convert to float
      if(isnan(lon))
        lon=0; % safeguard if not found
      end
%       lon = str2num(regexprep(line, ['year *[0-9]+ *day *[0-9]+ *seconds *[0-9]+\.[0-9]+ *lat *-?[0-9]+\.[0-9]+ *lon'], ''))
      info.lon=lon;
      
      % prepare lon lat str human readable
      if(lon<=180)
        lonstr=[sprintf('%1.1f',lon),'^\circ\mathrm{E}'];
      else
        lonstr=[sprintf('%1.1f',360-lon),'^\circ\mathrm{W}'];
      end
      if(lat>=0)
        latstr=[sprintf('%1.1f',lat),'^\circ\mathrm{N}'];
      else
        latstr=[sprintf('%1.1f',-lat),'^\circ\mathrm{S}'];
      end
      
      % Build position string.
%       posstr=['(lat., lon.) = (',latstr,', ',lonstr,')'];
      positionString=[planettt{PLA},'$\left(',latstr,', ',lonstr,'\right)$'];
      
      % Build date string.
      switch(PLA)
        case 3
          % Earth, easy peasy lemon squeezy% Deduce month and day.
          prevsum=0;
          sum=0;
          m=1;
          while(sum<daysincenewyear)
            prevsum=sum;
            sum=sum+eomday(year,m);
            m=m+1;
          end
          m=m-1;
          info.month=m;
          dom=daysincenewyear-prevsum;
          info.dom=dom;
          yyyymmdd=[num2str(year),'/',pad(num2str(m),2,'left','0'),'/',pad(num2str(dom),2,'left','0')];
          % Build hh:mm:ss.
          secondssincenewday=str2num(regexprep(regexprep(line, 'lat.*', ''), 'year *[0-9]+ *day *[0-9]+ *seconds', ''));
          if(isempty(secondssincenewday))
            secondssincenewday=0; % safeguard if not found
          end
          info.sod=secondssincenewday;
    %       hhmm=[pad(num2str(floor(secondssincenewday/3600)),2,'left','0'),':',pad(num2str(floor((secondssincenewday - floor(secondssincenewday/3600)*3600)/60)),2,'left','0')];
    %       hhmmss=[hhmm,':','00'];
          hhmmss=datestr(seconds(secondssincenewday),'hh:MM:SS');
          % Convert to local time using longitude.
          CoordinatedUniversalTimeStr=[yyyymmdd,' ',hhmmss];
          localsat=UTC2SolarApparentTime(CoordinatedUniversalTimeStr,lon); % See https://fr.mathworks.com/matlabcentral/fileexchange/32804-convert-utc-to-solar-apparent-time.
          splitlocalsat=split(localsat, ' ');
          % Build date string.
    %       datestr=[num2str(daysincenewyear), '$^{th}$ of ', num2str(year), ', ',num2str(floor(secondssincenewday/3600)),':',num2str(floor((secondssincenewday - floor(secondssincenewday/3600)*3600)/60)), ' UT'];
          dateString=[yyyymmdd, ', ', hhmmss, ' UTC, ', splitlocalsat{2}, ' LT'];
        case 4
          % mars case
          dateString = ['LS ',sprintf('%.3f',LS), '$^\circ$, LT ',sprintf('%.3f',LT), 'h'];
        otherwise
          dateString = ['[unknown planet, cannot dedude datestring]'];
      end
      
    % Build secondary info
    elseif(count==2)
      switch(PLA)
        case 3
          f107a=str2num(regexprep(regexprep(line, '.*F107A', ''), 'F107.*', ''));
          info.f107a=f107a;
          f107=str2num(regexprep(regexprep(line, '.*F107', ''), 'AP.*', ''));
          info.f107=f107;
          aptab=str2num(regexprep(line, '.*AP', ''));
          if(isempty(aptab))
            aptab=0; % safeguard if not found
          end
          ap=aptab(1);
          info.ap=ap;
          secondaryInfoString=['F10.7 avg. = ', sprintf('%.1f',f107a), ', F10.7 = ', sprintf('%.1f',f107), ', AP = ', sprintf('%.1f',ap)];
        case 4
          % nothing
          secondaryInfoString=[];
        otherwise
          error('should not have reached here');
      end
    else
      stop=1;
    end
  end
  fclose('all');
end

function val = regExpQuantity(chararray, Q, expr)
  reggg = [Q,' +',expr];
  val = regexp(chararray, reggg, 'match'); % match form
  if(isempty(val))
    disp(['[',mfilename,'] regexp not found (expr = ''',reggg,''', string = ''',chararray,''')']);
    val = nan;
  else
    val = regexp(val{1}, expr, 'match'); % extract float from form
    val = str2num(val{1}); % convert to float
  end
end