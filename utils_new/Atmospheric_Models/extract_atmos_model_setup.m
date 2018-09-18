% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [datestr, posstr, year, daysincenewyear, secondssincenewday, lat, lon, f107a, f107, ap] = extract_atmos_model_setup(DATAFILE)
  fid = fopen(DATAFILE);
  if(fid==-1)
    error(strcat("Cannot open file ", DATAFILE,').'))
  end
  stop=0;
  count=0;
  while(stop==0)
    count=count+1;
    line = fgetl(fid);
    if(count==1)
      year=str2num(regexprep(regexprep(line, 'day.*', ''), 'year', ''));
      daysincenewyear=str2num(regexprep(regexprep(line, 'seconds.*', ''), 'year *[0-9]+ day', ''));
      lat=str2num(regexprep(regexprep(line, 'year *[0-9]+ day *[0-9]+ seconds *[0-9]+\.[0-9]+ *lat', ''), 'lon.*', ''));
      lon=str2num(regexprep(line, 'year *[0-9]+ day *[0-9]+ seconds *[0-9]+\.[0-9]+ *lat *-?[0-9]+\.[0-9]+ *lon', ''));
      if(lon<=180)
        lonstr=[sprintf('%1.1f',lon),'$^\circ$E'];
      else
        lonstr=[sprintf('%1.1f',360-lon),'$^\circ$W'];
      end
      if(lat>=0)
        latstr=[sprintf('%1.1f',lat),'$^\circ$N'];
      else
        latstr=[sprintf('%1.1f',-lat),'$^\circ$S'];
      end
      
      % Build position string.
%       posstr=['(lat., lon.) = (',latstr,', ',lonstr,')'];
      posstr=['(',latstr,', ',lonstr,')'];
      
      % Deduce month and day.
      prevsum=0;
      sum=0;
      m=1;
      while(sum<daysincenewyear)
        prevsum=sum;
        sum=sum+eomday(year,m);
        m=m+1;
      end
      m=m-1;
      dom=daysincenewyear-prevsum;
      yyyymmdd=[num2str(year),'/',pad(num2str(m),2,'left','0'),'/',pad(num2str(dom),2,'left','0')];
      % Build hh:mm:ss.
      secondssincenewday=str2num(regexprep(regexprep(line, 'lat.*', ''), 'year *[0-9]+ day *[0-9]+ seconds', ''));
      hhmm=[pad(num2str(floor(secondssincenewday/3600)),2,'left','0'),':',pad(num2str(floor((secondssincenewday - floor(secondssincenewday/3600)*3600)/60)),2,'left','0')];
      hhmmss=[hhmm,':','00'];
      % Convert to local time using longitude.
      CoordinatedUniversalTimeStr=[yyyymmdd,' ',hhmmss];
      localsat=UTC2SolarApparentTime(CoordinatedUniversalTimeStr,lon); % See https://fr.mathworks.com/matlabcentral/fileexchange/32804-convert-utc-to-solar-apparent-time.
      splitlocalsat=split(localsat, ' ');
      % Build date string.
%       datestr=[num2str(daysincenewyear), '$^{th}$ of ', num2str(year), ', ',num2str(floor(secondssincenewday/3600)),':',num2str(floor((secondssincenewday - floor(secondssincenewday/3600)*3600)/60)), ' UT'];
      datestr=[yyyymmdd, ', ', hhmmss, ' UT, ', splitlocalsat{2}, ' LT'];
      
    elseif(count==2)
      f107a=str2num(regexprep(regexprep(line, '.*F107A', ''), 'F107.*', ''));
      f107=str2num(regexprep(regexprep(line, '.*F107', ''), 'AP.*', ''));
      aptab=str2num(regexprep(line, '.*AP', ''));
      ap=aptab(1);
    else
      stop=1;
    end
  end
  fclose('all');
end

