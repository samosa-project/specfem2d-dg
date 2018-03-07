% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [datestr, posstr, year, daysincenewyear, secondssincenewday, lat, lon, f107a, f107, ap] = extract_data_setup(DATAFILE)
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
      secondssincenewday=str2num(regexprep(regexprep(line, 'lat.*', ''), 'year *[0-9]+ day *[0-9]+ seconds', ''));
      datestr=[num2str(daysincenewyear), 'th day of ', num2str(year), ' at ' num2str(floor(secondssincenewday/3600)),':',num2str(floor((secondssincenewday - floor(secondssincenewday/3600)*3600)/60)), ' UT'];

      lat=str2num(regexprep(regexprep(line, 'year *[0-9]+ day *[0-9]+ seconds *[0-9]+\.[0-9]+ *lat', ''), 'lon.*', ''));
      lon=str2num(regexprep(line, 'year *[0-9]+ day *[0-9]+ seconds *[0-9]+\.[0-9]+ *lat *-?[0-9]+\.[0-9]+ *lon', ''));
      posstr=['lat. ',num2str(lat), ', lon. ',num2str(lon)];
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

