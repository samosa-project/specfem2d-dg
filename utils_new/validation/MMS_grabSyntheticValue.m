function [Xq, Yq, value] = MMS_grabSyntheticValue(OFd)
  debug = 0;
  
  OFd=char(OFd);
  if(not(OFd(end)==filesep))
    OFd = [OFd, filesep];
  end
  
  [Xq, Yq, ~, ~] = loadStations(OFd);
  
  station = 1;
  unknown = 'BXZ';
  typeDisplay = readExampleFiles_extractParam([OFd, 'input_parfile'],'seismotype','int');
  if(not(typeDisplay==2))
    error('kek');
  end
  
  [Ztime, Zamp, ~, ~, ~] = gji2020_loadSomeSynthetics(OFd, station, typeDisplay, unknown, 1, 0, 1, 0, -1);
%   value = Zamp(end);
  value = mean(Zamp(end-999:end)); % more stable value, since oscillates alot
  
  if(debug)
    figure();
    plot(Zamp);
  end
end

