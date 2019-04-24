% Author:        Léo Martire.
% Description:   Prepares stations for a specfic type of run.
% Notes:         TODO.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function [xminmax, zminmax, interface, Xsource, debfin, d, name] = AboveAndBelowGround_Periodic_WithTilt__MarsInSIGHT(simulationfolder)
  addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools');
  parfile        = [simulationfolder, 'parfile_input'];
  sourcefile     = [simulationfolder, 'source_input'];
  interfacesfile = [simulationfolder, 'interfaces_input'];
%   xminmax = [extractParamFromInputFile(parfile, 'xmin', 'float'), extractParamFromInputFile(parfile, 'xmax', 'float')];
%   Xsource = [extractParamFromInputFile(sourcefile, 'xs', 'float'), extractParamFromInputFile(sourcefile, 'zs', 'float')];
%   [zmin, zmax] = extractZminZmaxFromInterfacesFile(interfacesfile);
%   zminmax = [zmin,zmax];
  [xminmax, zminmax, Xsource] = readExampleFiles(parfile, sourcefile, interfacesfile);
  
  xminmaxNEARSOURSSS=[-15, 15]*1e3;
  spacingstations(1) = 4000;
  spacingstations(2) = 400;
  interface = [-1e9, 1e9;0, 0];
  ground_clearance = 10; % altitude/depth of the ground stations.
  shift_for_tiltcomputation = 5; % horizontal shift for stations used for tilt computation
  
  isitok=-1;
  while(not(ismember(isitok,[0,1])))
    isitok=input(['[', mfilename, '] Stations over whole domain planned to be spaced by ',num2str(spacingstations(1)),' [m]. Is that ok (0 for no, 1 for yes)? > ']);
  end
  if(isitok==0)
    error(['[',mfilename,', ERROR] Parameter was not ok, re-choose parametrisation in script.']);
  end
  isitok=-1;
  while(not(ismember(isitok,[0,1])))
    isitok=input(['[', mfilename, '] Stations near source       planned to be spaced by ',num2str(spacingstations(2)),' [m]. Is that ok (0 for no, 1 for yes)? > ']);
  end
  if(isitok==0)
    error(['[',mfilename,', ERROR] Parameter was not ok, re-choose parametrisation in script.']);
  end
  
  isitok=-1;
  while(not(ismember(isitok,[0,1])))
    isitok=input(['[', mfilename, '] Ground clearance (above and below) planned to be ',num2str(ground_clearance),' [m]. Is that ok (0 for no, 1 for yes)? > ']);
  end
  if(isitok==0)
    error(['[',mfilename,', ERROR] Parameter was not ok, re-choose parametrisation in script.']);
  end
  disp(['[', mfilename, '] Horizontal shift for tilt computation planned to be ',num2str(shift_for_tiltcomputation),' [m].']);

  lid = 0;
  
  lid = lid+1; idhoriz = lid;
  d(lid) = 0;
  debfin(lid, 1, :) = Xsource(1)*[1, 1]; % xdeb xfin
  debfin(lid, 2, :) = Xsource(2)*[1, 1]; % zdeb zfin
  name{lid} = ['source'];
  
  for LINE=1:numel(spacingstations)
    switch(LINE)
      case 1
        % whole domain
        locminmax = xminmax; locname='[over whole domain]';
      case 2
        % near source
        locminmax = xminmaxNEARSOURSSS; locname='[near source]';
      otherwise
        error('KOK');
    end

    lid = lid+1; idhoriz = lid;
    d(lid) = spacingstations(LINE);
    debfin(lid, 1, :) = [locminmax(1)+d(lid), locminmax(2)-d(lid)]; % xdeb xfin
    debfin(lid, 2, :) = [1, 1]*ground_clearance; % zdeb zfin
    name{lid} = ['horizontal over ground, ', locname];
    lid = lid+1; idhorizvz = lid;
    d(lid) = d(idhoriz);
    debfin(lid, 1, :) = debfin(idhoriz, 1, :); % xdeb xfin same as horizontal over ground
    debfin(lid, 2, :) = -1*[1, 1]*ground_clearance; % zdeb zfin
    name{lid} = ['horizontal under ground, ', locname];
    lid = lid+1;
    d(lid) = d(idhoriz);
    debfin(lid, :, :) = debfin(idhorizvz, :, :); % all same as horizontal under ground
    debfin(lid, 1, :) = debfin(lid, 1, :)-shift_for_tiltcomputation; % but shifted to the left
    name{lid} = ['horizontal under ground shifted -', num2str(shift_for_tiltcomputation), 'm for tilt computation, ', locname];
    lid = lid+1;
    d(lid) = d(idhoriz);
    debfin(lid, :, :) = debfin(idhorizvz, :, :); % all same as horizontal under ground
    debfin(lid, 1, :) = debfin(lid, 1, :)+shift_for_tiltcomputation; % but shifted to the right
    name{lid} = ['horizontal under ground shifted ', num2str(shift_for_tiltcomputation), 'm for tilt computation, ', locname];
    
    if(0)
      % on periodic boundary condition
      lid = lid+1; d(lid) = 2*ground_clearance;
      debfin(lid, 1, :) = [xminmax(2), xminmax(2)]; debfin(lid, 2, :) = [1, -1]*ground_clearance;
      name{lid} = ['over & under ground @ right boundary'];
      lid = lid+1; d(lid) = diff(xminmax)-2*shift_for_tiltcomputation;
      debfin(lid, 1, :) = [xminmax(1)+shift_for_tiltcomputation,xminmax(2)-shift_for_tiltcomputation]; debfin(lid, 2, :) = -1*[1, 1]*ground_clearance;
      name{lid} = ['left&right of right boundary for tilt computation'];
    end
  end
end