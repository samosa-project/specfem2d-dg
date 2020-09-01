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

function [xminmax, zminmax, interface, Xsource, debfin, d, name] = AboveGroundSimple_Microbaroms(simulationfolder)
  addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools'); 
  parfile        = [simulationfolder, 'parfile_input'];
  sourcefile     = [simulationfolder, 'source_input'];
  interfacesfile = [simulationfolder, 'interfaces_input'];
  [xminmax, zminmax, Xsource] = readExampleFiles(parfile, sourcefile, interfacesfile);
  STAT_DX = 15e3;
  STAT_Z = [1, 15, 30]*1e3; % altitude of stations.
  groundclearance = 3;
  interface = [-1e9, 1e9;0, 0];
  isitok=-1;
  while(not(ismember(isitok,[0,1])))
    isitok=input(['[', mfilename, '] Stations'' DX planned to be ',num2str(STAT_DX),' [m]. Is that ok (0 for no, 1 for yes)? > ']);
  end
  if(isitok==0)
    error(['[',mfilename,', ERROR] Parameter was not ok, re-choose parametrisation in script.']);
  end
  isitok=-1;
  while(not(ismember(isitok,[0,1])))
    isitok=input(['[', mfilename, '] Stations'' Z planned to be ',num2str(STAT_Z),' [m]. Is that ok (0 for no, 1 for yes)? > ']);
  end
  if(isitok==0)
    error(['[',mfilename,', ERROR] Parameter was not ok, re-choose parametrisation in script.']);
  end
  lid = 0;
%   lid = lid+1;
%   d(lid) = 0;
%   debfin(lid, 1, :) = Xsource(1)*[1, 1]; % xdeb xfin
%   debfin(lid, 2, :) = Xsource(2)*[1, 1]; % zdeb zfin
%   name{lid} = ['source'];
  for i=1:numel(STAT_Z)
    lid = lid+1;
    d(lid) = STAT_DX;
    debfin(lid, 1, :) = [xminmax(1)+d(lid), Xsource(1)]; % xdeb xfin
    debfin(lid, 2, :) = [1, 1]*STAT_Z(i); % zdeb zfin
    name{lid} = ['horizontal altitude ',num2str(STAT_Z),' [m]'];
  end
  lid = lid+1;
  d(lid) = 1000;
  debfin(lid, 1, :) = [-1, 1]*8e3; % xdeb xfin
  debfin(lid, 2, :) = [1, 1]*groundclearance; % zdeb zfin
  name{lid} = ['microbarom monitoring'];
end