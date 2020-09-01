% Author:        Léo Martire.
% Description:   Identifies the folder structure, addpath-genpath all utils, and returns the full path to the EXAMPLES folder.
% Notes:         N. A.
%
% Usage:
%   [SPCFMEXloc] = setup_overall()
% with:
%   N. A.
% yields:
%   SPCFMEXloc the full path to the EXAMPLES folder.

function [SPCFMEXloc] = setup_overall()
  % Identify the folder containing this script.
  thisFolder = [regexprep(mfilename('fullpath'),mfilename,'')];
  
  % Try and go up to root of SPECFEM2D-DG, defined as the folder containing the '.git' folder.
  curFolder = thisFolder;
  found = 0;
  while(not(found))
    if(exist([curFolder,filesep,'.git'],'dir'))
      found = 1;
      break;
    end
    spl = split(curFolder, filesep);
    spl(end) = [];
    curFolder = join(spl, filesep);
    curFolder = curFolder{1};
  end
  rootFolder = curFolder;
  
  addpath(genpath([rootFolder,'utils']));
  
  SPCFMEXloc = [rootFolder,filesep,'EXAMPLES',filesep];
end

