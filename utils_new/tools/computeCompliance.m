% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function [CompZ, CompH, V] = computeCompliance(atmfile, parfile, speedType)
  if(not(exist('speedType')))
    speedType = '+w';
  end
  
  [~, ~, ~, SOUNDSPEED, ~, ~, ~, ~, ~, ~, ~, ~, ~, WIND] = extract_atmos_model(atmfile,3,0,-1);
  % get wind at ground
  SOUNDSPEED = SOUNDSPEED(1);
  WIND = WIND(1);
  
  model = readExampleFiles_extractParfileModels(parfile); % get parfile models;
  model = model(2, :); disp(['[',mfilename,'] From parfile, assume first ground model is numbered 2.']);
  % extract rho, vp, and vs from this model
  rho = model(3);
  vp = model(4);
  vs = model(5);
  
  switch(speedType)
    case '+w'
      signC = 0;
      signW = 1;
    case '-w'
      signC = 0;
      signW = -1;
    case 'c'
      signC = 1;
      signW = 0;
    case 'c+w'
      signC = 1;
      signW = 1;
    case 'c-w'
      signC = 1;
      signW = -1;
    otherwise
      error(['[',mfilename,', ERROR] speedType not implemented']);
  end
  V = signC*SOUNDSPEED + signW*WIND;
  
  [CompZ, CompH] = compliance(1, V, rho, vp, vs);
  
  CompZ = abs(CompZ);
  CompH = abs(CompH);
  disp(['[',mfilename,'] Computed compliances, with']);
  disp([blanks(length(mfilename)+2),'     V   = ',sprintf('%6.1f',V),' [m/s] (',speedType,' at ground),']);
  disp([blanks(length(mfilename)+2),'     rho = ',sprintf('%6.1f',rho),' [kg/m^3],']);
  disp([blanks(length(mfilename)+2),'     vp  = ',sprintf('%6.1f',vp),' [m/s],']);
  disp([blanks(length(mfilename)+2),'     vs  = ',sprintf('%6.1f',vs),' [m/s],']);
  disp([blanks(length(mfilename)+2),' are C_Z = ',sprintf('%12.3e',CompZ),' [??],']);
  disp([blanks(length(mfilename)+2),'     C_H = ',sprintf('%12.3e',CompH),' [??].']);
end

