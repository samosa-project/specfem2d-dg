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

function [CompZ, CompH, V] = computeCompliance(atmfile, parfile, speedType, groundModelNumber)
  if(not(exist('speedType')))
    speedType = '+w';
  end
  
  if(not(exist('groundModelNumber')))
    groundModelNumber='all';
  end
  
  if(strcmp(groundModelNumber,'all'))
    allGroundModels0_oneGroundModel1 = 0;
  else
    allGroundModels0_oneGroundModel1 = 1;
  end
  
  [~, ~, ~, SOUNDSPEED, ~, ~, ~, ~, ~, ~, ~, ~, ~, WIND] = extract_atmos_model(atmfile,3,0,-1);
  % get wind at ground
  SOUNDSPEED = SOUNDSPEED(1);
  WIND = WIND(1);
  
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
  
  model = readExampleFiles_extractParfileModels(parfile); % get parfile models;
  if(allGroundModels0_oneGroundModel1)
    model = model(groundModelNumber, :); disp(['[',mfilename,'] From parfile, select ground model numbered ',num2str(groundModelNumber),'.']);
  end
  
  for m = 1:size(model,1)
    % loop through models
    % extract rho, vp, and vs from this model
    rho = model(m, 3);
    vp = model(m, 4);
    vs = model(m, 5);
    
    [CompZ(m), CompH(m)] = compliance(1, V, rho, vp, vs);

    CompZ(m) = abs(CompZ(m));
    CompH(m) = abs(CompH(m));
    disp(['[',mfilename,'] Computed compliances for model ',num2str(groundModelNumber),' (rho = ',sprintf('%.1e',rho),' [kg/m^3], vp = ',sprintf('%6.1f',vp),' [m/s], vs = ',sprintf('%6.1f',vs),' [m/s]), with']);
    disp([blanks(length(mfilename)+2),'     V   = ',sprintf('%6.1f',V),' [m/s] (',speedType,' at ground),']);
    disp([blanks(length(mfilename)+2),' are C_Z = ',sprintf('%12.3e',CompZ(m)),' [??],']);
    disp([blanks(length(mfilename)+2),'     C_H = ',sprintf('%12.3e',CompH(m)),' [??].']);
  end
end

