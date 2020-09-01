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

function [Zs,Zp] = computeImpedanceGround(parfile)
  rhovpvs = readExampleFiles_extractParfileModels(parfile); % get parfile models;
  rhovpvs = rhovpvs(2, :); disp(['[',mfilename,'] From parfile, assume first ground model is numbered 2.']);
  rhovpvs = rhovpvs(3:5); % extract rho, vp, and vs from this model
%   Zs = 1/prod(rhovpvs([1, 3])); % 1/(rho*vs)
%   Zp = 1/prod(rhovpvs([1, 2])); % 1/(rho*vp)
  Zs = impedance(rhovpvs(1), rhovpvs(3));
  Zp = impedance(rhovpvs(1), rhovpvs(2));
  
  disp(['[',mfilename,'] Computed ground impedances, with']);
  disp([blanks(length(mfilename)+2),'     rho = ',sprintf('%6.1f',rhovpvs(1)),' [kg/m^3],']);
  disp([blanks(length(mfilename)+2),'     vp  = ',sprintf('%6.1f',rhovpvs(2)),' [m/s],']);
  disp([blanks(length(mfilename)+2),'     vs  = ',sprintf('%6.1f',rhovpvs(3)),' [m/s],']);
  disp([blanks(length(mfilename)+2),' are Z_S = ',sprintf('%12.3e',Zs),' [s.m^2/kg],']);
  disp([blanks(length(mfilename)+2),'     Z_P = ',sprintf('%12.3e',Zp),' [s.m^2/kg].']);
end

