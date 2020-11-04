function [Z] = computeImpedanceAir(atmfile)
  [~,rho,~,c]=extract_atmos_model(atmfile,3,0,-1);
  
  % get values at ground
  rho=rho(1);
  c=c(1);
  Z = impedance(rho, c);
  
  disp(['[',mfilename,'] Computed air impedance, with']);
  disp([blanks(length(mfilename)+2),'    rho = ',sprintf('%12.3e',rho),' [kg/m^3],']);
  disp([blanks(length(mfilename)+2),'    c   = ',sprintf('%6.1f',c),' [m/s],']);
  disp([blanks(length(mfilename)+2),' is Z   = ',sprintf('%12.3e',Z),' [s.m^2/kg].']);
end