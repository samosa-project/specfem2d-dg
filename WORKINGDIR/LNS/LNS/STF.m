function [STFvals] = STF(tt)
  perio = 0.5;
  a  = pi*pi*((1/perio).^2);
  to = 0.5;
%   dt = tt(2)-tt(1);
  STFvals = zeros(4,length(tt));

  STFvals = -2.d0*a*(tt-to).*exp(-a*(tt-to).^2); % Gaussian derivative in time.
  
%   is_backward=false;
%   if(is_backward)
%     data_exact = load(strcat('./sisvy_c340.dat'));
%     t_loc = data_exact(:,1);
%     rho_loc = data_exact(:,1);
%     vx_loc  = data_exact(:,2);
%     vy_loc  = data_exact(:,3);
%     p_loc   = data_exact(:,4);
% 
%     data_num   = load(strcat('./sisvy_c140.dat'));
%     t_loc_num = data_num(:,1);
%     rho_loc_num = data_num(:,1);
%     vx_loc_num  = data_num(:,2);
%     vy_loc_num  = data_num(:,3);
%     p_loc_num   = data_num(:,4);
%     for i=1:length(tt)
%       loc = find(abs(t_loc-tt(i)) <= dt/2);
%       loc = loc(1);
% 
%       rho=1;
%       c = 340;
%       p     = (c.^2).*rho./gamma_atmos(1);
%       E1     = p./(gamma_atmos(1)-1.);
%       c = 140;
%       p     = (c.^2).*rho./gamma_atmos(1);
%       E2     = p./(gamma_atmos(1)-1.);
% 
%       source_term(1,loc) = rho_loc(loc)-rho_loc_num(loc);
%       source_term(2,loc) = vx_loc(loc)-vx_loc_num(loc);
%       source_term(3,loc) = vy_loc(loc)-vy_loc_num(loc);
%       source_term(4,loc) = (p_loc(loc)-E1)-(p_loc_num(loc)-E2);
%     end
%   else
%     source_term = -2.d0*a*(tt-to).*exp(-a*(tt-to).^2); % Gaussian in time.
%   end
end