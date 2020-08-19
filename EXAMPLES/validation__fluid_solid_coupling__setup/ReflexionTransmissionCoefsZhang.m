
function [R, T] = ReflexionTransmissionCoefsZhang(fts0_stf1, alpha__1, rho__1, alpha__2, beta__2, rho__2, i_inc)
  % Coefficients from Zhang (2018).
  % Zhang, G., Hao, C., & Yao, C. (2018). Analytical study of the reflection and transmission coefficient of the submarine interface. Acta Geophysica, 66(4), 449–460. https://doi.org/10.1007/s11600-018-0153-y
  
  % Conversion from fluid P to solid S is strangely low...
  
  normalise = 0; % add the sqrt() factor to transmissions
  
  Z_1 = alpha__1*rho__1;
  Z_2P = alpha__2*rho__2;
  Z_2S = beta__2*rho__2;
  
  if(fts0_stf1)
    i__2 = i_inc; % i_inc is i2 in this case
    i__1 = snells(alpha__2, alpha__1, i__2); % deduce i1
    j__2 = snells(alpha__2, beta__2, i__2); % deduce j2
    
    a = alpha__1*beta__2*rho__1*rho__2;
    p = sin(i__2)/alpha__2; % slowness 
    b = -alpha__2*beta__2*rho__2^2 + 4*p^2*alpha__2*beta__2^3*rho__2^2 - 4*p^2*cos(i__2)*cos(j__2)*beta__2^4*rho__2^2 - 4*p^4*alpha__2*beta__2^5*rho__2^2;
    c = 2*p^2*alpha__1*beta__2^3*rho__1*rho__2;
    d = 2*p*cos(i__2)*beta__2^2*rho__2;
    E = b*cos(i__1)-a*cos(i__2);
    F = p*alpha__2*(c/(rho__1*alpha__1) - rho__2*beta__2);
    
    T = (b - beta__2*rho__2^2*alpha__2) / E;
    R(1) = - (rho__1*(-beta__2*rho__2^2 + a - c)*alpha__1 - c*rho__2*alpha__2) / (E*rho__1*alpha__1);
    R(2) = (rho__1*alpha__1 - rho__2*alpha__2) * d / E;
    
    disp(['[',mfilename,', WARNING] Not implemented!']);
  else
    % FLUID-TO-SOLID
    i__1 = i_inc; % i_inc is i1 in this case
    i__2 = snells(alpha__1, alpha__2, i__1); % deduce i2
    j__2 = snells(alpha__1, beta__2, i__1); % deduce j2
    
    G = 2*sin(j__2)*sin(2*j__2)*cos(i__1)*(Z_2P*cos(j__2) - Z_2S*cos(i__2)) - Z_1*cos(i__2) - Z_2P*cos(i__1);
    
    R = 1 + 2*Z_1*cos(i__2) / G;
    
    T(1) = -2*Z_1*cos(i__1)*cos(2*j__2) / G;
    
    T(2) = 4*Z_1*sin(j__2)*cos(i__1)*cos(i__2) / G;
    
    if(normalise)
      T(1) = T(1) * sqrt(Z_2P*cos(i__2)/(Z_1*cos(i__1))); % normalisation factor which I don't trust
      T(2) = T(2) * sqrt(Z_2S*cos(j__2)/(Z_1*cos(i__1))); % normalisation factor which I don't trust
    end
  end
end