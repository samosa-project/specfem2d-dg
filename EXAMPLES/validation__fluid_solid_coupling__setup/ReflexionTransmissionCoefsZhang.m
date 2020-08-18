
function [R, T] = ReflexionTransmissionCoefsZhang(fts0_stf1, vp_1, Z_1, vp_2, vs_2, Z_2P, Z_2S, i_inc)
  % Coefficients from Zhang (2018).
  % Zhang, G., Hao, C., & Yao, C. (2018). Analytical study of the reflection and transmission coefficient of the submarine interface. Acta Geophysica, 66(4), 449–460. https://doi.org/10.1007/s11600-018-0153-y
  
  % Conversion from fluid P to solid S is strangely low...
  
  normalise = 0; % add the sqrt() factor to transmissions
  
  if(fts0_stf1)
    R = [0.4, 0.3];
    T = 0.2;
    disp(['[',mfilename,', WARNING] Not implemented!']);
  else
    % FLUID-TO-SOLID
    i1 = i_inc; % i_inc is i1 in this case
    i2 = snells(vp_1, vp_2, i1); % deduce i2
    j2 = snells(vp_1, vs_2, i1); % deduce j2
    
    G = 2*sin(j2)*sin(2*j2)*cos(i1)*(Z_2P*cos(j2) - Z_2S*cos(i2)) - Z_1*cos(i2) - Z_2P*cos(i1);
    
    R = 1 + 2*Z_1*cos(i2) / G;
    
    T(1) = -2*Z_1*cos(i1)*cos(2*j2) / G;
    
    T(2) = 4*Z_1*sin(j2)*cos(i1)*cos(i2) / G;
    
    if(normalise)
      T(1) = T(1) * sqrt(Z_2P*cos(i2)/(Z_1*cos(i1))); % normalisation factor which I don't trust
      T(2) = T(2) * sqrt(Z_2S*cos(j2)/(Z_1*cos(i1))); % normalisation factor which I don't trust
    end
  end
end