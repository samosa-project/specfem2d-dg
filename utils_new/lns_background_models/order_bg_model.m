function [order, lab] = order_bg_model()
  i = 1;
  order(i,:) = 'xx'; lab{i} = 'x[m]'; i=i+1;
  order(i,:) = 'zz'; lab{i} = 'z[m]'; i=i+1;
  order(i,:) = 'rh'; lab{i} = 'rho0[kg.m^{-3}]'; i=i+1; % density
  order(i,:) = 'vx'; lab{i} = 'v0x[m.s^{-1}]'; i=i+1; % x wind
  order(i,:) = 'vz'; lab{i} = 'v0z[m.s^{-1}]'; i=i+1; % z wind
  order(i,:) = 'pr'; lab{i} = 'p0[Pa]'; i=i+1; % pressure
  order(i,:) = 'gr'; lab{i} = 'g[m.s^{-2}]'; i=i+1; % gravity
  order(i,:) = 'ga'; lab{i} = 'gamma[1]'; i=i+1; % gamma
  order(i,:) = 'mu'; lab{i} = 'mu[kg.m^{-1}.s^{-1}]'; i=i+1; % viscosity
  order(i,:) = 'ka'; lab{i} = 'kappa[m.kg.s^{-3}.K{-1}]'; i=i+1; % thermal conductivity
end

