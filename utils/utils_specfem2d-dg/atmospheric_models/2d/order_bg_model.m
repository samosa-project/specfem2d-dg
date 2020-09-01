% Author:        Léo Martire.
% Description:   Defines how quantities are ordered in the 2D background model files used by the LNS module.
% Notes:         Works both for writing files (to be read by SPECFEM2D-DG) and reading files (written by SPECFEM2D-DG).
%                This order should always agree with the relevant routines in 'src/specfem2D/lns_load_background_model.F90':
%                - 'lns_read_background_model' for the reading, and
%                - 'output_lns_interpolated_model' and 'output_lns_interpolated_model_binary' for the outputting.
%
% Usage:
%   [order, lab, tex, unit] = order_bg_model()
% with:
%   N. A.
% yields:
%   order a char matrix encoding the ordering of physical quantities,
%   lab   the corresponding human-readable label,
%   tex   the corresponding TeX-readable label,
%   unit  the corresponding unit.

function [order, lab, tex, unit] = order_bg_model()
  i = 1;
  order(i,:) = 'xx'; unit{i}='m'; lab{i} = 'x[m]'; tex{i}='x'; i=i+1;
  order(i,:) = 'zz'; unit{i}='m'; lab{i} = 'z[m]'; tex{i}='z'; i=i+1;
  order(i,:) = 'rh'; unit{i}='kg/m$^3$'; lab{i} = 'rho0[kg.m^{-3}]'; tex{i}='\rho_0'; i=i+1; % density
  order(i,:) = 'vx'; unit{i}='m/s'; lab{i} = 'v0x[m.s^{-1}]'; tex{i}='v_{0,x}'; i=i+1; % x wind
  order(i,:) = 'vz'; unit{i}='m/s'; lab{i} = 'v0z[m.s^{-1}]'; tex{i}='v_{0,z}'; i=i+1; % z wind
  order(i,:) = 'pr'; unit{i}='Pa'; lab{i} = 'p0[Pa]'; tex{i}='p_0'; i=i+1; % pressure
  order(i,:) = 'gr'; unit{i}='m/s$^2$'; lab{i} = 'g[m.s^{-2}]'; tex{i}='g'; i=i+1; % gravity
  order(i,:) = 'ga'; unit{i}='1'; lab{i} = 'gamma[1]'; tex{i}='\gamma'; i=i+1; % gamma
  order(i,:) = 'mu'; unit{i}='kg/(m$\cdot$s)'; lab{i} = 'mu[kg.m^{-1}.s^{-1}]'; tex{i}='\mu'; i=i+1; % viscosity
  order(i,:) = 'ka'; unit{i}='m$\cdot$kg/(s$^{3}\cdot$K)'; lab{i} = 'kappa[m.kg.s^{-3}.K{-1}]'; tex{i}='\kappa'; i=i+1; % thermal conductivity
end

