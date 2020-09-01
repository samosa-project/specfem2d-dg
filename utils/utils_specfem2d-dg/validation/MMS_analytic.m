function [drho_th, dvx_th, dvz_th, dE_th, dp_th] = MMS_analytic(testCase, Xe, Ye, RHO0, V0_x, E0, P0, GAM)
  % Build analytic solution (cf. LNS_manufactured_solutions.mw).
  [RHO_cst, VX_cst, VZ_cst, E_cst, dRHO_x, dRHO_z, dVX_x, dVX_z, dVZ_x, dVZ_z, dE_x, dE_z] = MMS_constants(testCase);
  % Constitutive variables.
  drho_th = RHO_cst * (sin(dRHO_x*pi*Xe) + sin(dRHO_z*pi*Ye));
  dvx_th  = VX_cst  * (sin(dVX_x*pi*Xe)  + sin(dVX_z*pi*Ye));
  dvz_th  = VZ_cst  * (sin(dVZ_x*pi*Xe)  + sin(dVZ_z*pi*Ye));
  dE_th   = E_cst   * (sin(dE_x*pi*Xe)   + sin(dE_z*pi*Ye));
  % Pressure.
  %dp_th = (GAM-1)*(-0.5*V0_x^2*drho_th + dE_th); % Compute directly.
  dp_th = (GAM-1)*(E0+dE_th - 0.5*(RHO0+drho_th).*((V0_x+dvx_th).^2 + dvz_th.^2)) - P0;
end

