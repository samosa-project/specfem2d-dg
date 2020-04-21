function [RHO_cst, VX_cst, VZ_cst, E_cst, dRHO_x, dRHO_z, dVX_x, dVX_z, dVZ_x, dVZ_z, dE_x, dE_z, colourPlot, markerPlot, markerSize] = MMS_constants(testCase)
  % Parameters for analytic solution (cf. LNS_manufactured_solutions.mw).
  % Should be the same as in source code.
  switch(testCase)
    case 'inviscid'
      colourPlot = [1,0,0]; markerPlot = 'o'; markerSize = 10;
      RHO_cst = 0.001;
      VX_cst = 0.;
      E_cst = 0.05;
    case 'kappa'
      colourPlot = [0,0.4,0]; markerPlot = 'diamond'; markerSize = 8;
      RHO_cst = 0.;
      VX_cst = 0.;
      E_cst = 0.05;
    case 'mu'
      colourPlot = [0,0,1]; markerPlot = 'square'; markerSize = 10;
      RHO_cst = 0.;
      VX_cst = 0.001;
      E_cst = 0.;
    otherwise
      error(['[, ERROR]']);
  end
  dRHO_x = 1;
  dRHO_z = 2;
  dVX_x = 2;
  dVX_z = 0;
  VZ_cst = 0.;
  dVZ_x = 1;
  dVZ_z = 2;
  dE_x = 3;
  dE_z = 4;
end