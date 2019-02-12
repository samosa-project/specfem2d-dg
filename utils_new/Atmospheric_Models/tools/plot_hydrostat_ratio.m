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

function [] = plot_hydrostat_ratio(Z, RHO, P, G)
  D=differentiation_matrix(Z, 0);
  hydrostatic_ratio = (D * P) ./ (-RHO .* G);
  
  figure('units','normalized','outerposition',[0 0 0.5 1]);
  plot(hydrostatic_ratio, Z, ones(size(Z)), Z, 'k:');
  xlabel('$\partial_zp/(-\rho g)$');
  ylabel('$z$ [m]');
  title({['Hydrostatic Ratio'],['max. error: ',sprintf('%.2g',max(abs(hydrostatic_ratio-1))*100),'\%']});
  set(gca,'ticklabelinterpreter','latex');
  grid on;
  box on;
end

