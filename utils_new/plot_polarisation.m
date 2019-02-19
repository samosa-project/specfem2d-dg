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

function plot_polarisation(sig_x, sig_z, titlefig)
  if(not(exist('titlefig')))
    titlefig='';
  end
  
  figure();
  
  plot(sig_x, sig_z, 'white', 'linewidth', 1);
  daspect([1 1 1]);
  
  minmax=[min(sig_x), max(sig_x), ...
          min(sig_z), max(sig_z)];
  tenPow_x=10.^floor(log10(abs(minmax)));
  minmax_nextint=sign(minmax).*ceil(abs(minmax./tenPow_x)).*tenPow_x;
  
  axis(minmax_nextint);
  set(gca,'Color','k');
  set(gca,'GridColor','white');
  set(gca, 'TickLabelInterpreter', 'latex');
  set(gca,'TickDir','both');
  grid on;
  box on;
  xlabel(['$v_x$ [m/s]']);
  ylabel(['$v_z$ [m/s]']);
  title(titlefig);
end