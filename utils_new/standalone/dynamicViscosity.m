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

function [mu] = dynamicViscosity(rho, T, p)
  rho = reshape(rho, numel(rho), 1);
  T = reshape(T, numel(T), 1);
  p = reshape(p, numel(p), 1);
  
  [~, k, dN2] = thermodynamicalConstants();
  
  mu = (2 * k / (3 * pi^(1.5) * dN2^2)) .* T .* (rho ./ p).^(0.5);
  
  disp(['[',mfilename,'] Dynamic viscosity computed using N2 collision diameter, k, rho, T, and p.']);
end