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

function [kappa] = thermalConductivity(T)
  T = reshape(T, numel(T), 1);
  
  kappa = 2.64638e-3 * T.^(1.5) ./ (T + 245.4 * 10.^(-12 ./ T));
  
  disp(['[',mfilename,'] Thermal conductivity computed using T and USSA76''s emprical formula.']);
end