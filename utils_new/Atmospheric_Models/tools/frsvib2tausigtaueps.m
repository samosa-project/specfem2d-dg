% Author:        Léo Martire.
% Description:   Converts relaxation parameters (f_r, S_{vib}) to relaxation times (tau_\sigma, tau_\epsilon). Implements 10.1007/s11214-016-0324-6, equation (11).
% Notes:         See 10.1007/s11214-016-0324-6.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function [TAUSIG, TAUEPS] = frsvib2tausigtaueps(FR, SVIB)
  one_over_twopifr = 1./(2*pi*FR);
  TAUSIG = 0.5*one_over_twopifr.*(-SVIB + sqrt(SVIB.^2. + 4.));
  TAUEPS = TAUSIG + SVIB.*one_over_twopifr;
end