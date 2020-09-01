% Author:        Léo Martire.
% Description:   Converts relaxation parameters (f_r, S_{vib}) to relaxation times (tau_\sigma, tau_\epsilon).
%                Implements Garcia et al. (2017)'s Equation (11).
% Notes:         Garcia, R. F., Brissaud, Q., Rolland, L. M., Martin, R., Komatitsch, D., Spiga, A., et al. (2017). Finite-Difference Modeling of Acoustic and Gravity Wave Propagation in Mars Atmosphere: Application to Infrasounds Emitted by Meteor Impacts. Space Science Reviews, 211(1–4), 547–570. https://doi.org/10.1007/s11214-016-0324-6.
%
% Usage:
%   [TAUSIG, TAUEPS] = frsvib2tausigtaueps(FR, SVIB)
% with:
%   FR     [Hz] the relaxation frequency,
%   SVIB   [1] the relaxation strength,
% yields:
%   TAUSIG [s] the stress relaxation time,
%   TAUEPS [s] the deformation relaxation time.

function [TAUSIG, TAUEPS] = frsvib2tausigtaueps(FR, SVIB)
  one_over_twopifr = 1 ./ (2 * pi * FR);
  TAUSIG = 0.5 * one_over_twopifr .* (-SVIB + sqrt(SVIB.^2. + 4.));
  TAUEPS = TAUSIG + SVIB .* one_over_twopifr;
end