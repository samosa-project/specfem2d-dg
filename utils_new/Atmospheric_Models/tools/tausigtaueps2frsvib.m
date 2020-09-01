% Author:        Léo Martire.
% Description:   Converts relaxation times (tau_\sigma, tau_\epsilon) to relaxation parameters (f_r, S_{vib}).
%                Implements Garcia et al. (2017)'s Equation (10).
% Notes:         Garcia, R. F., Brissaud, Q., Rolland, L. M., Martin, R., Komatitsch, D., Spiga, A., et al. (2017). Finite-Difference Modeling of Acoustic and Gravity Wave Propagation in Mars Atmosphere: Application to Infrasounds Emitted by Meteor Impacts. Space Science Reviews, 211(1–4), 547–570. https://doi.org/10.1007/s11214-016-0324-6.
%
% Usage:
%   [FR, SVIB] = tausigtaueps2frsvib(TAUSIG, TAUEPS)
% with:
%   TAUSIG [s] the stress relaxation time,
%   TAUEPS [s] the deformation relaxation time,
% yields:
%   FR     [Hz] the relaxation frequency,
%   SVIB   [1] the relaxation strength.

function [FR, SVIB] = tausigtaueps2frsvib(TAUSIG, TAUEPS)
  TAUR = sqrt(TAUSIG .* TAUEPS);
  FR = 1 ./ (2 * pi * TAUR);
  SVIB = (TAUEPS - TAUSIG) ./ TAUR;
end