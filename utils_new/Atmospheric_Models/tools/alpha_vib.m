% Author:        Léo Martire.
% Description:   Computes the vibrationnal attenuation coefficient.
%                Implements Garcia et al. (2017)'s Equation (10).
% Notes:         Garcia, R. F., Brissaud, Q., Rolland, L. M., Martin, R., Komatitsch, D., Spiga, A., et al. (2017). Finite-Difference Modeling of Acoustic and Gravity Wave Propagation in Mars Atmosphere: Application to Infrasounds Emitted by Meteor Impacts. Space Science Reviews, 211(1–4), 547–570. https://doi.org/10.1007/s11214-016-0324-6.
%
% Usage:
%   [alpha_vib] = alpha_vib(freq, C, FR, SVIB)
% with:
%   freq      [Hz] the queried frequencies,
%   C         [m/s] the speed of sound,
%   FR        [Hz] the vibrational relaxation frequency,
%   SVIB      [1] the vibrational relaxation strength,
% yields:
%   alpha_vib [m^{-1}] the vibrationnal attenuation coefficient

function [alpha_vib] = alpha_vib(freq, C, FR, SVIB)
  % Make column.
  freq = reshape(freq, 1, max(size(freq)));
  C = reshape(C, max(size(C)), 1);
  FR = reshape(FR, max(size(FR)), 1);
  SVIB = reshape(SVIB, max(size(SVIB)), 1);
  
  % Check provided dimensions.
  dimensions = [size(C); size(FR); size(SVIB)];
  if(any(dimensions(:,1)/dimensions(1,1)~=1) || any(dimensions(:,2)/dimensions(1,2)~=1))
    disp(['  [',mfilename,', ERROR] Dimensions of inputs:']);
    dimensions
    error(['  [',mfilename,', ERROR] Dimensions of all atmospheric datasets must agree.']);
  end
  
  % Compute.
  alpha_vib = (((pi*SVIB)./C) .* ((freq.^2)./FR)) ./ (1 + (freq./FR).^2);
end

