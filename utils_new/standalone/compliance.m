% Author:        Léo Martire.
% Description:   TODO.
% Notes:         [Sorrells, 1971, Kenda et al., 2017].
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

% Sorrells, G. G. (1971). A Preliminary Investigation into the Relationship between Long‐Period Seismic Noise and Local Fluctuations in the Atmospheric Pressure Field. Geophysical Journal of the Royal Astronomical Society, 26(1–4), 71–82. https://doi.org/10.1111/j.1365-246X.1971.tb03383.x
% Kenda, B., Lognonné, P., Spiga, A., Kawamura, T., Kedar, S., Banerdt, W. B., et al. (2017). Modeling of Ground Deformation and Shallow Surface Waves Generated by Martian Dust Devils and Perspectives for Near-Surface Structure Inversion. Space Science Reviews, 211(1–4), 501–524. https://doi.org/10.1007/s11214-017-0378-0
% Mavko, G., Mukerji, T., & Dvorkin, J. (2009). The Rock Physics Handbook (2nd ed.). Cambridge University Press.

function [CompZ, CompH] = compliance(freq, SPEED, rho, vp, vs)
  % TODO: no denpendency on freq?
  
  % 1st lamé parameter and shear modulus [Mavko et al, 2009]
  lambda = rho * (vp^2 - 2*vs^2);
  mu = rho * vs^2;

  % elastic coefficients in expressions [Sorrells, 1971, Kenda et al., 2017]
  Coef1 = (lambda+2*mu)/(2*mu*(lambda+mu));
  Coef2 = 1/(2*(lambda+mu));
  
  CompZ = 0.0*freq - SPEED*Coef1*complex(0.0, 1.0);
  
  CompH = 0.0*freq - SPEED*Coef2;
  
  disp(['[',mfilename,'] Computed compliance coefficients from [Kenda et al., 2017]''s formulas, using (\rho, v_p, v_s) to get (\lambda, \mu).']);
end
