% Author:        Léo Martire.
% Description:   A useful struct, with thermodynamical parameters.
% Notes:         N/A.
%
% Usage:
%   theStruct = thermodynamicalParameters()
% with:
%   N/A,
% yields:
%   theStruct the useful structure.
%
% .mu.idealN2   computes using using ideal gas hypothesis, and a full N2
%                 composition
%               must be called with (rho, T, p), density, temperature,
%                 and pressure
%
% .mu.Suther    computes using Sutherland's empirical formula
%               must be called with (T, BETA, S), temperature, and
%                 the two Sutherland parameters (BETA, S):
%                   for Earth, choose (1.458e-6, 110.4) [NASA, 1976, Eq. (51)
%                   for Mars,  choose (1.490e-6, 217.0) [NASA, 1976, Eq. (51)
% 
% .kappa.Eucken computes using the Eucken expression
%               must be called with (MU, CV), dynamic viscosity, and
%                 isochoric specific heat
%
% .kappa.ussa76 computes using USSA76's empirical formula
%               must be called with T, temperature

function theStruct = thermodynamicalParameters()
  theStruct.mu.idealN2   = @mu__Ideal_N2_Gas; % air dynamic viscosity [kg s^{-1} m^{-1}] using the ideal gas hypothesis
  theStruct.mu.Suther    = @mu__Sutherland; % air dynamic viscosity [kg s^{-1} m^{-1}] using the ideal gas hypothesis
  theStruct.kappa.ussa76 = @kappa__USSA76; % air thermal conductivity [kg m s^{-3} K^{-1}] using USSA76 empirical formula
  theStruct.kappa.Eucken = @kappa__Eucken; % Mars, air thermal conductivity [kg m s^{-3} K^{-1}] using Eucken's expression
  
  disp(['[',mfilename,'] Structure loaded. Run ''help ',mfilename,''' to get help.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MU                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mu] = mu__Ideal_N2_Gas(rho, T, p)
  rho = reshape(rho, numel(rho), 1);
  T = reshape(T, numel(T), 1);
  p = reshape(p, numel(p), 1);
  TC = thermodynamicalConstants();
  k = TC.k;
  dN2 = TC.N2.d;
  mu = (2 * k / (3 * pi^(1.5) * dN2^2)) .* T .* (rho ./ p).^(0.5);
  % TODO: upgrade this formula by taking into account polyatomic mixtures.
  disp(['[',mfilename,'] Dynamic viscosity ',dynamicViscosityUnit,' computed using ideal N2 gas hypothesis.']);
  disp([blanks(numel(mfilename)),'   Using N2 collision diameter, k, rho, T, and p.']);
  disp([blanks(numel(mfilename)),'   See [Chapman, 1991, Eq. (5.21, 3-4)].'])
end

function [mu] = mu__Sutherland(T, BETA, S)
  T = reshape(T, numel(T), 1);
  if(not(numel(BETA)==numel(S) && numel(S)==1))
    error(['[',mfilename,', ERROR] BETA and S must have only one value.']);
  end
  mu = BETA * T.^(0.5) ./ (1 + S./T);
  disp(['[',mfilename,'] Dynamic viscosity computed using Sutherland''s formula.']);
  disp([blanks(numel(mfilename)),'   Used only T.']);
  disp([blanks(numel(mfilename)),'   See [Sutherland, 1893].']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KAPPA                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kappa = kappa__USSA76(T)
  kappa = 2.64638e-3 * T.^(1.5) ./ (T + 245.4 * 10.^(-12 ./ T));
  disp(['[',mfilename,'] Thermal conductivity computed using T and USSA76''s emprical formula.']);
  disp([blanks(numel(mfilename)),'   Used only T.']);
  disp([blanks(numel(mfilename)),'   See [NASA, 1976, Eq. (53)].']);
end

function kappa = kappa__Eucken(MU, CV)
  R = thermodynamicalConstants(); R=R.R;
  kappa = 0.25*15*R*MU * (4*CV/(15*R) + 3/5);
  disp(['[',mfilename,'] Thermal conductivity computed using the Eucken expression, for a pure gas, composed of polyatomic molecules.']);
  disp([blanks(numel(mfilename)),'   See [Bass and Chambers, 2001], and [Hirschfelder et al., 2009]’s (8.2-33), page 534.']);
end