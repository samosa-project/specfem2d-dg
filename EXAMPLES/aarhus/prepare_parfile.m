clear all;
close all;
clc;

R=8.31446261815324;
k=1.380649e-23;

p = 2000;
T = 273.15+21.6;
sound_velocity = 271.594; sound_velocity_desc = 'experimental setup';

% Molar fractions (CO2 attenuation).
Xco2 = 0.99;
Xn2 = 0;
Xar = 0;
Xh2o = 0.01;


% (cP, cV) in m^2.s^{-2}.K^{-1} for pure CO2 at 20°C from https://www.engineeringtoolbox.com/specific-heat-capacity-gases-d_159.html
cP = 844;
cV = 655;
gamma = cP/cV;
cpcv_desc = ['from pure CO2 at 20°C (yields gamma=',sprintf('%.5f',gamma),') from https://www.engineeringtoolbox.com/specific-heat-capacity-gases-d_159.html'];

dco2 = 3.30e-10; % kinetic diameter of CO2 from [https://en.wikipedia.org/wiki/Kinetic_diameter]=Ismail, Ahmad Fauzi; Khulbe, Kailash; Matsuura, Takeshi, Gas Separation Membranes: Polymeric and Inorganic, Springer, 2015 ISBN 3319010956.

rho = p*(cP/cV)/(sound_velocity^2); rho_desc = ['p*(cP/cV)/(sound_velocity^2);'];
mu = ( k * T * sqrt(rho/p)) / (pi^(3/2) * dco2^2); mu_desc = ['(k*T*sqrt(rho/p))/(pi^(3/2)*dco2^2) with dco2=',sprintf('%.2e', dco2),''];
kappa = 0.25*(15*R*mu)*(((4*cV)/(15*R)) + 3/5); kappa_desc = ['from Eucken expression 0.25*(15*R*mu)*(((4*cv)/(15*R)) + 3/5)'];


% CO2 ATTENUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mostly from Raphaël's MCD script.
% All C's must be in L2MT-2TH-1N-1 (cf. thèse).
C_V_inf = (5./2.)*R;
C_P_inf = C_V_inf*gamma;
C_V_0 = (5./2.)*R;
C_P_0 = C_V_0*gamma;
theta = 960.0; % Degenerate bending mode temperature (Bass, 2001).
Cprime = R * Xco2 * (theta/T)^2 * exp(-theta/T) / ((1-exp(-theta/T))^2); % (Bass, 2001, Eq. (7)).
% Rate of energy transfer from CO_2 during collisions with CO_2.
% Don't now where these formulas come from.
kco2 = ((0.219*p)/(mu))*(exp(-60.75/((T^(1./3.)))));
kn2 = ((1.44*p)/(mu))*exp(-78.29/((T^(1./3.))));
kar = kn2;
kh2o = ((6.e-2)*p)/mu;
kk = (Xco2*kco2)+(Xn2*kn2)+(Xar*kar)+(Xh2o*kh2o);
tau_VT = 1.0 / (kk*(1.0-exp(-theta/T)));
tau_VS = (C_P_inf/C_P_0)*tau_VT;
FR = 1./(2.*pi*tau_VS);
SVIB = (Cprime*R)/(C_P_inf*(C_V_inf+Cprime));
addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/Atmospheric_Models/tools'));
[tau_SIG, tau_EPS] = frsvib2tausigtaueps(FR, SVIB);

freq = 2090;
a_vib = (((pi*SVIB)./sound_velocity) .* ((freq.^2)./FR)) ./ (1 + (freq./FR).^2);
tauepssig_desc = ['CO2 relaxation modelling an alpha_vib(f=',sprintf('%.0f',freq),')=',sprintf('%.3e',a_vib),', see Bass (2001, 10.1121/1.1365424) and Garcia (2017, 10.1007/s11214-016-0324-6).'];

Zrot = 61.1*exp(-16.8*T^(-1/3));
arot_factor = 1+3*gamma*(gamma-1)*R*Zrot/(4*1.25*C_P_0);
mu_withrot = mu * arot_factor;
mu_desc = ['classical mu=',sprintf('%.8e',mu),' (defined as ',mu_desc,'), but modified by a factor (1+3*gamma*(gamma-1)*R*Zrot/(4*1.25*C_P_0)) to account for rotational attenuation (see these)'];

% PRINT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['# Fluid model (if MODEL\=''external''). Generated from ',mfilename('fullpath'),'.']);
disp(['USE_ISOTHERMAL_MODEL            = .false.        # Isobaric model for the fluid part.']);
disp(['SCALE_HEIGHT                    = 0.             # Only used if USE_ISOTHERMAL_MODEL==.true..']);
disp(['surface_density                 = ',sprintf('%.12f',rho),' # $\rho_0$ [kg.m^{-3}]: ', rho_desc, '.']);
disp(['sound_velocity                  = ',sprintf('%.10f',sound_velocity),' # $c$ [m.s^{-1}]: ', sound_velocity_desc, '.']);
disp(['wind                            = 0.             # No wind.']);
disp(['gravity                         = 0.             # Only used if USE_ISOTHERMAL_MODEL==.true., if USE_ISOTHERMAL_MODEL==.false. the gravity field is anyhow forced to 0.']);
disp(['dynamic_viscosity               = ',sprintf('%.8e',mu_withrot),' # $\mu$, dynamic viscosity [kg.s^{-1}.m^{-1}]: ',mu_desc,'.']);
disp(['thermal_conductivity            = ',sprintf('%.8e',kappa),' # $\kappa$, thermal conductivity [kg.m.s^{-3}.K^{-1}]: ',kappa_desc,'.']);
disp(['tau_epsilon                     = ',sprintf('%.8e',tau_EPS),' # $\tau_\epsilon$, strain relaxation time [s]: ',tauepssig_desc,'.']);
disp(['tau_sigma                       = ',sprintf('%.8e',tau_SIG),' # $\tau_\sigma$,   stress relaxation time [s]: ',tauepssig_desc,'.']);
disp(['constant_p                      = ',sprintf('%.10f',cP),' # $c_p$, isobaric  specific heat capacity [m^2.s^{-2}.K^{-1}]: ',cpcv_desc,'.']);
disp(['constant_v                      = ',sprintf('%.10f',cV),' # $c_v$, isochoric specific heat capacity [m^2.s^{-2}.K^{-1}]: ',cpcv_desc,'.']);
disp(' ');
disp(' ');