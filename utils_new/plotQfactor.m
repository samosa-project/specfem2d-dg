% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   The computations are based on equations (29), (30), (31),
%                and (A6) of [Carcione et al., 1988]. Those having an
%                error (a factor 1/L), we also implement the correction
%                that was presented in equation (12) of
%                [Moczo and Kristek, 2005].
% Last modified: See file metadata.
% Usage:         Fill section "Parameters" according to the calculation to
%                be performed.
% Notes:         Relaxation times can be found by uncommenting the
%                relevant "write" instructions in the "attenuation_model"
%                subroutine in the "attenuation_model.f90" source file
%                under "/src/specfem".

% [Carcione et al., 1988] Carcione, J. M., Kosloff, D., & Kosloff, R. (1988). Wave propagation simulation in a linear viscoelastic medium. Geophysical Journal International, 95(3), 597-611.
% [Moczo and Kristek, 2005] Moczo, P., & Kristek, J. (2005). On the rheological models used for time‐domain methods of seismic wave propagation. Geophysical Research Letters, 32(1).

% clear all;
% close all;
clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 8); % Default at 6.
set(0, 'defaultTextFontSize', 22);
set(0, 'defaultAxesFontSize', 22); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 2; % Spatial dimension.

% Frequency range.
% f0=20; f=logspace(log10(f0)-1,log10(f0)+1,1000);
% f=logspace(-2,2,1000);
f=0:0.01:200;

% Relaxation times.
% tenu1 = [ 0.68618123823217914       0.11484954959686124       2.52182139147004644E-002  4.61038115634616451E-003]; tsnu1 = [ 0.59955862363404411       0.10472199817747203       2.29507982576383947E-002  3.98011786537642853E-003]; tenu2 = [ 0.76541298418943360       0.12452976261503652       2.75274529794950519E-002  5.33944660760637090E-003]; tsnu2 = [ 0.58414505070274858       0.10242051624315561       2.24560942265705048E-002  3.87816212327913917E-003];
tenu1=[0.25313396195527554,1.09163768103025405E-002]; tsnu1=[0.22962286350984851,9.85685340323495633E-003]; tenu2=[0.26508359588818803,1.16512236085833335E-002]; tsnu2=[0.21678194278647300,9.34322053214898390E-003];

% Material parameters.
rho = 1500; vp = 2000; vs = 2206;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of quality      %
% factors.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=length(tenu1); % Number of Zener solids.
omega=f*2*pi; % Pulsation.

% Corrected complex modulus from [Moczo and Kristek, 2005], equation (12).
M1C = 0;
M2C = 0;
for ii=1:L
  M1C = M1C + (1+1i*omega*tenu1(ii))./(1+1i*omega*tsnu1(ii));
  M2C = M2C + (1+1i*omega*tenu2(ii))./(1+1i*omega*tsnu2(ii));
end
M1C = M1C/L;
M2C = M2C/L;
M1 = rho * (2*vp^2 - 2*(d-1)*vs^2);
M2 = 2*rho * vs^2;
M1C = M1C*M1;
M2C = M2C*M2;

invQ1=imag(M1C+(d-1)*M2C)./real(M1C+(d-1)*M2C); % Inverse quality factor for P-waves. [Carcione et al., 1988], equation (29).
invQ2=imag(M2C)./real(M2C); % Inverse quality factor for S-waves. [Carcione et al., 1988], equation (30).
invQb=imag(M1C)./real(M1C); % Inverse quality factor for bulk. [Carcione et al., 1988], equation (31).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots.                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
loglog(f,invQ1,'displayname','P-wave attenuation'); hold on;
loglog(f,invQ2,'displayname','S-wave attenuation'); hold on;
loglog(f,invQb,'displayname','bulk attenuation'); hold on;
grid on;
xlabel("$f$ (Hz)");
ylabel("$Q^{-1}$");
title({"Inverse $Q$ factor",[num2str(L),' SLS, $\rho=', num2str(rho),'$ kg/m$^3$, $v_p=',num2str(vp),'$ m/s, $v_s=',num2str(vs),'$ m/s']});
xlim([f(1),f(end)]);
legend('location','best');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%