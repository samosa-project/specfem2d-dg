% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         Based on equations (29), (30), (31), and (A6) of [Carcione et al., 1988].

% [Carcione et al., 1988] Carcione, J. M., Kosloff, D., & Kosloff, R. (1988). Wave propagation simulation in a linear viscoelastic medium. Geophysical Journal International, 95(3), 597-611.

clear all;
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
f=logspace(-2,2,1000);

% Relaxation times.
% tenu1 = [ 0.68618123823217914       0.11484954959686124       2.52182139147004644E-002  4.61038115634616451E-003]; tsnu1 = [ 0.59955862363404411       0.10472199817747203       2.29507982576383947E-002  3.98011786537642853E-003]; tenu2 = [ 0.76541298418943360       0.12452976261503652       2.75274529794950519E-002  5.33944660760637090E-003]; tsnu2 = [ 0.58414505070274858       0.10242051624315561       2.24560942265705048E-002  3.87816212327913917E-003];
tenu1=[0.25313396195527554,1.09163768103025405E-002]; tsnu1=[0.22962286350984851,9.85685340323495633E-003]; tenu2=[0.26508359588818803,1.16512236085833335E-002]; tsnu2=[0.21678194278647300,9.34322053214898390E-003];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of quality      %
% factors.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=length(tenu1); % Number of Zener solids.
omega=f*2*pi; % Pulsation.

% Method from straight reading of [Carcione et al., 1988].
% TODO: check this method.
M1C = 1 - L;
M2C = 1 - L;
for ii=1:L
  M1C = M1C + (1+1i*omega*tenu1(ii))./(1+1i*omega*tsnu1(ii));
  M2C = M2C + (1+1i*omega*tenu2(ii))./(1+1i*omega*tsnu2(ii));
end
M1C = M1C/L;
M2C = M2C/L;

% Alternative method, yielding different results (by Quentin Brissaud).
% TODO: check this method.
M1C = 0;
M2C = 0;
for ii=1:L
  M1C = M1C + (1+1i*omega*tenu1(ii))./(1+1i*omega*tsnu1(ii));
  M2C = M2C + (1+1i*omega*tenu2(ii))./(1+1i*omega*tsnu2(ii));
end
M1C = M1C/L;
M2C = M2C/L;

invQ1=imag(M1C+(d-1)*M2C)./real(M1C+(d-1)*M2C); % Inverse quality factor for P-waves. [Carcione et al., 1988], equation (29).
invQ2=imag(M2C)./real(M2C); % Inverse quality factor for S-waves. [Carcione et al., 1988], equation (30).
invQb=imag(M1C)./real(M1C); % Inverse quality factor for bulk. [Carcione et al., 1988], equation (31).

figure();
semilogx(f,invQ1,'displayname','P-waves');
semilogx(f,invQ2,'displayname','S-waves');
semilogx(f,invQb,'displayname','bulk');
grid on;
xlabel("$f$ (Hz)");
ylabel("$Q^{-1}$");
title("Inverse $Q$ factor");