% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

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

tenu1 = [ 0.68618123823217914       0.11484954959686124       2.52182139147004644E-002  4.61038115634616451E-003]; tsnu1 = [ 0.59955862363404411       0.10472199817747203       2.29507982576383947E-002  3.98011786537642853E-003]; tenu2 = [ 0.76541298418943360       0.12452976261503652       2.75274529794950519E-002  5.33944660760637090E-003]; tsnu2 = [ 0.58414505070274858       0.10242051624315561       2.24560942265705048E-002  3.87816212327913917E-003];
% tenu1=[0.25313396195527554,1.09163768103025405E-002]; tsnu1=[0.22962286350984851,9.85685340323495633E-003]; tenu2=[0.26508359588818803,1.16512236085833335E-002]; tsnu2=[0.21678194278647300,9.34322053214898390E-003];

% f0=20; f=logspace(log10(f0)-1,log10(f0)+1,1000);
f=logspace(-2,2,1000);

w=f*2*pi;
mcnu1 = 1 - length(tenu1);
mcnu2 = 1 - length(tenu2);
for ii=1:length(tenu2)
  mcnu1 = mcnu1 + (1+1i*w*tenu1(ii))./(1+1i*w*tsnu1(ii));
  mcnu2 = mcnu2 + (1+1i*w*tenu2(ii))./(1+1i*w*tsnu2(ii));
end

figure();
% plot(f,imag(mcnu2)./real(mcnu2));
semilogx(f,imag(mcnu2)./real(mcnu2));
xlabel("$f$ (Hz)");
ylabel("$Q^{-1}$");
title("Inverse $Q$ factor");