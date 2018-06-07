% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         Configure the 3 steps according to your needs.
% Notes:         Hardcoded bottom forcings (DG extension) tests were moved to 'forcings_tests.m'.

clear all;
% close all;
clc;
format compact;
set(0, 'DefaultLineLineWidth', 3); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 8); % Default at 6.
set(0, 'defaultTextFontSize', 22);
set(0, 'defaultAxesFontSize', 22); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

set(groot, 'defaultSurfaceEdgeColor', 'none');

% 3 steps:
% 1) Prepare forcing in spacetime.
% 2) Interpolate forcing on SPECFEM timesteps and space mesh.
% 3) Export to well formatted file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Prepare forcing in       %
% spacetime.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here, prepare (however you want) a meshgrid [T, X], a forcing matrix
% FORCING which rely on [X,T], a minimum abscissa MINX (under which forcing
% is 0), a maximum abscissa MAXX (above which forcing is 0), and a maximum
% time MAXTIME (above which forcing is 0).

% Microbaroms. %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.         %
%%%%%%%%%%%%%%%%%%%%%%%
N = 21; % Points per period.
T0 = 7; % Temporal period.
nT0 = 10.5; % Number of temporal periods to represent.
L0 = 200; % Spatial period.
nL0 = 30; % Number of spatial periods to represent.
%%%%%%%%%%%%%%%%%%%%%%%
% Ranges.             %
%%%%%%%%%%%%%%%%%%%%%%%
dt = T0 / N; % Temporal resolution.
dx = L0 / N; % Spatial resolution.
% Temporal frequency range.
f_max = 0.5 / dt;
df = 1 / (nT0 * T0);
f = - f_max:df:f_max;
% Spatial frequency range.
k_max = 0.5 / dx;
dk = 1 / (nL0 * L0);
k = - k_max:dk:k_max;
% Temporal and spatial spaces.
t = 0:dt:nT0 * T0;
x = - 0.5 * nL0 * L0:dx:0.5 * nL0 * L0;
% 2D ranges.
[T, X] = meshgrid(t, x);
[F, K] = meshgrid(f, k);
%%%%%%%%%%%%%%%%%%%%%%%
% First draft of      %
% signal.             %
%%%%%%%%%%%%%%%%%%%%%%%
k_0 = 2 * pi / L0;
w_0 = 2 * pi / T0;
s = sin(k_0 * X - w_0 * T) + sin(- k_0 * X - w_0 * T); % Interfering monofrequency waves.
%%%%%%%%%%%%%%%%%%%%%%%
% Send to frequency   %
% range.              %
%%%%%%%%%%%%%%%%%%%%%%%
S = fftshift(fft2(s));
SR = real(S);
SI = imag(S);
%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.          %
%%%%%%%%%%%%%%%%%%%%%%%
% Convolve reals with wide Gaussian, replace imaginaries with random phase.
% sig=1; % Spread of Gaussian mask (should be <2);
% GMask=fspecial('gaussian', (floor(min((w_0/(2*pi))/df,(k_0/(2*pi))/dk))-1)*[1,1],0.5*sig); % Create Gaussian mask which covers f=]0,2/T0[ x k=]0,2/L0[.
% GMask=fspecial('gaussian', [floor((w_0/(2*pi))/df),floor((k_0/(2*pi))/dk)]-1,0.5*sig); % Create Gaussian mask which covers f=]0,2/T0[ x k=]0,2/L0[.
% GMask=fspecial('gaussian', [floor((k_0/(2*pi))/dk),floor((w_0/(2*pi))/df)]-1,0.5*sig); % Create Gaussian mask which covers f=]0,2/T0[ x k=]0,2/L0[.
sig_T = ((1 / T0) / 3) / 5;
sig_X = ((1 / L0) / 3) / 1;
[Fgm, Kgm] = meshgrid(0:df:2 / T0, 0:dk:2 / L0);
GMask = exp(- ((Fgm - 1 / T0) .^ 2 / (2 * sig_T ^ 2) + (Kgm - 1 / L0) .^ 2 / (2 * sig_X ^ 2)));
GMask = GMask / max(max(GMask)); % Normalise, not to mess with maximum amplitude.
% disp(0.5./(size(GMask).*[dk,df])) % Display span of mask in physical quantities (how much around base periods is now considered).
% Option 1: full random phase, relying on real part of spectrum.
RPhase = exp(1j * 2 * pi * rand(size(S)));
% RPhase = 1; % No phase.
Sn = conv2(SR, GMask, 'same') .* RPhase;
% % Option 2: separate real (convolve) and imaginary (add noise) parts.
% SnR=conv2(SR,GMask,'same');
% SnI=SI+(2*rand(size(S))-1)*0.01*(max(max(SI))-min(min(SI)))*0.1;
% Sn=SnR+1j*SnI;
% Mirror and convolve imaginary part in frequency.
SnR = real(Sn); SnI = imag(Sn);
SnI = SnI + conj(fliplr(SnI));
Sn = SnR + 1j * SnI;
SnR = real(Sn);
SnI = imag(Sn);
%%%%%%%%%%%%%%%%%%%%%%%
% Send back to time   %
% range.              %
%%%%%%%%%%%%%%%%%%%%%%%
sn = real(ifft2(fftshift(Sn)));
% Filter in time.
fcut = 10 / T0; % Must be < f_max=0.5/dt;
for ix = 1:length(x)
  [signal_LP, ~] = custom_filter(t, sn(ix, :), fcut);
  sn(ix, :) = signal_LP;
end
clear('ix');
% Filter in space (mainly safeguard).
fcut = 10 / L0; % Must be < k_max=0.5/dx;
for it = 1:length(t)
  [signal_LP, ~] = custom_filter(x, sn(:, it), fcut);
  sn(:, it) = signal_LP;
end
clear('it');
clear('fcut');
%%%%%%%%%%%%%%%%%%%%%%%
% Plot process.       %
%%%%%%%%%%%%%%%%%%%%%%%
% figure();
% subplot(231); surf(T,X,s,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(t), max(t), min(x), max(x)]);
% subplot(232); surf(F,K,SR,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(f), max(f), min(k), max(k)]);
% subplot(233); surf(F,K,SI,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(f), max(f), min(k), max(k)]);
% subplot(234); surf(T,X,sn,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(t), max(t), min(x), max(x)]);
% subplot(235); surf(F,K,SnR,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(f), max(f), min(k), max(k)]);
% subplot(236); surf(F,K,SnI,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(f), max(f), min(k), max(k)]);
%%%%%%%%%%%%%%%%%%%%%%%
% Apodisation.        %
%%%%%%%%%%%%%%%%%%%%%%%
MINX = - 0.5 * nL0 * L0;
MAXX = 0.5 * nL0 * L0;
n = 5; apox = 0.25 .* (1. - erf((x - MAXX + 0.5 * n * L0) / (0.25 * n * L0))) .* (1 + erf((x - MINX - 0.5 * n * L0) / (0.25 * n * L0))); % Apodisation over n spatial period on each side.
MAXTIME = nT0 * T0;
n = 1; apot0 = 0.5 .* (1 + erf((t - 0.5 * n * T0) / (0.25 * n * T0))); % Apodisation over n temporal periods at beginning.
n = 1.5; apot = 0.5 .* (1. - erf((t - MAXTIME + 0.5 * n * T0) / (0.25 * n * T0))); % Apodisation over n temporal periods at end.
apot = apot0 .* apot;
apot(1) = 0; % Tweak to make sure forcing starts at 0.
FORCING = sn .* apox'.*apot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Interpolate on the       %
% SPECFEM mesh.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here, it is sufficient to configure the next 2 sections (SPECFEM time 
% and SPECFEM mesh).

%%%%%%%%%%%%%%%%%%%%%%%
% SPECFEM  time.      %
%%%%%%%%%%%%%%%%%%%%%%%
NSTEPSPCFM = 8000;
dtSPCFM = 1.5d-2;
istageSPCFM = 1;
t0SPCFM = - 2 * dtSPCFM;
tSPCFM = 0:dtSPCFM / istageSPCFM:NSTEPSPCFM * dtSPCFM;
%%%%%%%%%%%%%%%%%%%%%%%
% Reproduce SPECFEM x %
% mesh.               %
%%%%%%%%%%%%%%%%%%%%%%%
% This array must reproduce exactly mesh at z=0 (with GLL points).
xSPCFM = linspace(- 45e3, 45.d3, 1500 + 1);
GLL = [0, 0.345346329292028 / 2, 0.5, 1.654653670707980 / 2]; % UGLY METHOD, I'M SORRY
xSPCFMwGLL = [];
for i = 1:length(xSPCFM) - 1
  xSPCFMwGLL = [xSPCFMwGLL, xSPCFM(i) + (xSPCFM(i + 1) - xSPCFM(i)) * GLL];
end
% If an external mesh is used, the positions of points at z=0 must be entered here instead.
%%%%%%%%%%%%%%%%%%%%%%%
% Prepare meshgrid    %
% and interpolate.    %
%%%%%%%%%%%%%%%%%%%%%%%
[TSPCFM, XSPCFM] = meshgrid(tSPCFM, xSPCFMwGLL);
snapo_query = interp2(T, X, FORCING, TSPCFM, XSPCFM); % Linear interpolation.
snapo_query(isnan(snapo_query)) = 0; % Set zeros instead of NaNs outside of microbarom zone.
%%%%%%%%%%%%%%%%%%%%%%%
% Plot result (for    %
% a visual inspection %
% of the quality of   %
% interpolation).     %
%%%%%%%%%%%%%%%%%%%%%%%
if (not(any(size(FORCING) > 5000) || any(size(snapo_query) > 5000)))
  figure();
  subplot(121); surf(T, X, FORCING, 'edgecolor', 'interp', 'facecolor', 'interp'); view([0, 0, 1]); axis([min(t), max(t), min(x), max(x)]);
  subplot(122); surf(TSPCFM, XSPCFM, snapo_query, 'edgecolor', 'interp', 'facecolor', 'interp'); view([0, 0, 1]); axis([min(t), max(t), min(x), max(x)]);
end
%%%%%%%%%%%%%%%%%%%%%%%
% Output automatic    %
% warnings.           %
%%%%%%%%%%%%%%%%%%%%%%%
dxSPCFM = mean(diff(xSPCFM));
disp(['dxSPCFM/dx = ', num2str(dxSPCFM / dx)]);
disp(['dtSPCFM/dt = ', num2str(dtSPCFM / dt)]);
if (dxSPCFM / dx > 5)
  disp(['[WARNING] dxSPCFM/dx = ', num2str(dxSPCFM / dx), ' is high, subsampling can generate unwanted behaviour.']);
end
if (dtSPCFM / dt > 5)
  disp(['[WARNING] dtSPCFM/dt = ', num2str(dtSPCFM / dt), ' is high, subsampling can generate unwanted behaviour.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Export to file.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here, format shoud be compatible with the reading which is done in
% the subroutine 'prepare_external_forcing' in 'boundary_terms_DG.f90'.

%%%%%%%%%%%%%%%%%%%%%%%
% Path to file.       %
%%%%%%%%%%%%%%%%%%%%%%%
EXPORTFILEDIR = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200/';
EXPORTFILENAME = [EXPORTFILEDIR, 'external_bottom_forcing.dat'];
%%%%%%%%%%%%%%%%%%%%%%%
% Test data.          %
%%%%%%%%%%%%%%%%%%%%%%%
% EXPORTFILENAME='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_external_forcing/external_bottom_forcing.dat';
% tSPCFM=0:4e-4:1;
% xSPCFM=linspace(-50,50,51);
% GLL=[0, 0.345346329292028/2, 0.5, 1.654653670707976/2]; % UGLY METHOD, I'M SORRY
% xSPCFMwGLL=[];
% for i=1:length(xSPCFM)-1
%   xSPCFMwGLL=[xSPCFMwGLL,xSPCFM(i)+(xSPCFM(i+1)-xSPCFM(i))*GLL];
% end
% [TSPCFM,XSPCFM]=meshgrid(tSPCFM,xSPCFMwGLL);
% % snapo_query=(abs(XSPCFM)<25).*(T<0.5).*sin(XSPCFM/8).*sin(TSPCFM*12.5); % Test forcing function.
% snapo_query=(abs(XSPCFM)<25).*(TSPCFM<0.5).*sin(XSPCFM/(0.25*8)).*sin(TSPCFM*4*12.5); % Test forcing function.
% MAXTIME=0.5;
% MINX=-25;
% MAXX=30;
%%%%%%%%%%%%%%%%%%%%%%%
% Detect relevant     %
% indices, and print  %
% corresponding       %
% values to file.     %
%%%%%%%%%%%%%%%%%%%%%%%
itstop = find(abs(TSPCFM(1, :) - MAXTIME) == min(abs(TSPCFM(1, :) - MAXTIME))) + 1;
ixmin = max(find(abs(XSPCFM(:, 1) - MINX) == min(abs(XSPCFM(:, 1) - MINX))) - 1, 1);
ixmax = min(find(abs(XSPCFM(:, 1) - MAXX) == min(abs(XSPCFM(:, 1) - MAXX))) + 1, length(XSPCFM(:, 1)));
f_new = fopen(EXPORTFILENAME, 'w');
for it = 1:itstop
  for ix = ixmin:ixmax
     fprintf(f_new, '%.5e %.8e %.5e', TSPCFM(1, it), XSPCFM(ix, 1), snapo_query(ix, it));
    fprintf(f_new, "\n");
  end
  if (ismember(it, floor((1:10) * 0.1 * itstop)))
     disp(['Writing to file (', num2str(ceil(it / itstop * 100)), ' % complete).']);
  end
end
fclose('all');