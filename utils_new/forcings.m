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
T0 = 14; % Temporal period.
nT0 = 10.5; % Number of temporal periods.
L0 = 200; % Spatial period.
nL0 = 160; % Number of spatial periods (total).
A1 = 1;
A2 = 1;
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
s_displ = A1*cos(k_0 * X + w_0 * T + pi/2) + A2*cos(k_0 * X - w_0 * T - pi/2); % Interfering monofrequency waves, displacement (pi/2 phase shifts in order not to generate discontinuities).
% s = A1*sin(k_0 * X - w_0 * T) + A2*sin(-k_0 * X - w_0 * T); % Interfering monofrequency waves, velocity.
%%%%%%%%%%%%%%%%%%%%%%%
% Send to frequency   %
% range.              %
%%%%%%%%%%%%%%%%%%%%%%%
S_displ = fftshift(fft2(s_displ));
S_displ_R = real(S_displ);
S_displ_I = imag(S_displ);
%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.          %
%%%%%%%%%%%%%%%%%%%%%%%
% Convolve reals with wide Gaussian, replace imaginaries with random phase.
% sig=1; % Spread of Gaussian mask (should be <2);
% GMask=fspecial('gaussian', (floor(min((w_0/(2*pi))/df,(k_0/(2*pi))/dk))-1)*[1,1],0.5*sig); % Create Gaussian mask which covers f=]0,2/T0[ x k=]0,2/L0[.
% GMask=fspecial('gaussian', [floor((w_0/(2*pi))/df),floor((k_0/(2*pi))/dk)]-1,0.5*sig); % Create Gaussian mask which covers f=]0,2/T0[ x k=]0,2/L0[.
% GMask=fspecial('gaussian', [floor((k_0/(2*pi))/dk),floor((w_0/(2*pi))/df)]-1,0.5*sig); % Create Gaussian mask which covers f=]0,2/T0[ x k=]0,2/L0[.
sig_T = ((1 / T0) / 3) / 1.5; % Width of the Gaussian mask in time.
sig_X = ((1 / L0) / 3) / 1.5; % Width of the Gaussian mask in space.
[Fgm, Kgm] = meshgrid(0:df:2 / T0, 0:dk:2 / L0);
GMask = exp(- ((Fgm - 1 / T0) .^ 2 / (2 * sig_T ^ 2) + (Kgm - 1 / L0) .^ 2 / (2 * sig_X ^ 2)));
GMask = GMask / max(max(GMask)); % Normalise, not to mess with maximum amplitude.
% disp(0.5./(size(GMask).*[dk,df])) % Display span of mask in physical quantities (how much around base periods is now considered).
% Option 1: full random phase, relying on real part of spectrum.
RPhase = exp(1j * 2 * pi * rand(size(S_displ)));
% RPhase = 1; % No phase.
S_displ_R_convolved = conv2(S_displ_R, GMask, 'same');
S_displ_new = S_displ_R_convolved .* RPhase;
% % Option 2: separate real (convolve) and imaginary (add noise) parts.
% SnR=conv2(SR,GMask,'same');
% SnI=SI+(2*rand(size(S))-1)*0.01*(max(max(SI))-min(min(SI)))*0.1;
% Sn=SnR+1j*SnI;
% Mirror and convolve imaginary part in frequency.
S_displ_new_R = real(S_displ_new); % Save real part.
S_displ_new_I = imag(S_displ_new); % Extract imaginary part.
S_displ_new_I = S_displ_new_I + conj(fliplr(S_displ_new_I)); % Mirroring.
S_displ_new = S_displ_new_R + 1j * S_displ_new_I; % To be sent back to spacetime range.
S_displ_new_R = real(S_displ_new); % For plotting purposes only.
S_displ_new_I = imag(S_displ_new); % For plotting purposes only.
%%%%%%%%%%%%%%%%%%%%%%%
% Send back to        %
% spacetime range.    %
%%%%%%%%%%%%%%%%%%%%%%%
s_displ_new = real(ifft2(fftshift(S_displ_new)));
% Low-pass in time.
fcut = min(20 / T0,0.99*f_max); % Must be < f_max=0.5/dt;
for ix = 1:length(x)
  [low_pass, ~] = custom_filter(t, s_displ_new(ix, :), fcut);
  s_displ_new(ix, :) = low_pass;
end
clear('ix');
% Low-pass in space (mainly safeguard).
fcut = min(20 / L0,0.99*k_max); % Must be < k_max=0.5/dx;
for it = 1:length(t)
  [low_pass, ~] = custom_filter(x, s_displ_new(:, it), fcut);
  s_displ_new(:, it) = low_pass;
end
clear('it');
clear('fcut');
%%%%%%%%%%%%%%%%%%%%%%%
% Plot process.       %
%%%%%%%%%%%%%%%%%%%%%%%
if(0)
  figure();
  subplot(331); surf(T,X,s_displ,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(t), max(t), min(x), max(x)]);
  subplot(332); surf(F,K,S_displ_R,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(f), max(f), min(k), max(k)]);
  subplot(333); surf(F,K,S_displ_I,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(f), max(f), min(k), max(k)]);
  subplot(335); surf(F,K,S_displ_R_convolved,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(f), max(f), min(k), max(k)]);
  subplot(337); surf(T,X,s_displ_new,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(t), max(t), min(x), max(x)]);
  subplot(338); surf(F,K,S_displ_new_R,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(f), max(f), min(k), max(k)]);
  subplot(339); surf(F,K,S_displ_new_I,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(f), max(f), min(k), max(k)]);
end
%%%%%%%%%%%%%%%%%%%%%%%
% Convert signal in   %
% displacement to     %
% signal in velocity. %
%%%%%%%%%%%%%%%%%%%%%%%
[s_vel_new,~]=gradient(s_displ_new,t,x);
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
apox([1,end])=0; % Tweak to make sure forcing starts at 0.
FORCING = s_vel_new .* apox'.*apot;
% Ask user to verify forcing is ok.
forcok=-1;
disp(['  Forcing spans [',num2str(MINX), ', ', num2str(MAXX), '] m and [0, ',num2str(MAXTIME),'] s.']);
while(not(ismember(forcok,[0,1])))
  forcok=input('  Is that ok (0 for no, 1 for yes)? > ');
end
if(forcok==0)
  error('  Forcing was not ok, re-chose parametrisation in script.');
end
clear('forcok');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Interpolate on the       %
% SPECFEM mesh.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here, it is sufficient to configure the next 2 sections (SPECFEM time 
% and SPECFEM mesh).

%%%%%%%%%%%%%%%%%%%%%%%
% SPECFEM  time.      %
%%%%%%%%%%%%%%%%%%%%%%%
dtSPCFM = 1.5d-2;
nstageSPCFM = 1;
t0SPCFM = - 2 * dtSPCFM;
tSPCFM = 0:dtSPCFM / nstageSPCFM:MAXTIME*1.1;
%%%%%%%%%%%%%%%%%%%%%%%
% Reproduce SPECFEM x %
% mesh.               %
%%%%%%%%%%%%%%%%%%%%%%%
% This array must reproduce exactly mesh at z=0 (with GLL points).
nx = 1060;
xSPCFM = linspace(- 45e3, 45e3, nx + 1);
GLL = [0, 0.345346329292028 / 2, 0.5, 1.654653670707980 / 2]; % UGLY METHOD, I'M SORRY
xSPCFMwGLL = [];
for i = 1:length(xSPCFM) - 1
  xSPCFMwGLL = [xSPCFMwGLL, xSPCFM(i) + (xSPCFM(i + 1) - xSPCFM(i)) * GLL];
end
xSPCFMwGLL = [xSPCFMwGLL, xSPCFM(end)];
% Ask user to verify mesh is ok.
meshok=-1;
disp(['  Interpolating mesh (final mesh) spans [',num2str(min(xSPCFMwGLL)), ', ', num2str(max(xSPCFMwGLL)), '] m with ',num2str(nx), ' points. This mesh HAS TO match SPECFEM''s mesh.']);
while(not(ismember(meshok,[0,1])))
  meshok=input('  Is that ok (0 for no, 1 for yes)? > ');
end
if(meshok==0)
  error('  Mesh was not ok, re-chose parametrisation in script.');
end
clear('meshok');
% If an external mesh is used, the positions of points at z=0 must be entered here instead.
%%%%%%%%%%%%%%%%%%%%%%%
% Prepare meshgrid    %
% and interpolate.    %
%%%%%%%%%%%%%%%%%%%%%%%
[TSPCFM, XSPCFM] = meshgrid(tSPCFM, xSPCFMwGLL);
FORCING_INTERP = interp2(T, X, FORCING, TSPCFM, XSPCFM); % Linear interpolation.
FORCING_INTERP(isnan(FORCING_INTERP)) = 0; % Set zeros instead of NaNs outside of microbarom zone.
%%%%%%%%%%%%%%%%%%%%%%%
% Plot result (for    %
% a visual inspection %
% of the quality of   %
% interpolation).     %
%%%%%%%%%%%%%%%%%%%%%%%
if (not(any(size(FORCING) > 5000) || any(size(FORCING_INTERP) > 5000)))
  figure();
  subplot(121); surf(T, X, FORCING, 'edgecolor', 'interp', 'facecolor', 'interp'); view([0, 0, 1]); axis([min(t), max(t), min(x), max(x)]);
  subplot(122); surf(TSPCFM, XSPCFM, FORCING_INTERP, 'edgecolor', 'interp', 'facecolor', 'interp'); view([0, 0, 1]); axis([min(t), max(t), min(x), max(x)]);
end
if(0==1)
  pcolor(TSPCFM, XSPCFM, FORCING_INTERP); shading interp; axis([min(t), max(t), min(x), max(x)]);
end
%%%%%%%%%%%%%%%%%%%%%%%
% Output automatic    %
% warnings.           %
%%%%%%%%%%%%%%%%%%%%%%%
dxSPCFM = mean(diff(xSPCFMwGLL));
disp('  Interpolation resolution:');
disp(['    dxSPCFM/dx = ', num2str(dxSPCFM / dx)]);
disp(['    dtSPCFM/dt = ', num2str(dtSPCFM / dt)]);
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
%%%%%%%%%%%%%%%%%%%%%%%
% Test data.          %
%%%%%%%%%%%%%%%%%%%%%%%
% EXPORTFILEDIR = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_external_forcing/';
EXPORTFILEDIR = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/stratobaro_test_EBF/';
% tSPCFM=0:4e-4:1;
% xSPCFM=linspace(-50,50,51);
% GLL=[0, 0.345346329292028/2, 0.5, 1.654653670707976/2]; % UGLY METHOD, I'M SORRY
% xSPCFMwGLL=[];
% for i=1:length(xSPCFM)-1
%   xSPCFMwGLL=[xSPCFMwGLL,xSPCFM(i)+(xSPCFM(i+1)-xSPCFM(i))*GLL];
% end
% xSPCFMwGLL = [xSPCFMwGLL, xSPCFM(end)];
% [TSPCFM,XSPCFM]=meshgrid(tSPCFM,xSPCFMwGLL);
% MAXTIME=0.5;
% MINX=-25;
% MAXX=25;
% LAPO=7.5; apox = 0.25 .* (1. - erf((XSPCFM - MAXX + 0.5 * LAPO) / (0.25 * LAPO))) .* (1 + erf((XSPCFM - MINX - 0.5 * LAPO) / (0.25 * LAPO)));
% TAPO=0.1; apot0 = 0.5 .* (1 + erf((TSPCFM - 0.5 * TAPO) / (0.25 * TAPO))); apot = 0.5 .* (1. - erf((TSPCFM - MAXTIME + 0.5 * TAPO) / (0.25 * TAPO))); apot0(:,1)=0;
% % FORCING_INTERP=(abs(XSPCFM)<MAXX).*(TSPCFM<MAXTIME).*sin(XSPCFM/8).*sin(TSPCFM*12.5); % Test forcing function.
% % FORCING_INTERP=(abs(XSPCFM)<MAXX).*(TSPCFM<MAXTIME).*sin(XSPCFM/(0.25*8)).*sin(TSPCFM*4*12.5); % Test forcing function.
% % FORCING_INTERP=(abs(XSPCFM)<MAXX).*(TSPCFM<MAXTIME).*apox.*apot0.*apot; % Test forcing function.
% FORCING_INTERP=(abs(XSPCFM)<MAXX).*(TSPCFM<MAXTIME).*sin(TSPCFM*2*pi/0.25).*apox.*apot0.*apot; % Test forcing function.
%%%%%%%%%%%%%%%%%%%%%%%
% Detect relevant     %
% indices, and print  %
% corresponding       %
% values to file.     %
%%%%%%%%%%%%%%%%%%%%%%%
if(0)
  EXPORTFILENAME = [EXPORTFILEDIR, 'external_bottom_forcing.dat'];
  % Ask user to verify export path is ok.
  pathok=-1;
  disp(['  Export path is ''',EXPORTFILENAME,'''.']);
  while(not(ismember(pathok,[0,1])))
    pathok=input('  Is that ok (0 for no, 1 for yes)? > ');
  end
  if(pathok==0)
    error('  Export path was not ok, re-chose in script.');
  end
  clear('pathok');
  % Export.
  bytespervalue=12.956059264925035;
  expectedsize=prod(size(FORCING_INTERP))*bytespervalue;
  disp([' File will be ', num2str(expectedsize), ' bytes (',num2str(expectedsize/1024),' kB, ',num2str(expectedsize/1048576),' MB).']);
  itstop = find(abs(TSPCFM(1, :) - MAXTIME) == min(abs(TSPCFM(1, :) - MAXTIME))) + 1;
  ixmin = max(find(abs(XSPCFM(:, 1) - MINX) == min(abs(XSPCFM(:, 1) - MINX))) - 1, 1);
  ixmax = min(find(abs(XSPCFM(:, 1) - MAXX) == min(abs(XSPCFM(:, 1) - MAXX))) + 1, length(XSPCFM(:, 1)));
  f_new = fopen(EXPORTFILENAME, 'w');
  for it = 1:itstop
    for ix = ixmin:ixmax
      fprintf(f_new, '%.5e %.8e %.5e', TSPCFM(1, it), XSPCFM(ix, 1), FORCING_INTERP(ix, it));
      fprintf(f_new, "\n");
    end
    if (ismember(it, floor((1:10) * 0.1 * itstop)))
       disp(['Writing to file (', num2str(ceil(it / itstop * 100)), ' % complete).']);
    end
  end
  fclose('all');
end