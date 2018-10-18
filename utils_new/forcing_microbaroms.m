% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         Configure the 3 steps according to your needs.
% Notes:         Hardcoded bottom forcings (DG extension) tests were moved to 'forcings_tests.m'.

clear all;
close all;
clc;
format compact;
set(0, 'DefaultLineLineWidth', 3); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 8); % Default at 6.
set(0, 'defaultTextFontSize', 22);
set(0, 'defaultAxesFontSize', 22); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

set(groot, 'defaultSurfaceEdgeColor', 'none');

plot_process=0;
plot_forcing=0;

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
% Forcing related.
% T0 = 14; % Temporal period.
% T0 = 12.5; % Temporal period.
T0 = 13.54; % Temporal period from WAVEWATCH3 multi_1.partition.glo_30m.20160523060000-48500182000.
nT0 = 100; % Number of temporal periods (must be integer to prevent FFT misbehaving).
% nT0 = 2; % Number of temporal periods (must be integer to prevent FFT misbehaving).
% L0 = 200; % Spatial period.
L0 = 286.2; % Spatial period from WAVEWATCH3 multi_1.partition.glo_30m.20160523060000-48500182000.
% nL0 = 160; % Number of spatial periods (total, must be integer to prevent FFT misbehaving).
nL0 = 120; % Number of spatial periods (total, must be integer to prevent FFT misbehaving).
% nL0 = 5; % Number of spatial periods (total, must be integer to prevent FFT misbehaving).
MINX = - 0.5 * nL0 * L0; % Set forcing maximum x here.
MAXX = 0.5 * nL0 * L0; % Set forcing maximum x here.
MAXTIME = nT0 * T0; % Set forcing maximum time here.

% Frequency range related.
% invspread = 1.5; % Inverse spread factor for Gaussian convolution in frequency range (the higher the less spread).
% invspread = 2; % Inverse spread factor for Gaussian convolution in frequency range (the higher the less spread).
% invspread = 3; % Inverse spread factor for Gaussian convolution in frequency range (the higher the less spread).
% invspread = 4.5; % Inverse spread factor for Gaussian convolution in frequency range (the higher the less spread).
invspread = 5; % Inverse spread factor for Gaussian convolution in frequency range (the higher the less spread).
% invspread = 5.5; % Inverse spread factor for Gaussian convolution in frequency range (the higher the less spread).
% invspread = 6; % Inverse spread factor for Gaussian convolution in frequency range (the higher the less spread).
% invspread = 1e2; % Inverse spread factor for Gaussian convolution in frequency range (the higher the less spread).

% Amplitude related.
% A = 1; % Amplitude of waves displacement for draft signal (peak-to-peak is 2 times that).
% A_aimed=A/2; % Aimed amplitude for displacement (peak-to-peak is 2 times that).
% A_aimed=1e-2; % Aimed amplitude for velocity (peak-to-peak is 2 times that). Velocity of 1m/s results in 50 Pa amplitude (100 Pa p2p). Microbaroms are typically 0.1-0.5 Pa amplitude (0.2-1 Pa p2p), which is 1-5 microbar.
% A_aimed=1; % Aimed amplitude for velocity (peak-to-peak is 2 times that). Velocity of 1m/s results in 50 Pa amplitude (100 Pa p2p). Microbaroms are typically 0.1-0.5 Pa amplitude (0.2-1 Pa p2p), which is 1-5 microbar.
% A_aimed=7; % Aimed amplitude for displacement (peak-to-peak is 2 times that). [m].
% A_aimed=7.08; % Amplitude from WAVEWATCH3 multi_1.partition.glo_30m.20160523060000-48500182000.
A_aimed=7.08/20; % Amplitude from WAVEWATCH3 multi_1.partition.glo_30m.20160523060000-48500182000, divided by 20.

% Apodisation related.
activate_apo = 1; % Activate apodisations.
% activate_apo = 0; % Deactivate apodisations.
napoxlr=10; % Number of periods for space apodisation on left/right sides.
% napoxlr=2; % Number of periods for space apodisation on left/right sides.
napot0=1; % Number of periods for time apodisation at t=0.
napotend=napot0; % Number of periods for time apodisation at t=MAXTIME.

% SPECFEM-DG related.
% If an external mesh is used, the positions of points at z=0 must be implemented in the section "SPECFEM x mesh" below, instead.
dt = 1.5d-2; % Set DT as in parfile.
% nx = 1060; % Set nx as in parfile.
% xmin = -45e3; xmax = 45e3; % Set as in parfile.
% nx = 377; % Set nx as in parfile.
nx = 8060; % Set nx as in parfile.
% xmin = MINX; xmax = MAXX; % Set as in parfile.
xmin = -52e3; xmax = 512e3; % Set as in parfile.
periodise = 0; % Set this if microbaroms do not touch lateral edges of mesh.
% periodise = 1; % Set this if microbaroms do touch lateral edges of mesh.

% Path to file for export (do not forget the "/" at the end).
% EXPORTFILEDIR = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200/';
% EXPORTFILEDIR = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/microbaroms_periodic/';
% EXPORTFILEDIR = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/microbaroms_patch/';
EXPORTFILEDIR = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mb_huge/';

disp(['[',mfilename,'] Preparing meshgrids.']);
%%%%%%%%%%%%%%%%%%%%%%%
% Time, SPECFEM and   %
% forcing.            %
%%%%%%%%%%%%%%%%%%%%%%%
t = 0:dt:MAXTIME+5*dt; % Span a bit more than what is needed.
% if(mod(size(t,2),2)==0)
%   t=[t,t(end)+dt]; % In order to be sure to have odd number of points.
% end
%%%%%%%%%%%%%%%%%%%%%%%
% Space, SPECFEM.     %
%%%%%%%%%%%%%%%%%%%%%%%
x_specfem = linspace(xmin, xmax, nx + 1); % Set x span. This array must reproduce exactly mesh at z=0 (with GLL points).
max_dx_specfem=max(diff(x_specfem)); % Find back dx.
n=2; x_specfem=x_specfem(x_specfem>=MINX-n*max_dx_specfem & x_specfem<=MAXX+n*max_dx_specfem); % Remove useless points on the side, but span n elements more on each side.
NGLL = 5; GLL = fliplr([lglnodes(NGLL-1)/2+0.5]'); GLL = GLL(1:end-1); % Better method, and flexible w.r.t NGLL, but needs function lglnodes from https://fr.mathworks.com/matlabcentral/fileexchange/4775-legende-gauss-lobatto-nodes-and-weights.
% GLL = [0., 0.172673164646011457, 0.5, 0.827326835353988543]; % Ugly method, but exact down to eps since values were extracted from function lglnodes.
x_specfem_with_GLL = [];
for i = 1:length(x_specfem) - 1
  x_specfem_with_GLL = [x_specfem_with_GLL, x_specfem(i) + (x_specfem(i + 1) - x_specfem(i)) * GLL];
end
x_specfem_with_GLL = [x_specfem_with_GLL, x_specfem(end)];
%%%%%%%%%%%%%%%%%%%%%%%
% Space, forcing.     %
%%%%%%%%%%%%%%%%%%%%%%%
% Build space grid which has to be uniform for FFTs.
max_dx_specfem = max(diff(x_specfem_with_GLL));
x=linspace(x_specfem_with_GLL(1),x_specfem_with_GLL(end),floor((x_specfem_with_GLL(end)-x_specfem_with_GLL(1))/max_dx_specfem)+1);
dx=mean(diff(x));
if(periodise)
  % When periodising, we want to force extremities to have the same value. This is done later, below.
  % Because of this, we remove for all FFT computations the last SPECFEM point.
  x = x(1:end-1);
end
%%%%%%%%%%%%%%%%%%%%%%%
% Check steps.        %
%%%%%%%%%%%%%%%%%%%%%%%
meshok=-1;
while(not(ismember(meshok,[0,1])))
  meshok=input(['[',mfilename,'] > SPECFEM grid is xmin=',num2str(xmin),', xmax=',num2str(xmax),', nx=',num2str(nx),'. This has to match parfile. Is that so (0 for no, 1 for yes)? > ']);
end
if(meshok==0)
  error(['[',mfilename,', ERROR] > Mesh was not ok, re-chose parametrisation in script.']);
end
clear('meshok');
disp(['[',mfilename,'] > dt=',num2str(dt),', that is ',num2str(floor(T0/dt)),' iterations per main time period.']);
disp(['[',mfilename,'] > dx=',num2str(dx),', that is ',num2str(floor(L0/dx)),' elements per main spatial period.']);
disp(['[',mfilename,'] > dx=',num2str(dx),', max(dxSPCFM)=',num2str(max_dx_specfem),'. One should have dx>=dxSPCFM to prevent aliasing.']);
% Ask user to verify mesh is ok.
meshok=-1;
% disp(['[',mfilename,'] > Interpolating mesh (final mesh) spans [',num2str(min(xSPCFMwGLL)), ', ', num2str(max(xSPCFMwGLL)), '] m with ',num2str(nx), ' points. This mesh HAS TO match SPECFEM''s mesh.']);
while(not(ismember(meshok,[0,1])))
  meshok=input(['[',mfilename,'] > Is that ok (0 for no, 1 for yes)? > ']);
end
if(meshok==0)
  error(['[',mfilename,', ERROR] Mesh was not ok, re-chose parametrisation in script.']);
end
clear('meshok');

% % Ask user to verify forcing is ok.
% forcok=-1;
% disp(['[',mfilename,'] > Forcing spans [',num2str(MINX), ', ', num2str(MAXX), '] m and [0, ',num2str(MAXTIME),'] s.']);
% while(not(ismember(forcok,[0,1])))
%   forcok=input('  Is that ok (0 for no, 1 for yes)? > ');
% end
% if(forcok==0)
%   error('  Forcing was not ok, re-chose parametrisation in script.');
% end
% clear('forcok');

% dt = T0 / Nt; % Temporal resolution.
% dx = L0 / Nx; % Spatial resolution.
% Ask user to verify if creation steps are ok.
% dtxok=-1;
% dxSPCFM=min(diff(xSPCFM));
% disp(['[',mfilename,'] > dt=',num2str(dt),', dtSPCFM=',num2str(dtSPCFM),'. One should have dt>=dtSPCFM to prevent aliasing.']);
% disp(['[',mfilename,'] > dx=',num2str(dx),', dxSPCFM=',num2str(dxSPCFM),'. One should have dx>=dxSPCFM to prevent aliasing.']);
% while(not(ismember(dtxok,[0,1])))
%   dtxok=input('  Continue (0 for no, 1 for yes)? > ');
% end
% if(dtxok==0)
%   error('  dt was not ok, re-chose parametrisation in script.');
% end
% clear('dtxok');

% Temporal frequency range.
f_max = 0.5 / dt;
df = 1 / (t(end)-t(1));
f = - f_max:df:f_max;
% Spatial frequency range.
k_max = 0.5 / dx;
dk = 1 / (x(end)-x(1));
k = - k_max:dk:k_max;
% Temporal and spatial spaces.
% t = 0:dt:nT0 * T0;
% x = - 0.5 * nL0 * L0:dx:0.5 * nL0 * L0;
% 2D ranges.
disp(['[',mfilename,'] > Meshgrid would be (', num2str(numel(x)), ', ', num2str(numel(t)),'), that is ',sprintf('%.0e',numel(x)*numel(t)),' reals. Too heavy computations might crash your computer. 9e7 was ok, 1e8 chugged a bit, 2e8 chugged hard.']);
proceed=-1;
% disp(['[',mfilename,'] > Interpolating mesh (final mesh) spans [',num2str(min(xSPCFMwGLL)), ', ', num2str(max(xSPCFMwGLL)), '] m with ',num2str(nx), ' points. This mesh HAS TO match SPECFEM''s mesh.']);
while(not(ismember(proceed,[0,1,2])))
  proceed=input(['[',mfilename,'] > Stop (0), or proceed as is (1, at your own risk)? > ']);
end
switch(proceed)
  case 0
    error(['[',mfilename,', ERROR] User chose to stop script.']);
  case 2
    command=['sudo renice -100 -p ', num2str(feature('getpid'))];
    disp(['[',mfilename,'] Executing system command ''',command,''' in order to reduce priority.']);
    system(command);
  otherwise
    disp(['[',mfilename,'] Continuing...']);
end
clear('proceed');

[T, X] = meshgrid(t, x);
[F, K] = meshgrid(f, k);
clear('f', 'k'); % Clear variables as soon as possible, as they are not needed after this point, in order to free memory.

%%%%%%%%%%%%%%%%%%%%%%%
% Frequency range.    %
%%%%%%%%%%%%%%%%%%%%%%%
% disp(['[',mfilename,'] Creating draft of signal in spacetime range and sending it to frequency range.']);
% k_0 = 2 * pi / L0;
% w_0 = 2 * pi / T0;
% % s_displ = A*cos(k_0 * X + w_0 * T + pi/2) + A*cos(k_0 * X - w_0 * T - pi/2); % Interfering monofrequency waves, displacement (pi/2 phase shifts in order not to generate discontinuities).
% draft_displ = 2*A*cos(k_0 * X).*cos(w_0 * T); % Interfering monofrequency waves, displacement.
% % s = A1*sin(k_0 * X - w_0 * T) + A2*sin(-k_0 * X - w_0 * T); % Interfering monofrequency waves, velocity.
% spectrum_displ = fftshift(fft2(draft_displ));

disp(['[',mfilename,'] Creating spectrum in frequency range.']);
spectrum_displ=0*F;
spectrum_displ(abs(abs(F)-1/T0)==min(min(abs(abs(F)-1/T0))) & abs(abs(K)-1/L0)==min(min(abs(abs(K)-1/L0))))=1; % Diracs at pertinent places.
clear('F', 'K'); % Clear variables as soon as possible, as they are not needed after this point, in order to free memory.

spectrum_displ_R = real(spectrum_displ);
clear('spectrum_displ'); % Clear variables as soon as possible, as they are not needed after this point, in order to free memory.
% spectrum_displ_I = imag(spectrum_displ);

%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.          %
%%%%%%%%%%%%%%%%%%%%%%%
% Convolve reals with wide Gaussian, replace imaginaries with random phase.
% sig=1; % Spread of Gaussian mask (should be <2);
% GMask=fspecial('gaussian', (floor(min((w_0/(2*pi))/df,(k_0/(2*pi))/dk))-1)*[1,1],0.5*sig); % Create Gaussian mask which covers f=]0,2/T0[ x k=]0,2/L0[.
% GMask=fspecial('gaussian', [floor((w_0/(2*pi))/df),floor((k_0/(2*pi))/dk)]-1,0.5*sig); % Create Gaussian mask which covers f=]0,2/T0[ x k=]0,2/L0[.
% GMask=fspecial('gaussian', [floor((k_0/(2*pi))/dk),floor((w_0/(2*pi))/df)]-1,0.5*sig); % Create Gaussian mask which covers f=]0,2/T0[ x k=]0,2/L0[.
sig_T = ((1 / T0) / 3) / invspread; % Width of the Gaussian mask in time.
sig_X = ((1 / L0) / 3) / invspread; % Width of the Gaussian mask in space.
disp(['[',mfilename,'] Convolving spectrum in frequency range and adding random phase. sigma_T=',num2str(sig_T),', sigma_X=',num2str(sig_X),'.']);
[Fgm, Kgm] = meshgrid(0:df:2 / T0, 0:dk:2 / L0);
GMask = exp(- ((Fgm - 1 / T0) .^ 2 / (2 * sig_T ^ 2) + (Kgm - 1 / L0) .^ 2 / (2 * sig_X ^ 2)));
GMask = GMask / max(max(GMask)); % Normalise, not to mess with maximum amplitude.
clear('Fgm', 'Kgm'); % Clear variables as soon as possible, as they are not needed after this point, in order to free memory.
% disp(0.5./(size(GMask).*[dk,df])) % Display span of mask in physical quantities (how much around base periods is now considered).
% Option 1: full random phase, relying on real part of spectrum.
RPhase = exp(1j * 2 * pi * rand(size(spectrum_displ_R)));
% RPhase = 1; % No phase.
spectrum_displ_R_convolved = conv2(spectrum_displ_R, GMask, 'same');
clear('GMask', 'spectrum_displ_R'); % Clear variables as soon as possible, as they are not needed after this point, in order to free memory.
spectrum_displ_new = spectrum_displ_R_convolved .* RPhase;
if(not(plot_process))
  clear('spectrum_displ_R_convolved'); % Clear variables as soon as possible, as they are not needed after this point, in order to free memory.
end
clear('RPhase'); % Clear variables as soon as possible, as they are not needed after this point, in order to free memory.
% % Option 2: separate real (convolve) and imaginary (add noise) parts.
% SnR=conv2(SR,GMask,'same');
% SnI=SI+(2*rand(size(S))-1)*0.01*(max(max(SI))-min(min(SI)))*0.1;
% Sn=SnR+1j*SnI;

% % Mirror and convolve imaginary part in frequency.
% spectrum_displ_new_R = real(spectrum_displ_new); % Save real part.
% spectrum_displ_new_I = imag(spectrum_displ_new); % Extract imaginary part.
% spectrum_displ_new_I = spectrum_displ_new_I + conj(fliplr(spectrum_displ_new_I)); % Mirroring imaginary part.
% if(periodise)
%   spectrum_displ_new_I = spectrum_displ_new_I + conj(flipud(spectrum_displ_new_I)); % Mirroring imaginary part.
% end
% spectrum_displ_new = spectrum_displ_new_R + 1j * spectrum_displ_new_I; % To be sent back to spacetime range.
% % Mirror whole spectrum.
% spectrum_displ_new = spectrum_displ_new + conj(fliplr(spectrum_displ_new)); % Mirroring imaginary part.
% if(periodise)
%   spectrum_displ_new = spectrum_displ_new + conj(flipud(spectrum_displ_new)); % Mirroring imaginary part.
% end
% spectrum_displ_new_R = real(spectrum_displ_new); % Save real part.
% spectrum_displ_new_I = imag(spectrum_displ_new); % Extract imaginary part.

if(plot_process)
  spectrum_displ_new_R = real(spectrum_displ_new); % For plotting purposes only.
  spectrum_displ_new_I = imag(spectrum_displ_new); % For plotting purposes only.
end

%%%%%%%%%%%%%%%%%%%%%%%
% Send back to        %
% spacetime range.    %
%%%%%%%%%%%%%%%%%%%%%%%
disp(['[',mfilename,'] Sending back to spacetime range.']);
displ = ifft2(ifftshift(spectrum_displ_new), 'symmetric');
clear('spectrum_displ_new'); % Clear variables as soon as possible, as they are not needed after this point, in order to free memory.
displ = displ - mean(mean(displ)); % Ensure underlying displacement is zero-mean.
% if(0)
%   disp(['[',mfilename,'] Filtering noise.']);
%   % Low-pass in time.
%   fcut = min(20 / T0,0.99*f_max); % Must be < f_max=0.5/dt;
%   for ix = 1:length(x)
%     [low_pass, ~] = custom_filter(t, displ(ix, :), fcut);
%     displ(ix, :) = low_pass;
%   end
%   clear('ix');
%   % Low-pass in space (mainly safeguard).
%   fcut = min(20 / L0,0.99*k_max); % Must be < k_max=0.5/dx;
%   for it = 1:length(t)
%     [low_pass, ~] = custom_filter(x, displ(:, it), fcut);
%     displ(:, it) = low_pass;
%   end
%   clear('it');
%   clear('fcut');
% end
%%%%%%%%%%%%%%%%%%%%%%%
% Plot process.       %
%%%%%%%%%%%%%%%%%%%%%%%
if(plot_process)
  figure();
%   subplot(331); surf(T,X,s_displ,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(t), max(t), min(x), max(x)]);
  subplot(332); surf(F,K,spectrum_displ_R,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(f), max(f), min(k), max(k)]);
%   subplot(333); surf(F,K,S_displ_I,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(f), max(f), min(k), max(k)]);
  subplot(335); surf(F,K,spectrum_displ_R_convolved,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(f), max(f), min(k), max(k)]);
  subplot(337); surf(T,X,displ,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(t), max(t), min(x), max(x)]);
  subplot(338); surf(F,K,spectrum_displ_new_R,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(f), max(f), min(k), max(k)]);
  subplot(339); surf(F,K,spectrum_displ_new_I,'edgecolor','interp','facecolor','interp'); view([0,0,1]); axis([min(f), max(f), min(k), max(k)]);
  pause;
  clear('spectrum_displ_R', 'spectrum_displ_R_convolved', 'spectrum_displ_new_R', 'spectrum_displ_new_I'); % Clear variables as soon as possible, as they are not needed after this point, in order to free memory.
end
%%%%%%%%%%%%%%%%%%%%%%%
% Convert signal in   %
% displacement to     %
% signal in velocity, %
% and adjust          %
% amplitude.          %
%%%%%%%%%%%%%%%%%%%%%%%
disp(['[',mfilename,'] Converting to velocity and adjusting amplitude.']);
% disp(['[',mfilename,'] Displacement amplitude is [',num2str(min(min(s_displ_new))),', ',num2str(max(max(s_displ_new))),']']);
% ampli=max(abs(min(min(s_displ_new))),abs(max(max(s_displ_new)))); % Maximum amplitude of forcing.
% s_displ_new=s_displ_new*A_aimed/ampli; % Rescale.
% disp(['[',mfilename,'] Displacement amplitude is now [',num2str(min(min(s_displ_new))),', ',num2str(max(max(s_displ_new))),']']);
% [s_vel_new,~]=gradient(s_displ_new,t,x);
% disp(['[',mfilename,'] Forcing amplitude in velocity is [',num2str(min(min(s_vel_new))),', ',num2str(max(max(s_vel_new))),']']);
disp(['[',mfilename,'] > Displacement is in [',num2str(min(min(displ))),', ',num2str(max(max(displ))),'].']);
ampli_displ=max(abs(min(min(displ))),abs(max(max(displ))));
displ=displ*A_aimed/ampli_displ; % Scale displacement to aimed wave height.
disp(['[',mfilename,'] > Rescaled displacement is in [',num2str(min(min(displ))),', ',num2str(max(max(displ))),'].']);
[veloc,~]=gradient(displ,t,x);
clear('displ'); % Clear variables as soon as possible, as they are not needed after this point, in order to free memory.
disp(['[',mfilename,'] > Velocity is in [',num2str(min(min(veloc))),', ',num2str(max(max(veloc))),'].']);
% if(periodise)
% %   disp(['[',mfilename,'] > Periodisation: ', num2str(100*mean(abs((veloc(1,:)-mean(veloc([1,end],:),1))./mean(veloc([1,end],:),1))), '% mean change.']);
%   % Make sure one side is equal to the other.
%   veloc(1,:)=mean(veloc([1,end],:),1);
%   veloc(end,:)=veloc(1,:);
% end

%%%%%%%%%%%%%%%%%%%%%%%
% Apodisation.        %
%%%%%%%%%%%%%%%%%%%%%%%
if(activate_apo)
  disp(['[',mfilename,'] Applying apodisation.']);
  n=napoxlr; apox = 0.25 .* (1. - erf((x - MAXX + 0.5 * n * L0) / (0.18 * n * L0))) .* (1 + erf((x - MINX - 0.5 * n * L0) / (0.25 * n * L0))); % Apodisation over n spatial period on each side.
  n=napot0; apot0 = 0.5 .* (1 + erf((t - 0.5 * n * T0) / (0.18 * n * T0))); % Apodisation over n temporal periods at beginning.
  % n = 1.5; apot = 0.5 .* (1. - erf((t - MAXTIME + 0.5 * n * T0) / (0.25 * n * T0))); % Apodisation over n temporal periods at end.
  n=napotend; apot = 0.5 .* (1. - erf((t - MAXTIME + 0.5 * n * T0) / (0.18 * n * T0))); % Apodisation over n temporal periods at end.
  apot = apot0 .* apot;
  apot([1,end]) = 0; % Tweak to make sure forcing starts at 0.
  apox([1,end]) = 0; % Tweak to make sure forcing starts at 0.
  
  if(periodise)
    disp("Deactivating spatial apodisation due to periodisation.");
    apox = 1; % Deactivate apodisation.
  end
  
  veloc = veloc - mean(mean(veloc)); % Make sure velocity is zero-mean in order to ensure underlying displacement is zero-mean.
  veloc = veloc .* apox'.*apot; % Apply apodisation to make sure forcing is zero at beginning and end, and eventually on sides.
  clear('apox', 'apot', 'apot0'); % Clear variables as soon as possible, as they are not needed after this point, in order to free memory.
end

disp(['[',mfilename,'] Verifying zero-mean, zero-start, and zero-end.']);
mmveloc=mean(mean(veloc));
disp(['[',mfilename,'] > Maximum of velocity forcing at t=0     and t=T_max: ',num2str(max(abs(veloc(:,[1,end])))),',']);
disp(['[',mfilename,']                               at x=x_min and x=x_max: ',num2str(max(abs(veloc([1,end],:))')),'.']);
disp(['[',mfilename,'] > Mean of velocity: ',num2str(mmveloc),'.']);
while(mmveloc>1e-18)
  disp(['[',mfilename,'] > Tweaking again.']);
  veloc=veloc-mean(mean(veloc));
  veloc(:,[1,end])=0;
  veloc([1,end],:)=0;
  mmveloc=mean(mean(veloc));
  disp(['[',mfilename,'] > Maximum of velocity forcing at t=0     and t=T_max: ',num2str(max(abs(veloc(:,[1,end])))),',']);
  disp(['[',mfilename,']                               at x=x_min and x=x_max: ',num2str(max(abs(veloc([1,end],:))')),'.']);
  disp(['[',mfilename,'] > Mean of velocity: ',num2str(mmveloc),'.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Interpolate on the       %
% SPECFEM mesh.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['[',mfilename,'] Interpolating on SPECFEM meshgrid.']);
% Here, it is sufficient to configure the next 2 sections (SPECFEM time 
% and SPECFEM mesh).
%%%%%%%%%%%%%%%%%%%%%%%
% Prepare meshgrid    %
% and interpolate.    %
%%%%%%%%%%%%%%%%%%%%%%%
% xSPCFMwGLL=xSPCFMwGLL(xSPCFMwGLL>=MINX-5*dxSPCFM & xSPCFMwGLL<=MAXX+5*dxSPCFM); % Remove useless points to the side, but span a bit more than what is needed.
if(periodise)
  % Add anew the last point to mesh, and make extremities of forcing match.
  x = [x, x_specfem(end)];
  veloc=[veloc;veloc(1,:)];
  [T,X]=meshgrid(t,x); % Update meshgrid.
  clear('x'); % Clear variables as soon as possible, as they are not needed after this point, in order to free memory.
end
[T_specfem, X_specfem] = meshgrid(t, x_specfem_with_GLL); % Meshgrid SPECFEM.
clear('x_specfem_with_GLL', 't'); % Clear variables as soon as possible, as they are not needed after this point, in order to free memory.
veloc_specfem = interp2(T, X, veloc, T_specfem, X_specfem); % Linear interpolation.
clear('T', 'X', 'veloc'); % Clear variables as soon as possible, as they are not needed after this point, in order to free memory.
veloc_specfem(isnan(veloc_specfem)) = 0; % Set zeros instead of NaNs outside of microbarom zone.
%%%%%%%%%%%%%%%%%%%%%%%
% Plot result (for    %
% a visual inspection %
% of the quality of   %
% interpolation).     %
%%%%%%%%%%%%%%%%%%%%%%%
if(plot_process && plot_forcing)
  figure();
  subplot(121); surf(T, X, veloc, 'edgecolor', 'flat', 'facecolor', 'flat'); view([0, 0, 1]); axis([min(t), max(t), min(x), max(x)]); title("Uniform grid");
  subplot(122); surf(T_specfem, X_specfem, veloc_specfem, 'edgecolor', 'flat', 'facecolor', 'flat'); view([0, 0, 1]); axis([min(min(T_specfem)), max(max(T_specfem)), min(min(X_specfem)), max(max(X_specfem))]); title("SPECFEM grid");
end
if(plot_forcing)
  pcolor(T_specfem, X_specfem, veloc_specfem); shading interp; axis([min(min(T_specfem)), max(max(T_specfem)), min(min(X_specfem)), max(max(X_specfem))]);
end
disp(['[',mfilename,'] Velocity forcing on SPECFEM meshgrid is ready.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Export to file.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['[',mfilename,'] Exporting to file.']);
% Here, format shoud be compatible with the reading which is done in
% the subroutine 'prepare_external_forcing' in 'boundary_terms_DG.f90'.
%%%%%%%%%%%%%%%%%%%%%%%
% Test data.          %
%%%%%%%%%%%%%%%%%%%%%%%
% % EXPORTFILEDIR = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_external_forcing/';
% EXPORTFILEDIR = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_EBF/';
% t=0:4e-4:1;
% x_specfem=linspace(-50,50,51);
% GLL=[0, 0.345346329292028/2, 0.5, 1.654653670707976/2]; % UGLY METHOD, I'M SORRY
% x_specfem_with_GLL=[];
% for i=1:length(x_specfem)-1
%   x_specfem_with_GLL=[x_specfem_with_GLL,x_specfem(i)+(x_specfem(i+1)-x_specfem(i))*GLL];
% end
% x_specfem_with_GLL = [x_specfem_with_GLL, x_specfem(end)];
% [T_specfem,X_specfem]=meshgrid(t,x_specfem_with_GLL);
% % MAXTIME=0.5; MINX=-25; MAXX=25;
% MAXTIME=1; MINX=-50; MAXX=50;
% % LAPO=7.5; apox = 0.25 .* (1. - erf((X_specfem - MAXX + 0.5 * LAPO) / (0.25 * LAPO))) .* (1 + erf((X_specfem - MINX - 0.5 * LAPO) / (0.25 * LAPO)));
% apox=1;
% TAPO=0.25; apot0 = 0.5 .* (1 + erf((T_specfem - 0.5 * TAPO) / (0.25 * TAPO))); apot = 0.5 .* (1. - erf((T_specfem - MAXTIME + 0.5 * TAPO) / (0.25 * TAPO))); apot0(:,1)=0;
% % FORCING_INTERP=(abs(XSPCFM)<MAXX).*(TSPCFM<MAXTIME).*sin(XSPCFM/8).*sin(TSPCFM*12.5); % Test forcing function.
% % FORCING_INTERP=(abs(XSPCFM)<MAXX).*(TSPCFM<MAXTIME).*sin(XSPCFM/(0.25*8)).*sin(TSPCFM*4*12.5); % Test forcing function.
% % FORCING_INTERP=(abs(XSPCFM)<MAXX).*(TSPCFM<MAXTIME).*apox.*apot0.*apot; % Test forcing function.
% veloc_specfem=(abs(X_specfem)<=MAXX).*(T_specfem<=MAXTIME).*sin(T_specfem*2*pi/0.25).*apox.*apot0.*apot; % Test forcing function.
% veloc_specfem=A_aimed*sin(2*pi*X_specfem./L0).*sin(2*pi*T_specfem./T0).*apox'.*apot;
%%%%%%%%%%%%%%%%%%%%%%%
% Detect relevant     %
% indices, and print  %
% corresponding       %
% values to file.     %
%%%%%%%%%%%%%%%%%%%%%%%
% Ask user.
exportok=-1;
% disp(['[',mfilename,'] > Interpolating mesh (final mesh) spans [',num2str(min(xSPCFMwGLL)), ', ', num2str(max(xSPCFMwGLL)), '] m with ',num2str(nx), ' points. This mesh HAS TO match SPECFEM''s mesh.']);
while(not(ismember(exportok,[0,1])))
  bytespervalue=25.335164059114234; % [4418901533/178467744 + 4624117614/178467744] / 2. Formerly 39.985880510267151.
  expectedsize=prod(size(veloc_specfem))*bytespervalue;
  disp(['[',mfilename,'] > File would be ~', num2str(expectedsize), ' bytes (',num2str(expectedsize/1000),' kB, ',num2str(expectedsize/1000000),' MB).']);
  exportok=input(['[',mfilename,'] > Export to file (0 for no, 1 for yes)? Path will be confirmed later. > ']);
end
if(exportok==0)
  error(['[',mfilename,', ERROR] > Exporting disabled.']);
else
  disp(['[',mfilename,'] Exporting is ready.']);
  datestr_of_now=datestr(now,'_yymmdd_hhMMss');
  EXPORTFILENAME = [EXPORTFILEDIR, 'external_bottom_forcing',datestr_of_now,'.dat'];
  % Ask user to verify export path is ok.
  pathok=-1;
  disp(['[',mfilename,'] > Export path is ''',EXPORTFILENAME,'''.']);
  while(not(ismember(pathok,[0,1])))
    pathok=input(['[',mfilename,'] > Is that ok (0 for no, 1 for yes)? > ']);
  end
  if(pathok==0)
    error(['[',mfilename,'] > Export path was not ok, re-chose in script.']);
  end
  clear('pathok');
  % Export.
  system(['cp ',mfilename('fullpath'),'.m ',EXPORTFILEDIR,mfilename,datestr_of_now,'.m']);
  disp(['[',mfilename,'] Copied current script to export folder for parameter saving.']);
%   itstop = find(abs(TSPCFM(1, :) - MAXTIME) == min(abs(TSPCFM(1, :) - MAXTIME))) + 1;
%   ixmin = max(find(abs(XSPCFM(:, 1) - MINX) == min(abs(XSPCFM(:, 1) - MINX))) - 1, 1);
%   ixmax = min(find(abs(XSPCFM(:, 1) - MAXX) == min(abs(XSPCFM(:, 1) - MAXX))) + 1, length(XSPCFM(:, 1)));
  f_new = fopen(EXPORTFILENAME, 'w');
%   for it = 1:itstop
%     for ix = ixmin:ixmax
  for it = 1:size(T_specfem,2)
    for ix = 1:size(X_specfem,1)
      fprintf(f_new, '%.5g %.8g %.4g', T_specfem(1, it), X_specfem(ix, 1), veloc_specfem(ix, it));
      fprintf(f_new, "\n");
    end
    if (ismember(it, floor((1:10) * 0.1 * size(T_specfem,2))))
       disp(['[',mfilename,'] > Writing to file (', num2str(ceil(it / size(T_specfem,2) * 100)), ' % complete).']);
    end
  end
  fclose('all');
end
