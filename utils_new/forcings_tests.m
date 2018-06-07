% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         This script can be used to provide illustration of the various bottom forcings implemented in the DG extension.

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

set(groot,'defaultSurfaceEdgeColor','none');

%% Microbarom (smart way).
if(0==1)
  % Parameters.
  N=50; % Points per period.
  T0=7; % Temporal period.
  nT0=5; % Number of temporal periods to represent.
  dt=T0/N; % Temporal resolution.
  L0=200; % Spatial period.
  nL0=5; % Number of spatial periods to represent.
  dx=L0/N; % Spatial resolution.

  % Temporal frequency range.
  f_max=0.5/dt;
  df=1/(nT0*T0);
  f=0:df:f_max;
  % Spatial frequency range.
  k_max=0.5/dx;
  dk=1/(nL0*L0);
  % f_X=0:df_X:f_X_max;
  k=-k_max:dk:k_max;
  % Temporal and spatial spaces.
  t=0:dt:nT0*T0;
  x=0:dx:nL0*L0;

  % 2D ranges.
  [T,X]=meshgrid(t,x);
  [F,K]=meshgrid(f,k);

  % Spectra.
  S=zeros(size(F));
  % S(find(abs(f_X-1/X)==min(abs(f_X-1/X))),find(f_T-1/T==min(abs(f_T-1/T))))=1; % Dirac(1/T, 1/X).
  % S(find(abs(f_X)==min(abs(f_X))),find(f_T-1/T==min(abs(f_T-1/T))))=1; % Dirac(1/T, 0).
  % S(find(abs(f_X-1/X)==min(abs(f_X-1/X))),find(f_T-1/T==min(abs(f_T-1/T))))=1; S(find(abs(f_X+1/X)==min(abs(f_X+1/X))),find(f_T-1/T==min(abs(f_T-1/T))))=1; % Dirac(1/T, \pm1/X).

  % spread=1.75;FT0=1/T;sigma_T=(FT0/3)/spread;FX0=1/X;sigma_X=(FX0/3)/spread; S=exp( -( (FT-FT0).^2/(2*sigma_T^2) + (FX-FX0).^2/(2*sigma_X^2) ) ); % Gaussian @(1/T,1/X).
  spread=4;FT0=1/T0;sig_T=(FT0/3)/spread;FX0=1/L0;sig_X=(FX0/3)/spread; S=exp( -( (F-FT0).^2/(2*sig_T^2) + (K-FX0).^2/(2*sig_X^2) ) )+exp( -( (F-FT0).^2/(2*sig_T^2) + (K+FX0).^2/(2*sig_X^2) ) ); % Gaussians @(1/T,1/X) & (1/T, -1/X).

  figure(1e5);
  subplot(121);
  surf(F,K,S,'edgecolor','interp','facecolor','interp'); xlim([0, 2/T0]); ylim([-2/L0, 2/L0]); view([0,0,1]);

  % Build phase.
  PHASE_p_p=exp(1j*2*pi*rand(floor(length(k)/2)+1,length(f))); % Random Phase_(+, +).
  % PHASE_p_m=exp(1j*2*pi*rand(floor(length(f_X)/2),length(f_T))); % Random Phase_(+, -).
  % PHASE_p_m=flipud(PHASE_p_p(2:end,:)); % Phase_(+, -) = Phase_(+, +).
  PHASE_p_m=-flipud(PHASE_p_p(2:end,:)); % Phase_(+, -) = -Phase_(+, +). Corresponds to +180° w.r.t. (+, +).
  PHASE=[PHASE_p_p;PHASE_p_m]; % Phase_(+, .) is the concatenation.
  % PHASE(find(abs(f_X)==min(abs(f_X))),:)=0; % Set no Phase at k=0.
  % PHASE(:,find(abs(f_T)==min(abs(f_T))))=0; % Set no Phase at w=0.

  % Add phase.
  S=-S.*PHASE;

  % Mirror around f_T=0 (frequency increasing from left to right).
  S=[conj(fliplr(S(:,2:end))),S];

  % fftshift and ifft.
  s=ifft2(fftshift(S));
  % s=ifft2(S);

  % Representation.
  subplot(122);
  surf(T,X,real(s),'edgecolor','interp','facecolor','interp');view([0,0,1]);axis([t(1),t(end),x(1),x(end)]);

  % Verification.
  sn=real(s);
  SN=fftshift(fft2(sn));
  disp(max(max(abs(real(S)-real(SN)))));
  disp(max(max(abs(imag(S)-imag(SN)))));
end

%% Microbarom (brutal way).
if(0==1)
  perio=7;
  lambdo=200;
  MICROBAROM_MAXTIME=8*perio;
  MICROBAROM_RANGE=15e3;
  dt=1.5d-2;
  % t=max(0,MICROBAROM_MAXTIME-4*perio):dt:MICROBAROM_MAXTIME+perio/4;
  t=0:dt:3*perio;
  x=linspace(MICROBAROM_RANGE-15*lambdo,MICROBAROM_RANGE+lambdo/4,1000);
  % x=linspace(-MICROBAROM_RANGE-lambo/4,MICROBAROM_RANGE+lambo/4,1000);

  % Build spatial and temporal phases based on normal random walks.
  Phix=zeros(length(t),1); Phix(1)=0;
  Phit=zeros(length(t),1); Phit(1)=0;
  for i=2:length(t)
    % Phase must at most shift pi/2 over one period (otherwise it would be
    % too chaotic). Thus, if we denote the step s, and if the random walk is
    % the most unlucky possible (is all steps done in only one direction),
    % (period/dt)*s must be < pi/2. In other words, d must be <
    % pi*dt/(2*period). If the walk is normal, we thus choose n*sigma = 
    % pi*dt/(2*period), that is sigma = pi*dt/(n*2*period). n=1 is generally
    % a good choice.
    n=0.2;
    sigma=pi*dt/(n*2*perio);
    Phix(i)=Phix(i-1)+normrnd(0,sigma);
    Phit(i)=Phit(i-1)+normrnd(0,sigma);
  end

  % Apodisation
  % n=10; apox=0.25.*(1.-erf((x-MICROBAROM_RANGE+n*lambdo)/(0.5*n*lambdo))).*(1+erf((x+MICROBAROM_RANGE-n*lambdo)/(0.5*n*lambdo)));
  n=10; apox=(1-((-abs(x)+MICROBAROM_RANGE)/(n*lambdo)-1).^2 .*(abs(x)>MICROBAROM_RANGE-n*lambdo)).*(abs(x)<MICROBAROM_RANGE);
  % n=1; apot=0.5.*(1.-erf((t-MICROBAROM_MAXTIME+n*perio)/(0.5*n*perio)));
  n=2; apot=(t<MICROBAROM_MAXTIME-n*perio)+ (1-((t-MICROBAROM_MAXTIME)/(n*perio)+1).^2) .*(t>MICROBAROM_MAXTIME-n*perio).*(t<MICROBAROM_MAXTIME);
  % n=1; apot0=0.5*(t/(n/2*perio)).^2 .* (t<n/2*perio) + (0.5+0.5*(1-(2-t/(n/2*perio)).^2)) .*(t>n/2*perio).*(t<n*perio);
  n=1; apot0=2*(t/(n*perio)).^2 .* (t<n/2*perio) + ((4*t)/(n*perio)-2*(t/(n^2*perio)).^2-1) .*(t>n/2*perio).*(t<n*perio) + (t>n*perio);
  apot=apot.*apot0;

  [L0,T0]=meshgrid(x,t);

  Fx=sin(2.*pi*L0/lambdo+Phix);
  Ft=sin(2.*pi*T0/perio+Phit);
  Fx=Fx.*(abs(L0)<MICROBAROM_RANGE);
  Ft=Ft.*(T0<MICROBAROM_MAXTIME);
  Fx=Fx.*apox;
  Ft=Ft.*apot';
  S=Fx.*Ft;

  figure();
  surf(L0,T0,S,'EdgeColor','none');
  view([0,0,1]);
  colormap('jet');
  xlim([min(x),max(x)]);
  ylim([min(t),max(t)]);
  xlabel('$x$');
  ylabel('$t$');
end