% Author:        Quentin Brissaud, Léo Martire.
% Description:   TODO.
% Notes:         [Fritts and Alexander, 2003, Brissaud et al., 2016].
%
% Usage:
%   Use this file to plot on seismograms traces exiting from specfem2D in 
%   ASCII filesand a synthetic signal waited at the station in case of 
%   simpler homogeneous medium with a variable density and a constant 
%   velocity.
%   You will choose the directory in which the ASCII files are.
%   Then choose the number of receivers with the velocity of the signal and 
%   finally run this m-file.
%   This program is made for a maximum number of 9 stations.
% with:
%   TODO.
% yields:
%   TODO.

% [Brissaud et al., 2016] Brissaud, Q., Martin, R., Garcia, R. F., & Komatitsch, D. (2016). Finite-difference numerical modelling of gravitoacoustic wave propagation in a windy and attenuating atmosphere. Geophysical Journal International, 206(1), 308–327. https://doi.org/10.1093/gji/ggw121

clear all;
% close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SPCFMloc = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/';
addpath(genpath([SPCFMloc,'utils_new/Atmospheric_Models']));
% factor_err = 100; % factor by which multiply difference in plots.
prefix = 'validation_lns_fk';
% fignames
savefigname = ['compare2FKanalytic'];
savefigpath = ['/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/',prefix,'_info/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rootd=[SPCFMloc,'EXAMPLES/',prefix,'_1.00dx_1.00dt/']; % EXAMPLE path
% rootd=[SPCFMloc,'EXAMPLES/',prefix,'_1.00dx_0.50dt/']; % EXAMPLE path
% rootd=[SPCFMloc,'EXAMPLES/',prefix,'_0.50dx_0.50dt/']; % EXAMPLE path
% OFd = [rootd,'OUTPUT_FILES/'];
OFd = [rootd,'OUTPUT_FILES_isobaric_LNS_st1']; % ISOBARIC CASE
% OFd = [rootd,'OUTPUT_FILES_isothermal_LNS_st1'];
subsample = 0; % Subsample synthetics? See sumsamble_dt below, near T_0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beging program.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pretreatment.
simuType = 1; % 1 for forcing and 2 for attenuation
if(not(OFd(end)=='/')); OFd = [OFd, '/']; end;
% atmmodelfile=[directory,'1D_iso_enhanced.txt'];
atmmodelfile = [OFd, 'input_atmospheric_model.dat'];
parfile = [OFd, 'input_parfile'];
int_file = [OFd, 'input_interfaces'];

% Automatically load parameters from simulation's OUTPUT_FILES directory.
disp([' ']); disp(['[',mfilename,', INFO] Loading simulation parameters.']);
[zmin, zmax] = readExampleFiles_zMinMaxInterfacesFile(int_file);
dt = readExampleFiles_extractParam(parfile, 'DT', 'float');
TYPE_FORCING = readExampleFiles_extractParam(parfile, 'TYPE_FORCING', 'int');
P_0 = readExampleFiles_extractParam(parfile, 'main_spatial_period', 'float'); % main spatial period
T_0 = readExampleFiles_extractParam(parfile, 'main_time_period', 'float'); % main time period
X_0 = readExampleFiles_extractParam(parfile, 'forcing_initial_loc', 'float'); % forcing initial location
t_0 = readExampleFiles_extractParam(parfile, 'forcing_initial_time', 'float'); % forcing initial time
seismotype = readExampleFiles_extractParam(parfile, 'seismotype', 'int');
xmin = readExampleFiles_extractParam(parfile, 'xmin', 'float');
xmax = readExampleFiles_extractParam(parfile, 'xmax', 'float');
ymin = 0; ymax = 0; % Setting this for 2D simulations. For 3D simulations, another call to readExampleFiles_extractParam would be needed.
nx = readExampleFiles_extractParam(parfile, 'nx', 'int');
subsample_dt = T_0/100; % if subsample == 1, wanted dt for subsampling: ~30 pts/period or higher is safe
% Atmospheric model.
switch(readExampleFiles_extractParam(parfile, 'MODEL', 'string'))
  case 'default'
    externalDGAtmosModel = 0;
  case 'external_DG'
    externalDGAtmosModel = 1;
  otherwise
    error('MODEL not implemented.');
end
% Deduce other useful parameters.
dx = (xmax-xmin)/nx;
dy = dx; disp(['[',mfilename,'] Setting dy=dx. Make sure this is the case for your simulation.']);
t0 = 0.0;
% Safeguards.
if(any([readExampleFiles_extractParam(parfile, 'dynamic_viscosity', 'float'), readExampleFiles_extractParam(parfile, 'thermal_conductivity', 'float')]~=0))
  error(['[',mfilename,', ERROR] Simulation was viscous (dynamic_viscosity and/or thermal_conductivity != 0), this case is not implemented yet.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choices
mult_tSpan = 1; % multiplier for time span (used for analytical solution), must be whole numbers
mult_xSpan = 1; % multiplier for x span (used for analytical solution), must be whole numbers
mult_ySpan = 1; % multiplier for y span (used for analytical solution), must be whole numbers, if 2D simulation leave this to 1
% fs   = 1/dt ; %not used anywhere
% fmax = 1; %not used ??
% fmin = -fmax; %not used ??

% Atmospheric Parameters.
disp([' ']); disp(['[',mfilename,', INFO] Preparing atmospheric model.']);
FK_atmModel;

% Loading Results.
disp([' ']); disp(['[',mfilename,', INFO] Loading synthetics.']);
% data = readAndSubsampleSynth(OFd, 2, 'BXZ', 'semv', 1, 1, 0); sig_t = data(:,1)'; sig_v = data(:,2)'; figure(); plot(sig_t, sig_v); % test
[~, signal_type, vname_long, vunit] = seismotypeNames(seismotype, readExampleFiles_extractParam(parfile, 'USE_DISCONTINUOUS_METHOD', 'bool'));
[x_stat, z_stat, y_stat, ~] = loadStations(OFd);
pos_stat = [y_stat, z_stat, x_stat]; % try to do something to stick with rest of code
istattab = 1:numel(x_stat); nstat = length(istattab); % Build istattab.
ystattab = pos_stat(istattab, 1); zstattab = pos_stat(istattab, 2); xstattab = pos_stat(istattab, 3);
if(0)
  % TODO: A switch for this test case.
  % Selection with respect to propagation angle. See [Brissaud et al.,
  % 2016, Section 6.1].
%   NBV             = sqrt(NSQ(1))*10;
  NBV           = sqrt(NSQ(1));
  beta          = acos( 2*pi*(1/T_0)/NBV );
  beta_stat_tot = atan(pos_stat(:, 2) ./ (pos_stat(:,1)-X_0));
  beta_stat     = beta_stat_tot(istattab);
  istattab_ok   = istattab .* (beta_stat<beta)';
end
% Read synthetics' files.
for locStatNumber = 1:nstat
    gloStatNumber = istattab(locStatNumber);
    [tmpDat, nt] = readAndSubsampleSynth(OFd, gloStatNumber, 'BXX', ['sem', signal_type], subsample, subsample_dt, locStatNumber);
    Xtime(locStatNumber,1:nt) = tmpDat(:, 1)';
    Xamp(locStatNumber,1:nt) = tmpDat(:, 2)';
    [tmpDat, nt] = readAndSubsampleSynth(OFd, gloStatNumber, 'BXZ', ['sem', signal_type], subsample, subsample_dt, locStatNumber);
    Ztime(locStatNumber, 1:nt) = tmpDat(:, 1)';
    Zamp(locStatNumber, 1:nt) = tmpDat(:, 2)';
end
clear('tmpDat');
subsampledDt = Ztime(locStatNumber, 2)-Ztime(locStatNumber, 1);
nSamp = nt;
%   dt = subsampledDt; % update dt
disp(['[',mfilename,'] ',num2str(size(Ztime,1)),' synthetics loaded: ', num2str(nSamp), ' samples each.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of the analytical solution.
disp([' ']);
disp(['[',mfilename,', INFO] Building analytical solution.']);
tmax_anal = t_0 + max(zstattab)/SOUNDSPEED(1) + 0.5*T_0 + 1*T_0; % go as far as 2 times the time necessary (roughly) for the signal to pass last station (t_0 + travel time + half period), + full period for safety
% dt_anal   = T_0 / (2*100); % ensure fmax=1/dt >> 2*f0=2/T_0 by setting 1/dt = (2*2.5)/T_0
dt_anal   = min(diff(Ztime(1,:))); % Choose a dt corresponding to the dt of the simulation.
dx_anal = dx;
dy_anal = dy;
% syn = zeros(nstat,length(Ztime(nstat,:))) ;
% synX = zeros(nstat,length(Ztime(nstat,:))) ;
% synZ = zeros(nstat,length(Ztime(nstat,:))) ;
% NFFT2 = 2^nextpow2(length(Ztime(nstat,:))); 
% syn = zeros(nstat, nSamp); %unused?
% synX = zeros(nstat, nSamp); %unused?
% synZ = zeros(nstat, nSamp); %unused?
% NFFT2 = 2^nextpow2(nSamp)*extent_time;
% NFFT2 = 2^nextpow2(nSamp)*mult_tSpan;
NFFT2 = 2^nextpow2(ceil(tmax_anal/dt_anal))*mult_tSpan;
% dt=256*dt; % ?????
% Fourier space
NFFT1 = 2^nextpow2((xmax-xmin)/dx_anal)*mult_xSpan;
NFFT3 = 2^nextpow2((ymax-ymin)/dy_anal)*mult_ySpan;
disp(['[',mfilename,'] FFT3D matrix size: ',sprintf('%.1e', NFFT1*NFFT2*NFFT3),' reals.']);
disp(['[',mfilename,'] Usually, 9e7 is ok, 1e8 chugs a bit, 2e8 chugs HARD.']);
k = zeros(NFFT1,NFFT2,NFFT3);
% t = zeros(1,NFFT2) ;
% t = dt * (0:1:NFFT2-1);
% t = subsampledDt * (0:1:NFFT2-1);
t = dt_anal * (0:1:NFFT2-1);
% maxwind = max(t);
x = dx_anal * (0:1:NFFT1-1) + xmin;
y = dy_anal * (0:1:NFFT3-1) + ymin;
% IN THIS SENSE? -> NO!
%omega = 2.0*pi()*(1.0/(dt*NFFT2))*[ [0.0-[NFFT2/2-1:-1:1]] [0:1:NFFT2/2]];
%kx = 2.0*pi()*(1.0/(dx*NFFT1))*[[0.0-[NFFT1/2-1:-1:1]] [0:1:NFFT1/2] ];
% Define the starting Matrix
%     tendsig=1500.0;
[X, Y, T] = meshgrid(x, y, t);
% [X, T] = meshgrid(x, t);
%     Xp = X+2500;
%     Yp = Y+2500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Brissaud et al., 2016, Section 5.1]: "Calculation of the forcing signal for the whole time domain along the forcing boundary or at the point source."
switch(TYPE_FORCING)
  % See expressions in 'boundary_terms_DG.f90'.
  case(1)
%     Mo = 2 * (   exp(-((T-(t_0-T_0/4)-t0)/(T_0/4)).^2) ...
%                - exp(-((T-(t_0+T_0/4)-t0)/(T_0/4)).^2) ); % not really pretty
    Mo = -0.5*(2.*(pi/T_0)^2) * (T-t_0) .* exp(-((T-t_0)*pi/T_0).^2); % neary equivalent, and an actual Gaussian derivative
    displayMo = squeeze(Mo(1, 1, :));
  case(2)
    Mo = 2 * (   exp(-((X-(X_0-P_0/4))/(P_0/4)).^2) ...
               - exp(-((X-(X_0+P_0/4))/(P_0/4)).^2) ) .* ...
             (   exp(-((T-(t_0-T_0/4)-t0)/(T_0/4)).^2) ...
               - exp(-((T-(t_0+T_0/4)-t0)/(T_0/4)).^2) ); %...
    %              .* exp(-(sqrt((X-xo).^2+(Y-yo).^2)/(lambdo/4)).^2);
  otherwise
    error(['[', mfilename, ',ERROR] Bad TYPE_FORCING.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Brissaud et al., 2016, Section 5.1]: "Calculation of the 3D (or 2D) Fourier transform (spatial and time transformations) of that function."
TFMo = fftn(Mo);
%%%%
%%%% IMPORTANT CORRECTION
%%%%
% fft2 of matlab perform the projection on the function basis of the type:
% exp[i(Kx*X+Ky*Y)] whereas we use exp[i(Kx*X-omega*t)] !!!
% so:
% * Ky = -omega
% * and positive frequencies are in the second part of the fft table
% this is corrected by inserting minus sign in the expression below and by
% using positive frequencies of the second part of the fft table to ensure
% symmetry of the fft (see below)
% the formulas themselves (for dispersion relation) are not modified
%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of wavenumbers kx = 2π/λx and ky = 2π/λy for all spatial wavelengths λ_{x,y} [Brissaud et al., 2016, Section 5.1].
% omega = 2.0*pi()*(1.0/(dt*NFFT2))*[[0:1:NFFT2/2] [0.0-[NFFT2/2-1:-1:1]]];
omega = 2.0*pi()*(1.0/(dt_anal*NFFT2))*[[0:1:NFFT2/2] [0.0-[NFFT2/2-1:-1:1]]];
kx =    2.0*pi()*(1.0/(dx     *NFFT1))*[[0:1:NFFT1/2] [0.0-[NFFT1/2-1:-1:1]]];
ky =    2.0*pi()*(1.0/(dy     *NFFT3))*[[0:1:NFFT3/2] [0.0-[NFFT3/2-1:-1:1]]];
[KX, KY, Omega] = meshgrid(kx, ky, omega);
Nsqtab = NSQ(1, 1) + 0.0*Omega; % Assuming constant Nsq
% KX = Omega/velocity(1);
% see occhipinti 2008 for analytical solution (appendix A1)
onestab = 0.0*Nsqtab + 1.0;
Omega_intrinsic = Omega - wind_x*KX - wind_y*KY;

% Calculation of kz from dispersion relations for all wavenumbers kx, ky and time frequencies [Brissaud et al., 2016, Appendix B].
% [Brissaud et al., 2016, Eq. (B4)] ?? See brouillons 190412
% [Fritts and Alexander, 2003, Eq. (22)] ?? Not  ( 1 - Nsqtab./(Omega.^2)) ?? See brouillons 190412
% Base formula:
% KZ = sqrt(Nsqtab.*(KX.*KX)./((Omega-(wi)*KX).*(Omega-(wi)*KX))-(KX.*KX) + (1 - Nsqtab./(Omega.*Omega) ).*((Omega-(wi)*KX).*(Omega-(wi)*KX))/(SOUNDSPEED(1)^2) );
% KZ = sqrt(Nsqtab.*(KX.*KX)./((Omega-(wind_x)*KX).*(Omega-(wind_x)*KX))-(KX.*KX) + (1 - Nsqtab./(Omega.*Omega) ).*((Omega-(wind_x)*KX).*(Omega-(wind_x)*KX))/(SOUNDSPEED(1)^2) );
% Factorised base formula.
% KZ = sqrt(   (Nsqtab.*KX.^2) ./ ((Omega-wind_x*KX).^2) ...
%            - KX.^2 ...
%            + ( 1 - Nsqtab./(Omega.^2)) ...
%              .* ((Omega-wind_x*KX)./SOUNDSPEED(1)).^2      );
% Factorised base formula + generalisation to stick to Fritts and Alexander.
% KZ = sqrt(   Nsqtab .* (KX.^2+KY.^2) ./ (Omega_intrinsic.^2) ...
%            - (KX.^2+KY.^2) ...
%            - onestab/(4*H^2) ...
%            + ( 1 - Nsqtab./(Omega.^2)) ... % THIS TERM INCREASES LOCAL PROPAGATION SPEED
%              .* (Omega_intrinsic ./ SOUNDSPEED(1)).^2            );
% Full Fritts and Alexander (22) with f=0.
% KZ = sqrt(   (KX.^2+KY.^2) .* (Nsqtab./(Omega_intrinsic.^2) - 1) ...
%            - onestab/(4*H^2) + (Omega_intrinsic ./ SOUNDSPEED(1)).^2 );
% Full Fritts and Alexander (22) with f=0 and no H term.
% KZ = sqrt(   (KX.^2+KY.^2) .* (Nsqtab./(Omega_intrinsic.^2) - 1) ...
%            + (Omega_intrinsic ./ SOUNDSPEED(1)).^2 );
% Full Fritts and Alexander (24) with f=0.
% KZ = sqrt(   (KX.^2+KY.^2) .* (Nsqtab./(Omega_intrinsic.^2) - 1) ...
%            - onestab/(4*H^2)                                         );
% Old formula.
% KZ = sqrt( Nsqtab.*(KX.*KX + KY.*KY)./(Omega_intrinsic.*Omega_intrinsic)-(KX.*KX + KY.*KY) - onestab/(4*H*H) + (Omega_intrinsic.*Omega_intrinsic)/(SOUNDSPEED(1)^2) );
% Test new formula
if(externalDGAtmosModel)
  error('no formula for external atmos model');
else
  if(USE_ISOTHERMAL_MODEL)
    [KZ, corrFact] = FK_KZ_isothermal(Omega_intrinsic, Omega, KX, SOUNDSPEED, H, GRA, GAM, wind_x);
  else
    [KZ, corrFact] = FK_KZ_isobaric(Omega_intrinsic, Omega, KX, SOUNDSPEED);
  end
end
indNanKZ = find(isnan(KZ));
indInfKZ = find(isinf(KZ));
indImagKZ = find(imag(KZ)<0);
disp(['[',mfilename,'] Number of NaNs in KZ: ',num2str(numel(indNanKZ)),'.']);
if(numel(indNanKZ))
  disp(['[',mfilename,']   Setting those NaNs to zeros.']);
  KZ(indNanKZ) = 0.0;
end
disp(['[',mfilename,'] Number of Infs in KZ: ',num2str(numel(indInfKZ)),'.']);
if(numel(indInfKZ))
  disp(['[',mfilename,']   Setting those Infs to zeros.']);
  KZ(indInfKZ) = 0.0;
end
disp(['[',mfilename,'] Number of KZ such that Im(KZ)<0: ',num2str(numel(indImagKZ)),'.']);
if(numel(indImagKZ))
  disp(['[',mfilename,']   Imaginary part of KZ should be positive in order to attenuate the signal. Setting those to their conjugate (only flips the sign of the imaginary part).']);
  KZ(indImagKZ) = conj(KZ(indImagKZ));
end
% real(KZ) should be positive for positive frequencies and negative for
% negative frequencies in order to shift signal in positive times
% restore the sign of KZ depending on Omega-wi*KX
%     KZnew=real(KZ).*sign((Omega-wind_x*KX)).*sign(KX)+1i*imag(KZ);
% !!! Why KZ should have a sign opposite to Omega for GW NOT UNDERSTOOD !!!
% => because vg perpendicular to Vphi ?
KZnew = 0.0 - real(KZ).*sign(Omega_intrinsic) + 1i*imag(KZ);
KZ = KZnew;

%%%%%%%%%%%%%%%%%%%%%%%
% CORRECTING FACTORS. %
%%%%%%%%%%%%%%%%%%%%%%%
switch corrFact
  case 1
    % Correcting factor 1 for isothermal case (use it with the Maple script configured with the isothermal decay = z/H).
    KZ = KZ -1j*(1/H + 2/(SOUNDSPEED^2));
    disp(['[] Used a correcting factor: kz = kz - i*( 1/H+2/(c^2) )']);
    % 1/H comes from the fact that isothDecay=1/H (developing K \cdot X makes this 1/H go outside the i and finds back its place). TBH I do not remember how I found the second factor, sheer luck is not excluded.
  case 2
    % Correcting factor 2 for isothermal case (use it with the Maple script configured without the isothermal decay).
    KZ = KZ -1j*(8.589e-05);
    disp(['[] Used a correcting factor: kz = kz - i*8.589e-5']);
    % Cannot explain why this value.
  case 3
    % Correcting factor 3 for isothermal case (use it with the Maple script configured with isothermal decay = z/(2H)).
%     KZ = KZ -1j*(6.449e-05);
    KZ = KZ -1j*(1/(2*H) + 2/(SOUNDSPEED^2));
    disp(['[] Used a correcting factor: kz = kz - i*( 1/(2H)+2/(c^2) )']);
    % Same rationale as case 1.
end
% KZ = KZ -1j*(2/(SOUNDSPEED^2)); % Test correcting factor for isothermal case.
% KZ = KZ -1j*(1.1643e+05); % Test correcting factor for isothermal case.
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% restore negative frequencies from postive ones:
%      KZ(:,NFFT2/2+2:NFFT2)=0.0-real(KZ(:,NFFT2/2:-1:2))+1i*imag(KZ(:,NFFT2/2:-1:2));
% restore negative frequencies from postive ones:
% Corrected for "%%%%%%    IMPORTANT CORRECTION   %%%%"
%      KZ(:,NFFT2/2:-1:2)=0.0-real(KZ(:,NFFT2/2+2:NFFT2))+1i*imag(KZ(:,NFFT2/2+2:NFFT2));
% Above not taken into account for "filt", because only positive computations 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Brissaud et al., 2016, Section 5.1]: "Multiplication, in the Fourier
% domain, of the forcing function with a complex filter based on the
% representation of the solution in the case of an harmonic source or
% forcing term (see Appendix B for more details)."

%         z_station = zstattab(istat);
%         shifttime=0.0-real(KZ*z_station)./Omega;
%         KZ(1,1,:) = 0;
% remove the waves moving with a shifttime < 0
%         filt=exp(i*(KZ*z_station));
%       filt=exp(i*(KZ*z_station)).*(shifttime>=0);
%         filt(1,1,:)=0.0;
% remove waves assumed to propagate downward
%       indw=find((real(KZ)>0.0));
%      indw=find((real(KZ)<0.0));
%       filt(indw)=0.0;
      % try to remove aliasing
%         ind1=find(shifttime>(maxwind/2));
%        ind=find(abs(shifttime)>(maxwind/2));
%        filt(ind1)=exp(-abs(shifttime(ind1))/(maxwind/2)).*filt(ind1);
%         Mz=1*exp(z_station/(2*H))*ifftn((filt.*TFMo));
% compute X (horizontal component) assuming UX(f)=-(KZ/KX)*UZ
% equations 2.8 and 2.16 of Nappo, 2002 (not affected by wind)
% but frequencies close to (wi*KX-Omega)=0 create very large KZ...
%         TFMx=TFMo.*(KZ./KX).*filt;
%         ind1=find(isnan(TFMx));
%         ind2=find(isinf(TFMx));
%         TFMx(ind1)=0.0+i*0.0;
%         TFMX(ind2)=0.0+i*0.0;
%         Mx=0.0-exp(z_station/(2*H))*ifft2(TFMx);
%        xstattab_analytic = size(xstattab);
%        ystattab_analytic = size(ystattab);
synf = [];
for locStatNumber = 1:nstat
  gloStatNumber = istattab(locStatNumber);
  x_station = xstattab(gloStatNumber);
  z_station = zstattab(gloStatNumber);

  % MAYBE THIS FILTER TO TO BRING BACK TO ACTUAL TIME SERIES, would make sense
  %filt = exp(1i*(KZ*z_station));
  if(externalDGAtmosModel)
    error('no filter implemented for external atmos model');
  else
    if(USE_ISOTHERMAL_MODEL)
      % isothermal
      filt = exp(1i*(KX*x_station + KZ*z_station));
%           filt = exp(1i*(KX*x_station + KZ*z_station) + z_station/H);
%           filt = exp(1i*(KX*x_station + KZ*z_station) + z_station/(2*H));
%           filt = exp(1i*(KX*x_station + KZ*z_station) - 2*z_station/H);
%           filt = exp(1i*(KX*x_station + KZ*z_station) - z_station/(2*H));
%           filt = exp(1i*(KX*x_station + KZ*z_station) + z_station/H + 2*z_station/(SOUNDSPEED^2));
%           filt = exp(1i*(KX*x_station + KZ*z_station) - z_station/H);
    else
      % isobaric
      filt = exp(1i*(KX*x_station + KZ*z_station));
    end
  end
%       filt = exp(1i*(KX*x_station + KZ*z_station));
%       filt = exp(1i*(KX*x_station + KZ*z_station) + z_station/H);
%       shifttime=0.0-real(KZ*z_station)./Omega;
%       filt = exp(1i*KZ*z_station).*(shifttime>=0);
  % what
%       filt(1, 1, :) = 0.0; % why
  % remove waves assumed to propagate downward
%       indw=find((real(KZ)>0.0));
%       indw=find((real(KZ)<0.0));
%       filt(indw)=0.0;
%       Mz = exp(-z_station/(2*H)) * ifftn(filt.*TFMo);
%       Mz = exp(z_station/(2*H)) * ifftn(filt.*TFMo);
%       Mz = exp(z_station/(H)) * ifftn(filt.*TFMo);
%       Mz = 0.5*exp(z_station/(2*H)) * ifftn(filt.*TFMo);
  %Mz = exp(z_station/69109) * ifftn(filt.*TFMo);
%       Mz = exp(z_station/68719) * ifftn(filt.*TFMo);
  Mz = ifftn(filt.*TFMo);
%       Mz = ifftn(filt.*TFMo) * exp(2*z_station/SOUNDSPEED);
%           xcoord(istat) = xstattab(istat);
  ix = round((xstattab(gloStatNumber)-xmin)/dx) + 1; % first guess for x
  iy = round((ystattab(gloStatNumber)-zmin)/dy) + 1; % first guess for y
%           iy=round((ystattab(istat)-ymin)/dy) + 1;
  % GET X LOCATION OF STATION
%       ACHECK = xstattab(gloStatNumber); positionInMeshgrid = X(1, ix, 1); positionInMeshgrid_p1 = X(1, ix+1, 1); positionInMeshgrid_m1 = X(1, ix-1, 1);
  toCheck = xstattab(gloStatNumber); positionInMeshgrid = x(ix); positionInMeshgrid_p1 = x(ix+1); positionInMeshgrid_m1 = x(ix-1);
  if(abs(positionInMeshgrid-toCheck) > abs(positionInMeshgrid_p1-toCheck))
    ix = ix + 1; % if real x is closer to ix+1 than to ix, choose this ix+1
  elseif(abs(positionInMeshgrid-toCheck) > abs(positionInMeshgrid_m1-toCheck));
    ix = ix - 1; % or if real x is closer to ix-1 than to ix, choose this ix-1
  end  
  xstattab_analytic(gloStatNumber) = x(ix); %X(1, ix, 1);

  % GET Y LOCATION OF STATION
  if(numel(y)>1)
    % This if is to prevent trying to find another possible y when
    % simulation is 2D (because in the 2D case, the array 'y' only
    % contains one element, 0).
%       ACHECK = ystattab(gloStatNumber); positionInMeshgrid = Y(iy, 1, 1); positionInMeshgrid_p1 = Y(iy+1, 1, 1); positionInMeshgrid_m1 = Y(iy-1, 1, 1);
    toCheck = ystattab(gloStatNumber); positionInMeshgrid = y(iy); positionInMeshgrid_p1 = y(iy+1); positionInMeshgrid_m1 = y(iy-1);
    if(abs(positionInMeshgrid-toCheck) > abs(positionInMeshgrid_p1-toCheck))
      iy = iy + 1; % if real y is closer to iy+1 than to iy, choose this iy+1
    elseif(abs(positionInMeshgrid-toCheck) > abs(positionInMeshgrid_m1-toCheck))
      iy = iy - 1; % or if real y is closer to iy-1 than to iy, choose this iy-1
    end
  end
  ystattab_analytic(gloStatNumber) = y(iy); %Y(iy,1,1)

  % If out of mesh, bring back to closest mesh point.
  if(ix<=0); ix = 1; end
  if(iy<=0); iy = 1; end

  % Take real part as analytic solution.
  synf(locStatNumber, :) = real(Mz(iy, ix, :));
end

% test to see what forms has the analytic solution.
disp([' ']);
disp(['[',mfilename,', INFO] Checking analytical solution.']);
figcheckhandle=figure();
% colours=jet(nstat);
% colours=prism(nstat);
colours=winter(nstat);
if(externalDGAtmosModel)
  error('no th value for external atmos model');
else
  if(USE_ISOTHERMAL_MODEL)
    th_vz_v0 = exp( zstattab(end)/(2*H) );
    
    % Look at agreement with exponential decay.
    factor = 1; facttxt = '1';
%     factor = exp(2*zstattab/(SOUNDSPEED^2)); facttxt = 'exp(2*z/c$^2$)';
%     factor = exp(1.446e-5 * zstattab); facttxt = 'exp(1.446e-5 * z)';
    isothDecay = [zstattab, factor.* peak2peak(synf(:,:),2)./peak2peak(synf(1,:)), exp( zstattab/(2*H) ), peak2peak(Zamp(:,:),2)/peak2peak(Zamp(1,:))];
    x=zstattab; y=isothDecay(:,3)./isothDecay(:,2); lel = fit(x,y,'exp1');
    if(log10(abs(lel.b))>=-6)
      % If a non-negligible multiplying factor fit is found, plot and display.
      lel
      sprintf('%.16f',lel.b)
      figure();
      plot(zstattab, isothDecay(:,3), 'displayname', ['exponential decay']); hold on;
      plot(zstattab, isothDecay(:,4), 'displayname', ['synthetics']); hold on;
      plot(zstattab, isothDecay(:,2), 'displayname', ['analytic $\times$',facttxt]); hold on;
      legend();
    end
  else
    th_vz_v0 = 1;
  end
end
figure(figcheckhandle);
for i = 1:nstat;
  plot(t,synf(i,:),'displayname',num2str(zstattab(istattab(i))),'color',colours(i,:)); hold on;
end
title({['This analytical solution: $v(z=z_{max})/v(z=0)$=',num2str(peak2peak(synf(end,:))/peak2peak(synf(1,:))),'.'],['Theoretical value: =$e^{z/(2H)}$=',num2str(th_vz_v0),'.'],['Current simulation: =',num2str(peak2peak(Zamp(end,:))/peak2peak(Zamp(1,:))),'.']},'fontsize',14);
legend();
prettyAxes(figcheckhandle);
isThisOk=-1;
while(not(ismember(isThisOk,[0,1])))
  isThisOk=input(['[',mfilename,'] Does this analytical solution look ok (0 for no, 1 for yes)? > ']);
end
if(not(isThisOk))
  error('stop kek');
end

% TEST ONLY
% synf = Zamp; 
% test only !!!!
% synf_X=synf;
% clear TFMx;
% clear Mx ;
% clear Mz ;
timef = t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finalise.
factor_err = 100;
allStations = 1:nstat;
stationsToPlot = allStations;
nstatToPlot = numel(stationsToPlot);

disp([' ']);
disp(['[',mfilename,'] Computing errors.']);
FK_computeError;

disp([' ']);
disp(['[',mfilename,'] Plotting everything.']);
FK_plot;