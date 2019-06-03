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
clc ;
format shortg
% clear DIFF
% clear synf

set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'defaultTextFontSize', 20); set(0, 'defaultAxesFontSize', 20);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
SPCFMloc = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/';
addpath([SPCFMloc,'utils_new/Atmospheric_Models']); % extract_atmos_model
addpath([SPCFMloc,'utils_new/tools']); % readExampleFiles_extractParam, seismotypeNames

factor_err = 100; % factor by which multiply difference in plots.

rootd=[SPCFMloc,'EXAMPLES/validation_lns_fk/']; % EXAMPLE path
% OFd = [rootd,'OUTPUT_FILES/'];
% OFd = [rootd,'OUTPUT_FILES_long/'];
% OFd = [rootd,'OUTPUT_FILES_1904171716_redone_velocity/'];
% OFd = [rootd,'OUTPUT_FILES_1904171808_vel_isobaric/'];
% OFd = [rootd,'OUTPUT_FILES_1904171832_vel_isobaric_LNS/'];
% OFd = [rootd,'OUTPUT_FILES_isobaric_LNS_190603_st2_morestations_corrected']; % this is p', we want v'
OFd = [rootd,'OUTPUT_FILES_isobaric_LNS_190603_st1'];
% OFd = [rootd,'OUTPUT_FILES_isobaric_FNS_190603_st1'];

data_test0_or_readRun1 = 1;
subsample = 1; % subsample synthetics? see sumsamble_dt below, near T_0

%%%%%%%%%
% Display
subtitle = strcat('Forcing with the whole surface moving at the bottom of a fluid medium');
subtitle2 = strcat([' with a variable density profile and a sound speed equal to ',num2str(639.52),' m.s^{-1}']) ;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory of semd files
% directory = strcat('/home/garcia/SATGRAVI/Brissaud/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4/') ;
% directory = strcat('/home/garcia/SATGRAVI/Brissaud/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_LDDRK/') ;
% directory2 = strcat('/home/garcia/SATGRAVI/Brissaud/test_solution_analytique_gravi/Results.semd/') ;
% directory = strcat('/home/garcia/SATGRAVI/Brissaud/files_wind_quentin/gravi_forcing_nowind/') ; 
%directory = strcat('/home/garcia/SATGRAVI/Brissaud/files_wind_quentin/forcing_gravi_1Diso_windcte_noatten/') ; 
%  directory = strcat('/home/garcia/SATGRAVI/Brissaud/files_wind_quentin/FG1_1Diso_wind_noatten/') ; 
% directory = strcat('/home/garcia/SATGRAVI/Brissaud/files_wind_quentin/FG1_1Diso_nowind_noatten/') ; 
% directory = strcat('/home/qbrissaud/Documents/Results/GJI_PAPER/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4_3D/');
% directory = strcat('/home/qbrissaud/Documents/Results/LAST_GRAVI_Roland/200PROCS/');
% directory = strcat('./');
% OFd = [rootd,'OUTPUT_FILES/'];
simuType = 1; % 1 for forcing and 2 for attenuation
% directory = strcat('/home/qbrissaud/Documents/FD/FD_14/');
% directory = strcat('/home/qbrissaud/Documents/Results/GJI_PAPER/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4_3D/');
% directory = strcat('/home/qbrissaud/Documents/Results/GJI_PAPER/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4/') ; 
% directory = strcat('/home/qbrissaud/Documents/Results/GJI_PAPER/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4/');

if(not(OFd(end)=='/')); OFd = [OFd,'/']; end;
% atmmodelfile=[directory,'1D_iso_enhanced.txt'];
atmmodelfile=[OFd,'input_atmospheric_model.dat'];
parfile=[OFd,'input_parfile'];
int_file=[OFd,'input_interfaces'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beginning of the program %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
% Stations parameters
% nstat      = 7 ;             % Number of station.
% dx_station = 60000; %58333.3333333;   % x-distance between stations
% z_station  = 100250.0;        % Height along z of stations
% % x_0        = 150000.0; % first station along x
% x_0        = 600250.0; % first station along x
% istattab=[6:7:48]
% xstattab=[450250.0:50000.0:750250.0]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time and space domain
% variables (parfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp([' ']);
disp(['[',mfilename,', INFO] Loading simulation parameters.']);
% Automatically load parameters from simulation's OUTPUT_FILES directory.
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
subsample_dt = T_0/30; % if subsample == 1, wanted dt for subsampling: ~30 pts/period or higher is safe
% Atmospheric model.
switch(readExampleFiles_extractParam(parfile, 'MODEL', 'string'))
  case 'default'
    externalDGAtmosModel = 0;
  case 'external_DG'
    externalDGAtmosModel = 1;
  otherwise
    error('kok');
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
mult_tSpan = 1; % multiplier for time span (used for analytical solution)
mult_xSpan = 1; % multiplier for x span (used for analytical solution)
mult_ySpan = 1; % multiplier for y span (used for analytical solution)
% fs   = 1/dt ; %not used anywhere
% fmax = 1; %not used ??
% fmin = -fmax; %not used ??

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atmospheric Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([' ']);
disp(['[',mfilename,', INFO] Preparing atmospheric model.']);
if(externalDGAtmosModel)
  % External DG atmospheric model, load it.
  disp(['[',mfilename,'] Loading atmospheric model from external file in ''',OFd,'''.']);
  [ALT, RHO, ~, SOUNDSPEED, ~, ~, GRA, NSQ, ~, ~, ~, ~, ~, WIND, ~, ~, ~] = extract_atmos_model(atmmodelfile, 3, 0, 0); % plot_model(atmmodelfile,'-','k',[]);
%   H        = (ALT(2)-ALT(1))/log(RHO(1)/RHO(2))
  % Set H.
  H        = (-ALT(2))/log(RHO(2)/RHO(1)); % NEEDED FOR NSQ AND ANALYTIC SOLUTION
  disp(['[',mfilename,'] Scale height computed,         = ',num2str(H),'.']);
%   rho_0    = RHO(1); % used nowhere
%   Nsqt      = -(GRA(1)*( -1/H - GRA(1)/(velocity(1)^2) )); % used nowhere

  % Set NSQ.
  error(['[',mfilename,', ERROR] N^2 computation for external models not implemented yet.']);
  
  % Set wind_x.
  if(max(abs(diff(WIND)))==0)
    wind_x = WIND(1);
    disp(['[',mfilename,'] Wind_x loaded,                 = ',num2str(wind_x),'. Constant with altitude.']);
  else
    wind_x = WIND(1);
    disp(['[',mfilename,'] Wind_x loaded,                 = ',num2str(wind_x),'. Not constant with altitude, setting wind_x as wind at altitude 0.']);
  end
  
else
  % Internal model, load it.
  disp(['[',mfilename,'] Building atmospheric model from parfile in ''',OFd,'''.']);
  
  % Set H.
  H = readExampleFiles_extractParam(parfile, 'SCALE_HEIGHT', 'float'); % NEEDED FOR NSQ AND ANALYTIC SOLUTION
  disp(['[',mfilename,'] Scale height loaded,           = ',num2str(H),'.']);
  
  USE_ISOTHERMAL_MODEL = readExampleFiles_extractParam(parfile, 'USE_ISOTHERMAL_MODEL', 'bool');
  if(USE_ISOTHERMAL_MODEL)
    % isothermal case
    disp(['[',mfilename,'] Isothermal case.']);
    error(['[',mfilename,', ERROR] Isothermal case not implemented yet.']);
    % need to produce SOUNDSPEED
    GRA = readExampleFiles_extractParam(parfile, 'gravity', 'float'); % ONLY NEEDED FOR NSQ
  else
    % isobaric case
    disp(['[',mfilename,'] Isobaric case.']);
    SOUNDSPEED = readExampleFiles_extractParam(parfile, 'sound_velocity', 'float');
    disp(['[',mfilename,'] Speed of sound loaded,         = ',num2str(SOUNDSPEED),'.']);
    GRA = 0; % ONLY NEEDED FOR NSQ
    disp(['[',mfilename,'] Gravity set for isobaric case, = ',num2str(GRA),'.']);
  end
  
  % Set NSQ.
  NSQ = -(GRA*( -1/H - GRA/(SOUNDSPEED^2) )); % NOT SURE ABOUT THAT FORMULA, CHECK
  disp(['[',mfilename,'] N^2 computed,                  = ',num2str(NSQ),'.']);
  
  % Set wind_x.
  wind_x = readExampleFiles_extractParam(parfile, 'wind', 'float');
  disp(['[',mfilename,'] Wind_x loaded,                 = ',num2str(wind_x),'.']);
end
disp(['[',mfilename,'] Setting wind_y = wind_z = 0.']);
wind_y = 0;
wind_z = 0;

% Prepare figure name.
savefigname = [OFd,'compare2FKanalytic'];
if(externalDGAtmosModel)
  savefigname = [savefigname, '_externalDG'];
else
  if(USE_ISOTHERMAL_MODEL)
    savefigname = [savefigname, '_isobaric'];
  else
    savefigname = [savefigname, '_isothermal'];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading Results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([' ']);
disp(['[',mfilename,', INFO] Loading synthetics.']);
% data = readAndSubsampleSynth(OFd, 2, 'BXZ', 'semv', 1, 1, 0); sig_t = data(:,1)'; sig_v = data(:,2)'; figure(); plot(sig_t, sig_v); % test
[~, signal_type, vname_long, vunit] = seismotypeNames(seismotype, readExampleFiles_extractParam(parfile, 'USE_DISCONTINUOUS_METHOD', 'bool'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical model variables, as in parfile.
% xo     = 750000.0;
% xo     = 250000.0;
% lambdo = 10000.0;
% lambdo = 80000.0;
% V      =  639.52000;       % must be in m/s
% A      = 1.0;
% rho_cst = 0.02 ;
% Time function parameters
% perio  = 1600.0;
% to     = 1400.0;
% perio  = 50.0;
% to     = 65.0;
% xo = 600000.0;
% wind in x direction
% wi=99.999
% %  wi=4.9
% wi = 10;
% wi = 0.0
% rho_cst = 0.02 ;
% A = 1;
% V = 639.52; % [m/s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading data from the atmospherical model
% atmos = load(atmmodelfile) ;
% atmos = load('1D_iso_enhanced.txt');
% % CONVERT OLD FILES TO NEW ONES
% atmos = load([rootd,'1D_iso_enhanced.txt']);
% alti_atm = atmos(:,1); rho = atmos(:,2); velocity = atmos(:,3); gravity = atmos(:,4); Nsq = atmos(:,5); Pressure = atmos(:,6);
% cpcst=1005; cvcst=717.3447537473;
% wind = 0;
% rewrite_atmos_model([rootd,'atmospheric_model.dat'], [], alti_atm, rho, rho*0, velocity, Pressure, rho*0, gravity, Nsq, rho*0, rho*0, rho*0, rho*0, rho*0, rho*0+wind, rho*0+cpcst, rho*0+cvcst, rho*0);
% alti_atm = atmos(:,1) ;
% rho      = atmos(:,2) ;
% velocity = atmos(:,3) ;
% gravity  = atmos(:,4) ;
% Nsq      = atmos(:,5) ;
% Pressure = atmos(:,6) ;



if(data_test0_or_readRun1 == 1)
  %%%%%%%%%%%%%%%%%%%%%%
  % Load stations' data.
%   pos_stat = load(strcat(directory,'STATIONS')) ; % first code version
  [x_stat, z_stat, y_stat, ~] = loadStations(OFd);
  pos_stat = [y_stat, z_stat, x_stat]; % try to do something to stick with rest of code
  
  % Build istattab.
  istattab = 1:numel(x_stat);
  nstat = length(istattab);
  
  % istattab = 1:size(pos_stat(:,1),1);
  % istattab = [1:18:81] ;
  % istattab = [37 38 39 40 41 42 43 44 45];
  % istattab = [36 37 44 45]+9;
  %coef = 9; alti = 4; istattab = [alti+coef alti+2*coef alti+7*coef alti+8*coef];
  % istattab = [1:5] + 1*9;
%   istattab = [2 4 6 8] + 2*9;
%   istattab = [1,2,3,4]; % test
  % istattab = 5;
  % istattab = [1 2 8 9] + 4*9
  % istattab = [1 10 19 28] + 8
  % istattab = [10 19 55 64] + 0
  % istattab = [19 28 46 55] + 0
  % istattab = [13 22 49 58]
  % istattab = [1 19 55 73]
  % istattab = [1 10 64 73]+1
  %istattab = 10:18:81;
  %%%%%%%%%%%
  % TEST GRAVI NO WIND
  % ystattab = pos_stat(istattab,1)+250;
  % zstattab = pos_stat(istattab,2)+250 ;
  % zstattab(1) = zstattab(1) + 20000;
  % zstattab(2) = zstattab(2) + 5000;
  % zstattab(3) = zstattab(3) + 20000;
  % zstattab(4) = zstattab(4) + 20000;
  % xstattab = pos_stat(istattab,3)+250 ;
  % xstattab(3:4) = xstattab(3:4) + 0;
  % xstattab(1) = xstattab(1) - 5000;
  % xstattab(2) = xstattab(2) - 5000;
  %%%%%%%%%%%
  % TEST GRAVI  WIND
  % ystattab = pos_stat(istattab,1);
  % zstattab = pos_stat(istattab,2) ;
  % xstattab = pos_stat(istattab,3) ;
  % zstattab(1) = zstattab(1) + 7200;
  % xstattab(1) = xstattab(1) - 0;
  % zstattab(2) = zstattab(2) + 5000;
  % xstattab(2) = xstattab(2) + 3000;
  % zstattab(3) = zstattab(3) + 11000;
  % xstattab(3) = xstattab(3) - 2500;
  % zstattab(4) = zstattab(4) + 13500;
  % xstattab(4) = xstattab(4) + 1500;
  %%%%%%%%%%%
  % TEST GRAVI  WIND
  ystattab = pos_stat(istattab, 1);
  zstattab = pos_stat(istattab, 2);
  xstattab = pos_stat(istattab, 3);
  % zstattab(1) = zstattab(1) + 10000;
%   xstattab = xstattab + 000;
  % xstattab(4) = xstattab(4) + 0;
%   ystattab = ystattab + 000;

  % ???
%   X_0 = (pos_stat(end,1) + pos_stat(1,1))/2;
%   X_0 = X_0 + 0;
%   yo = X_0;
  % yo = (pos_stat(end,3) + pos_stat(1,3))/2;
  
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
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % First loop over the stations to read synthetics' files.
  % Subsample the results to save memory space.
%   nsub=100;
  for locStatNumber = 1:nstat
      gloStatNumber = istattab(locStatNumber);
      if(gloStatNumber < 10)
         prefix = 'S000';
      elseif(gloStatNumber < 100)
         prefix = 'S00'; 
      else
         prefix = 'S0';
      end
      % reading of the horizontal component
%       file = strcat(directory,prefix,num2str(istat),'.AA.BXX.sem',signal_type) ;
%       file = strcat(OFd,'AA.',prefix,num2str(gloStatNumber),'.BXX.sem',signal_type) ;
%       data = load(file) ;
%       if (istat == 1)
%       nt = floor(max(size(data))/nsub) ;
%       end
%       Xtime(locStatNumber,1:nt) = data(1:nsub:max(size(data)),1)' ;
%       Xamp(locStatNumber,1:nt) = data(1:nsub:max(size(data)),2)' ;
      
      [tmpDat, nt] = readAndSubsampleSynth(OFd, gloStatNumber, 'BXX', ['sem', signal_type], subsample, subsample_dt, locStatNumber);
      Xtime(locStatNumber,1:nt) = tmpDat(:, 1)';
      Xamp(locStatNumber,1:nt) = tmpDat(:, 2)';
      % reading of the vertical component
%       file = strcat(directory,prefix,num2str(istat),'.AA.BXZ.sem',signal_type) ;
%       file = strcat(OFd,'AA.',prefix,num2str(gloStatNumber),'.BXZ.sem',signal_type) ;
%       data = load(file) ;
%       nt=floor(max(size(data))/nsub);
%       Ztime(locStatNumber,1:nt) = data(1:nsub:max(size(data)),1)' ;
%       Zamp(locStatNumber,1:nt) = data(1:nsub:max(size(data)),2)' ;
      [tmpDat, nt] = readAndSubsampleSynth(OFd, gloStatNumber, 'BXZ', ['sem', signal_type], subsample, subsample_dt, locStatNumber);
      Ztime(locStatNumber, 1:nt) = tmpDat(:, 1)';
      Zamp(locStatNumber, 1:nt) = tmpDat(:, 2)';
  end
  clear('tmpDat');
%   clear data
%   dt = dt;
  subsampledDt = Ztime(locStatNumber, 2)-Ztime(locStatNumber, 1);
  nSamp = nt;
%   dt = subsampledDt; % update dt
else % Else on data_test0_or_readRun1.
  % Test data.
  error('not implemented')
  % nSamp=10000;
  % xstattab = [5000];
  % zstattab = [10000];
  % xo = xmax/2;
  % nstat = 1;
  % Ztime(1,:) = linspace(0,nSamp,dt);
end % Endif on data_test0_or_readRun1.

disp(['[',mfilename,'] Synthetics loaded: ', num2str(nSamp), ' samples each.']);
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of the analytical solution.
disp([' ']);
disp(['[',mfilename,', INFO] Building analytical solution.']);
tmax_anal = t_0 + max(zstattab)/SOUNDSPEED(1) + 0.5*T_0 + 1*T_0; % go as far as 2 times the time necessary (roughly) for the signal to pass last station (t_0 + travel time + half period), + full period for safety
dt_anal   = T_0 / (2*100); % ensure fmax=1/dt >> 2*f0=2/T_0 by setting 1/dt = (2*2.5)/T_0
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
% [Brissaud et al., 2016, Section 5.1]: "Calculation of wavenumbers
% kx = 2π/λx and ky = 2π/λy for all spatial wavelengths λx,y."
% omega = 2.0*pi()*(1.0/(dt*NFFT2))*[[0:1:NFFT2/2] [0.0-[NFFT2/2-1:-1:1]]];
omega = 2.0*pi()*(1.0/(dt_anal*NFFT2))*[[0:1:NFFT2/2] [0.0-[NFFT2/2-1:-1:1]]];
kx =    2.0*pi()*(1.0/(dx     *NFFT1))*[[0:1:NFFT1/2] [0.0-[NFFT1/2-1:-1:1]]];
ky =    2.0*pi()*(1.0/(dy     *NFFT3))*[[0:1:NFFT3/2] [0.0-[NFFT3/2-1:-1:1]]];
[KX, KY, Omega] = meshgrid(kx, ky, omega);
% Assuming constant Nsq
Nsqtab = NSQ(1, 1) + 0.0*Omega;
% KX = Omega/velocity(1);
% see occhipinti 2008 for analytical solution (appendix A1)
onestab = 0.0*Nsqtab + 1.0;
% Omega_intrinsic = Omega-(wi)*KX;
Omega_intrinsic = Omega - wind_x*KX - wind_y*KY;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Brissaud et al., 2016, Section 5.1]: "Calculation of kz from dispersion
% relations for all wavenumbers kx, ky and time frequencies (see Appendix
% B)."
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

% Test new formula (WORKS QUITE GOOD ON ISOBARIC, WITH Mz = ifftn(filt.*TFMo);)
KZ = sqrt( (wind_x.^2 - SOUNDSPEED.^2).*KX.^2 - 2*KX.*Omega + Omega.^2) ./ SOUNDSPEED;

ind1=find(isnan(KZ));
ind2=find(isinf(KZ));
disp(['[',mfilename,'] Number of NaNs in KZ: ',num2str(numel(ind1)),'.']);
disp(['[',mfilename,'] Setting those NaNs to zeros.']);
KZ(ind1)=0.0;
disp(['[',mfilename,'] Number of Infs in KZ: ',num2str(numel(ind2)),'.']);
disp(['[',mfilename,'] Setting those Infs to zeros.']);
KZ(ind2)=0.0;
% imaginary part of KZ should be positive in order to attenuate the signal
indimag = find(imag(KZ)<0);
disp(['[',mfilename,'] Number of KZ such that Im(KZ)<0: ',num2str(numel(indimag)),'.']);
disp(['[',mfilename,'] Setting those to their conjugate (only flips the sign of the imaginary part).']);
KZ(indimag) = conj(KZ(indimag));
% real(KZ) should be positive for positive frequencies and negative for
% negative frequencies in order to shift signal in positive times
% restore the sign of KZ depending on Omega-wi*KX
%     KZnew=real(KZ).*sign((Omega-wind_x*KX)).*sign(KX)+1i*imag(KZ);
% !!! Why KZ should have a sign opposite to Omega for GW NOT UNDERSTOOD !!!
% => because vg perpendicular to Vphi ?
KZnew = 0.0 - real(KZ).*sign(Omega_intrinsic) + 1i*imag(KZ);
KZ=KZnew;

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


% if(seismotype == 1) % REMOVED FOR TESTS
%   if(simuType == 1)   % forcing% REMOVED FOR TESTS



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
      z_station = zstattab(gloStatNumber);
      filt = exp(1i*(KZ*z_station));
%       shifttime=0.0-real(KZ*z_station)./Omega;
%       filt = exp(1i*KZ*z_station).*(shifttime>=0);
%       
      % what
%       filt(1, 1, :) = 0.0; % why

      % remove waves assumed to propagate downward
%       indw=find((real(KZ)>0.0));
%       indw=find((real(KZ)<0.0));
%       filt(indw)=0.0;

%       Mz = exp(z_station/(2*H)) * ifftn(filt.*TFMo);
%       Mz = 0.5*exp(z_station/(2*H)) * ifftn(filt.*TFMo);
      Mz = ifftn(filt.*TFMo);
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
      
       % IF NOT PART OF THE MESH
      if (ix<=0)
        ix=1;
      end
      if (iy<=0)
        iy=1;
      end
%           if(istat < 3)
%               ix = ix-1;
%           else
%               ix = ix+1;
%           end
%           synf(istat,:)=real(Mz(iy,ix,:));
      synf(locStatNumber, :) = real(Mz(iy, ix, :));
%       synf(locStatNumber,:) = real(Mz(ix,:));
% remove offset a the beginning may be due to aliasing
%           synf(istat,:)=synf(istat,:)-synf(istat,1);
%           synf_X(istat,:)=synf_X(istat,:)-synf_X(istat,1);
    end

%   end% REMOVED FOR TESTS
% end% REMOVED FOR TESTS

% test to see what forms has the analytic solution.
disp([' ']);
disp(['[',mfilename,', INFO] Checking analytical solution.']);
f=figure();
% colours=jet(nstat);
% colours=prism(nstat);
colours=winter(nstat);
for i = 1:nstat;
  plot(t,synf(i,:),'displayname',num2str(zstattab(istattab(i))),'color',colours(i,:)); hold on;
end
legend();
prettyAxes(f);
isThisOk=-1;
while(not(ismember(isThisOk,[0,1])))
  isThisOk=input(['[',mfilename,'] Does this analytical solution look ok (0 for no, 1 for yes)? > ']);
end
if(not(isThisOk))
  error('stop kek');
end

% TEST ONLY
% synf = Zamp;
% 
% test only !!!!
% synf_X=synf;
%

% clear TFMx;
% clear Mx ;
% clear Mz ;

timef = t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure of the nstat vertical components and synthetic signals against 
% time
disp([' ']);
disp(['[',mfilename,'] Plotting everything.']);
f=figure('units','normalized','outerposition',[0 0 1 1]);
ax = [];
for locStatNumber = 1 : nstat
  ax(locStatNumber) = subplot(nstat,1,nstat-locStatNumber+1); 
  hold on
  if(seismotype == 1)
    ID_t_anal_end = find(timef<=Ztime(locStatNumber,end), 1, 'last');
%     plot(timef(1:timef_final),real(synf(locStatNumber,1:timef_final)),'Color',[0 0 0],'LineWidth',1)
    plot(timef(1:ID_t_anal_end),synf(locStatNumber, 1:ID_t_anal_end),'Color',[0 0 0],'LineWidth',1,'displayname','analytical');
%     plot(timef,real(synf(locStatNumber,:))/max(real(synf(locStatNumber,:))),'Color',[0 0 0],'LineWidth',1)
%     plot(Ztime(locStatNumber,:),real(synf(locStatNumber,1:length(Ztime(locStatNumber,:))))/max(real(synf(locStatNumber,1:length(Ztime(locStatNumber,:))))),'Color',[0 0 0],'LineWidth',1)
  else
%     plot(Ztime(locStatNumber,:),syn(locStatNumber,:),'Color',[0 0 0],'LineWidth',1)
    plot(t,synf(locStatNumber,:),'Color',[0 0 0],'LineWidth',1,'displayname','analytical')
  end
  
  % Plot synthetic.
%   Ztime_temp = Ztime(locStatNumber,1:floor(dt/subsampledDt):end);
%   Zamp_temp  = Zamp(locStatNumber,1:floor(dt/subsampledDt):end);
%   SIZE_R = size(Ztime_temp)
%   disp(['[',mfilename,'] Plotting a ',num2str(numel(Ztime_temp)),' long array.']);
%         size((Zamp_temp-real(synf(istat,1:length(Ztime_temp))))')
%   DIFF(locStatNumber,:) = (Zamp_temp-real(synf(locStatNumber,1:length(Ztime_temp))))';
  plot(Ztime(locStatNumber,1:nt),Zamp(locStatNumber,1:nt),'-.k','LineWidth',2,'displayname','synthetic');
  
  % Produce and plot difference.
  if(min(diff(Ztime(1,:))) < min(diff(t)))
    % if dt_sythetic < dt_analytic, interpolate analytic on synthetic t (Ztime)
    err_t = Ztime(locStatNumber,1:nt);
    err_synth = Zamp(locStatNumber,1:nt);
    err_anal = interp1(t, synf(locStatNumber,:), err_t);
  else
    % if dt_sythetic < dt_analytic, interpolate synthetic on analytic t (t)
    err_t = t;
    err_synth = interp1(Ztime(locStatNumber,1:nt), Zamp(locStatNumber,1:nt), err_t);
    err_anal = synf(locStatNumber,:);
  end
  err_v = factor_err * abs(err_synth - err_anal);
  plot(err_t,err_v,':k','LineWidth',1.5,'displayname',['$',num2str(factor_err),'{\times}|$anal.$-$synth.$|$']);
  
  % x-y labels.
  if (locStatNumber == round(nstat/2))
%     ylabel([variable,' along z-axis (m)']);
    ylabel(['$',vname_long,'_z$ ',vunit]);
  end
  if (locStatNumber == 1)
    xlabel('time [s]');
  else
    xticklabels([]);
  end
%   text(00.943*max(Ztime(locStatNumber,:)),0.95,['X position : $',num2str(xstattab(locStatNumber)),'km$ (along x)']) 
%   legend('Analytical','Modeled','10*(Mo-An)','Location','East')
%   legend('Location','East');
%   if(locStatNumber==1)
  if(locStatNumber==nstat)
    % legend only on first plot
    legend('location','northeast');
  end
end
linkaxes(ax, 'x');
xlim([max(min(min(Ztime)),min(t)), min(max(max(Ztime)),max(t))]);
% title(['Gravito-acoustic wave propagation. Stations at altitude : ', num2str(zstattab(locStatNumber)),'km (along z)']);
title(['Stations at altitudes ', sprintf('%.0f ',zstattab(istattab)),' [m]']);
addpath('/usr/local/matlab/r2018a/toolbox/tightfig'); tightfig; % eventually tighten fig
prettyAxes(f);
customSaveFig(savefigname,{'fig','eps','jpg'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure of the nstat horizontal components and synthetic signals against 
% time
% figure
% for istat = 1 : nstat
%     
%     ax(istat)=subplot(nstat,1,nstat-istat+1) ; 
%     hold on
%     
%     if(seismotype ==1)
%          plot(timef,real(synf_X(istat,:)),'Color',[0 0 0],'LineWidth',1)
% %        plot(timef,real(synf_X(istat,:))/max(real(synf_X(istat,:))),'Color',[0 0 0],'LineWidth',1)
% 
% %                plot(Ztime(istat,:),real(synf_X(istat,1:length(Xtime(istat,:))))/max(real(synf_X(istat,1:length(Xtime(istat,:))))),'Color',[0 0 0],'LineWidth',1)
%     else
%     plot(Xtime(istat,:),syn(istat,:),'Color',[0 0 0],'LineWidth',1)
%     end
%   
%      plot(Xtime(istat,1:nt),Xamp(istat,1:nt),'-.k','LineWidth',2)
% %    plot(Xtime(istat,:),Xamp(istat,:)/max(Xamp(istat,:)),'-.k','LineWidth',2)
%     
%     if (istat == round(nstat/2))
%         ylabel([variable,' x-axis (m)'] ,'FontSize',14.3)
%     end
%     if (istat == 1)
%         xlabel('time (s)','FontSize',14.3)
%     end
%     
%     text(00.943*max(Xtime(istat,:)),0.95,['X position : $',num2str((istat-1)*dx_station + x_0),'km$ (along x)']) 
%     
%     legend('Analytical','Modeled','10*(Mo-An)','Location','West')
% 
% end
% 
% linkaxes(ax,'x')
% % xlim([0 20000])
% 
% title(strcat('Gravito-acoustic wave propagation. Stations at altitude : ', num2str(z_station/1000),'km (along z)'),'FontSize',24)



