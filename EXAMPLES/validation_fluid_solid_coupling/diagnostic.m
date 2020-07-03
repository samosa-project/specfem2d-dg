% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%  TODO.

clear all;
% close all;
clc;
thisFilePath = mfilename('fullpath'); spl = split(thisFilePath,filesep); spl(end)=[]; spl = join(spl,filesep); thisFilePath = spl{1}; cd(thisFilePath);
addpath(genpath('../../utils_new'));
% set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
% set(0, 'defaultTextFontSize', 24); set(0, 'defaultAxesFontSize', 24);
% set(0, 'DefaultTextInterpreter', 'latex');
% set(0, 'DefaultLegendInterpreter', 'latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFD='OUTPUT_FILES';
% OFD='OUTPUT_FILES_LNS_S2F_plotvz';
OFD='OUTPUT_FILES_LNS_F2S_plotvz';

inline1_table0 = 1;
plot_timeseries = 0;
normalise_ylims = 0;

classical1_zhang0 = 1; % theoretical values choice

mindelayincrefl = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stations' geometry.
th_separation = 20; % we expect stations to be separated by 20 [m] across the interface
breakOnWrongSeparation = 0; % break script if stations are not rightly spaced?

% Mesh geometry.
slope_angle = 30 *pi/180;
triangle_xmin = -75; % stations to the left of the triangle should monitor orthogonal conversion, and to the right oblique conversion

% Deduce wave parameters.
incident_angle_from_slope = slope_angle;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['**************** ', OFD, ' ****************']);

% Filename for report.
reportname = [OFD,'_report'];
if(exist(reportname,'file'))
  delete(reportname);
  disp(['[',mfilename,', INFO] Deleted existing previous report.']);
end
diary(reportname); % Start logging to file.

% Get paths.
OFD=[OFD,filesep];
parfile = [OFD,'input_parfile'];
sourcefile = [OFD,'input_source'];

% Test method used.
if(not(readExampleFiles_extractParam(parfile,'USE_DISCONTINUOUS_METHOD','bool')))
  error(['[',mfilename,'] Must use DG, set ''USE_DISCONTINUOUS_METHOD = .true.''']);
end
% Test models and get parameters.
default_model = strcmp(readExampleFiles_extractParam(parfile,'MODEL','string'),'default');
isobaric_model = not(readExampleFiles_extractParam(parfile,'USE_ISOTHERMAL_MODEL','bool'));
if(not(default_model & isobaric_model))
  error('[',mfilename,'] DG model must be default isobaric, set ''MODEL = default'' and ''USE_ISOTHERMAL_MODEL = .false.''');
end
[rho_1, vp_1, Z_1, rho_2, vp_2, vs_2, Z_2P, Z_2S] = grab_models(parfile); % Get impedances from the models.
[s2f1_or_f2s0, fig_title, f0] = check_test_case(parfile, sourcefile); % Load interface (should be z=0), source altitude, frequency of the source. Return deduced test case and corresponding title.
spatial_wavelength_f = vp_1/f0; % maybe will be useful
spatial_wavelength_s = vp_2/f0; % maybe will be useful

disp([' ']);
disp(['> ',fig_title, ' coupling, Z_s = ',num2str(Z_2P), ', Z_f = ',num2str(Z_1),'.']);

% Load stations and synthetics.
STATIONS = importdata([OFD,'STATIONS']); STATPOS = STATIONS.data(:,1:2);
nstat = size(STATIONS.data,1); station_ids = 1:nstat;
% Match station couples.
% disp(['[] (ID, x_s, z_s):']); disp([(1:size(STATPOS,1))',STATPOS]); pause
disp(['> Association of stations by couples.']);
couples = [[(1:3)',(4:6)']; [(7:9)',(10:12)']];
% check distance between each station couple
distance_between_each_couple = sum((STATPOS(couples(:,1),:) - STATPOS(couples(:,2),:)).^2,2).^0.5;
if(breakOnWrongSeparation)
  if(any(abs(distance_between_each_couple-th_separation)>1e-9))
    error(['[',mfilename,',ERROR] Some station couple is separated by more than ',num2str(th_separation),' [m].']);
  end
end
[time, amp] = load_synthetics(OFD, parfile, station_ids); % Load synthetics.

% If plot asked, plot time series.
if(plot_timeseries)
  plotOneByOne(time, squeeze(amp(1,:,:)), station_ids, STATPOS(:,1), STATPOS(:,2), normalise_ylims, fig_title, '$\delta p$ [Pa] or $v_x$ [m/s]');
end

format_positions = ['%',num2str(floor(log10(max(abs(STATPOS(:,2)))))+6)','.1f'];
format_values = ['%11.4e'];
time_format = ['%.4f'];

% Loop over couples of stations and compute ratios.
disp([' ']);
disp(['> Computing ratios.']);

for i=1:size(couples,1)
% for i=1
%   ztotest=zlist(i);
%   correspondingIDs = find(abs(STATPOS(:,2))==ztotest);
%   % with sorting above, first of correspondingIDs is incoming and last is transmitted
%   incoming_reflected = Zamp(correspondingIDs(1),:);
%   incoming_reflected_t = Ztime(correspondingIDs(1),:);
%   transmitted = Zamp(correspondingIDs(2),:);
%   transmitted_t = Ztime(correspondingIDs(1),:);
  correspondingIDs = couples(i, :);
  
  z_of_this_couple = STATPOS(correspondingIDs, 2);
  if(s2f1_or_f2s0)
    % if s2f, incoming is lower z
    ID_inc = find(z_of_this_couple==min(z_of_this_couple));
    ID_out = find(z_of_this_couple==max(z_of_this_couple));
  else
    % if f2s, incoming is upper z
    ID_inc = find(z_of_this_couple==max(z_of_this_couple));
    ID_out = find(z_of_this_couple==min(z_of_this_couple));
  end
  
  % Distinguish between orthogonal and oblique incidence
  if(all(STATPOS(correspondingIDs, 1)<triangle_xmin))
    % orthogonal
    incident_angle = 0;
  elseif(all(STATPOS(correspondingIDs, 1)>triangle_xmin))
    % oblique
    incident_angle = incident_angle_from_slope;
  else
    error(['at exactly ',num2str(triangle_xmin), ' or if on either side of ',num2str(triangle_xmin), ', cannot compute']);
  end
  
  % Compute experimental ratios.
  if(s2f1_or_f2s0)
    error('not impl');
  else
    % FLUID-TO-SOLID
    incoming_reflected_v = squeeze(amp(1, correspondingIDs(ID_inc),:)); % pressure waves
    incoming_reflected_t = time(correspondingIDs(ID_inc),:);
    transmitted_vx = squeeze(amp(1,correspondingIDs(ID_out),:)); % solid waves
    transmitted_vz = squeeze(amp(2,correspondingIDs(ID_out),:)); % solid waves
    transmitted_t = time(correspondingIDs(ID_out),:);
    
    transmitted_v = (transmitted_vx.^2+transmitted_vz.^2).^0.5; % shit version
    
    i2 = snells(vp_1, vp_2, incident_angle); % deduce i2
    j2 = snells(vp_1, vs_2, incident_angle); % deduce j2
    transmitted_vp = sin(i2)*transmitted_vx + cos(i2)*transmitted_vz; % ground velocity along P-wave direction
    transmitted_vs = cos(j2)*transmitted_vx + sin(j2)*transmitted_vz; % ground velocity 90° from S-wave direction
  end
  
  % Obtain theoretical values.
  if(s2f1_or_f2s0)
    if(classical1_zhang0)
      [R_s2f, T_s2f] = ReflexionTransmissionCoefs(Z_2P, Z_1, incident_angle, snells(vp_2, vp_1, incident_angle));
    else
      [R_s2f, T_s2f] = ReflexionTransmissionCoefsZhang(s2f1_or_f2s0, vp_1, Z_1, vp_2, vs_2, Z_2P, Z_2S, incident_angle);
      T_s2f = sum(T_s2f); % assume we are measuring both the P and the S
    end
    R_s2f = abs(R_s2f); % work in magnitude only
    % Rationale:
    % pt = T*pi, pr = R*pi
    % pt = Zs*vi => pt = T*Zs*vi => vi/pt = 1/(T*Zs)
    % (vi+vr) = (1+R)*vi => (vi+vr)/pt = T*Zs/(1+R)
%     Vi_over_Pt_th = 1/(T_s2f*Z_1); % v incident over p transmitted
    Pt_over_Vi_th = T_s2f*Z_1; % p transmitted over v incident
  %   ViVr_over_Pt = (1+R_s2f)/(T_s2f*Z_s); % If too close to interface, v transmitted over (p incident + p reflected)
%     disp(['CONVERSION COEFFICIENT NOT SURE']);
%     disp([fig_title, ' coupling, Z_s = ',num2str(Z_s), ', Z_f = ',num2str(Z_f),', T = ',num2str(T_s2f), ', R = ',num2str(R_s2f), '.']);
%     spatial_wavelength = spatial_wavelength_s;
%     localT_th = Vi_over_Pt_th; localT_th_name = 'Vi/Pt';
    localT_th = Pt_over_Vi_th; localT_th_name = 'Pt/Vi';
    localR_th = R_s2f; localR_th_name = 'Vr/Vi';
  else
    % FLUID-TO-SOLID
    if(classical1_zhang0)
      [R_f2s, T_f2s] = ReflexionTransmissionCoefs(Z_1, Z_2P, incident_angle, snells(vp_1, vp_2, incident_angle));
    else
      [R_f2s, T_f2s] = ReflexionTransmissionCoefsZhang(s2f1_or_f2s0, vp_1, Z_1, vp_2, vs_2, Z_2P, Z_2S, incident_angle);
      T_f2s = sum(T_f2s); % assume we are measuring both the P and the S
    end
    
    % Rationale:
    % pt = Tpi, pr=Rpi
    % pt=vtZs => vtZs=Tpi => vt/pi=T/Zs
    % (pi+pr) = (1+R)pi => vt/(pi+pr)=T/(Zs*(1+R))
    Vt_over_Pi_th = T_f2s/Z_1; % v transmitted over p incident
  %   Vt_over_PiPr_th = T_f2s/(Z_s*(1+R_f2s)); % if too close to interface, v transmitted over (p incident + p reflected)
%     disp([fig_title, ' coupling, Z_s = ',num2str(Z_s), ', Z_f = ',num2str(Z_f),', T = ',num2str(T_f2s), ', R = ',num2str(R_f2s), '.']);
%     spatial_wavelength = spatial_wavelength_f;
    localT_th = Vt_over_Pi_th; localT_th_name = 'Vt/Pi';
    localR_th = R_f2s; localR_th_name = 'Pr/Pi';
  end
  
  maxit = 50;
  % Find incoming wave.
  inc_peak_time = findFirstPeak(incoming_reflected_t, incoming_reflected_v);
  inc_peak = incoming_reflected_v(incoming_reflected_t==inc_peak_time);
  % Find transmitted wave.
  transm_peak_time = findFirstPeak(transmitted_t, transmitted_v);
  transm_peak = transmitted_v(transmitted_t==transm_peak_time);
  % Find reflected wave.
  refl_peak_time = findFirstPeak(fliplr(incoming_reflected_t),fliplr(incoming_reflected_v), localR_th*0.6, 2, maxit); % Method below assumes there is exactly two peaks in the signal (tries to find 'first' peak of flipped array).
  refl_peak_time = findFirstPeak(incoming_reflected_t(incoming_reflected_t>inc_peak_time+mindelayincrefl),incoming_reflected_v(incoming_reflected_t>inc_peak_time+mindelayincrefl));
  if(refl_peak_time<inc_peak_time)
    error(['good job, you found the reflection before the incident wave']);
  end
  refl_peak = incoming_reflected_v(incoming_reflected_t==refl_peak_time);
  
  % Obtain experimental values.
  if(s2f1_or_f2s0)
    P_over_V_ex = transm_peak/inc_peak; % transmitted is pressure, incoming is velocity
  else
    % FLUID-TO-SOLID
    V_over_P_ex = transm_peak/inc_peak; % transmitted is velocity, incoming is pressure
  end
  RoverI = refl_peak/inc_peak;

%   threshold_1_peak = 0.1; % (* spatial_wavelength)
%   threshold_2_peaks = 0.5; % (* spatial_wavelength)

%   % Boundary cases.
%   if(ztotest<threshold_1_peak*spatial_wavelength)
%     % close enough to boundary, consider double peak (incoming + reflected)
%     if(f2s0_or_s2f1)
%       % s2f
%       localV_over_P_th = ViVr_over_Pt;
%       voverpname = '(Vi+Vr)/Pt';
%     else
%       % f2s
%       localV_over_P_th = Vt_over_PiPr_th;
%       voverpname = 'Vt/(Pi+Pr)';
%     end
%   elseif(ztotest<threshold_2_peaks*spatial_wavelength)
%     % far from boundary, but not far enough to separate the two peaks
%     disp(['@z=+/-',sprintf(format_positions,ztotest),': Skipping this station, because we can neither separate 2 peaks nor consider they form only 1. Incoming spatial wavelength = ',num2str(spatial_wavelength),' m.']);
%     continue;
%   else
%     % okay, change nothing (use first detected peak)
%   end

  % Display.
  if(inline1_table0)
%     disp(['  > Couple at (x_1, x_2, z_1, z_2)=(',sprintf(format_positions,STATPOS(correspondingIDs,:)),'):']);
    disp(['  > Incident angle = ',sprintf('%3.0f',incident_angle*180/pi),'°. (Incident wave found at t=',sprintf(time_format,inc_peak_time),'s.)']);
    disp(['  > (t=',sprintf(time_format,transm_peak_time),'s) T',fig_title,'_th = ',sprintf(format_values,localT_th),' m/s/Pa, T',fig_title,' = ',localT_th_name,' = ',sprintf(format_values,V_over_P_ex),' m/s/Pa, relative error = ',sprintf('%6.3f',100*abs(V_over_P_ex-localT_th)/abs(localT_th)),'%.']);
    disp(['  > (t=',sprintf(time_format,refl_peak_time),'s) R',fig_title,'_th = ',sprintf(format_values,localR_th),'       , R',fig_title,' = ',localR_th_name,' = ',sprintf(format_values,RoverI),'       , relative error = ',sprintf('%6.3f',100*abs(RoverI-localR_th)/abs(localR_th)),'%.']);
  end
  
  disp(' ');
  % TODO: save agreements.
end

% TODO: plot of saved agreements.

diary off; % End logging.

% Function for Snell's law

% Function for reflexion and transmission coefficients computation.
