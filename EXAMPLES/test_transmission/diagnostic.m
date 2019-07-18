% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%  TODO.

clear all;
% close all;
% clc;
addpath(genpath('../../utils_new'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFD='OUTPUT_FILES_FNS_F2S_1.0dx_1.0dt';
% OFD='OUTPUT_FILES_FNS_F2S_0.5dx_0.5dt';
% OFD='OUTPUT_FILES_LNS_F2S_1.0dx_1.0dt';
% OFD='OUTPUT_FILES_LNS_F2S_1.0dx_0.5dt';
% OFD='OUTPUT_FILES_LNS_F2S_0.5dx_0.5dt';
OFD='OUTPUT_FILES';

% Plots?
plot_timeseries = 0;
normalise_ylims = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['**************** ', OFD, ' ****************']);

% Filename for report.
reportname = [OFD,'_report'];
if(exist(reportname,'file'))
  delete(reportname);
  disp(['Deleted existing previous report.']);
end
diary(reportname); % Start logging to file.

% Get paths.
OFD=[OFD,filesep];
parfile = [OFD,'input_parfile'];
sourcefile = [OFD,'input_source'];

% Test method used.
if(not(readExampleFiles_extractParam(parfile,'USE_DISCONTINUOUS_METHOD','bool')))
  error(['must use DG, set ''USE_DISCONTINUOUS_METHOD = .true.''']);
end

% Test models used.
default_model = strcmp(readExampleFiles_extractParam(parfile,'MODEL','string'),'default');
isobaric_model = not(readExampleFiles_extractParam(parfile,'USE_ISOTHERMAL_MODEL','bool'));
if(not(default_model & isobaric_model))
  error('DG model must be default isobaric, set ''MODEL = default'' and ''USE_ISOTHERMAL_MODEL = .false.''');
end

% Get impedances from the models.
rho_f = readExampleFiles_extractParam(parfile,'surface_density','float');
vp_f = readExampleFiles_extractParam(parfile,'sound_velocity','float');
Z_f = rho_f*vp_f;
[~,mas]=readExampleFiles_extractParfileModels(parfile);
if(not(numel(mas)==2))
  error(['Must use exactly 2 models.']);
end
mas(readExampleFiles_extractParam(parfile,'id_region_DG','int'))=[]; % remove dg model
rho_s = mas.rho;
vp_s = mas.vp;
Z_s = rho_s*vp_s;

% Load interface (should be z=0), source altitude, frequency of the source.
% Compute wavelengths.
zinterface = readExampleFiles_extractParam(parfile,'coord_interface','float');
zs = readExampleFiles_extractParam(sourcefile,'zs','float');
if(numel(zs)>1)
  % If many sources found (in the case of S2F coupling), check if all same altitude and keep only one.
  if(numel(unique(zs))==1)
    zs = zs(1);
  else
    error(['Sources must be at same altitude.']);
  end
end
f0 = readExampleFiles_extractParam(sourcefile,'f0','float');
if(numel(f0)>1)
  % If many f0 found (in the case of S2F coupling), check if all same and keep only one.
  if(numel(unique(f0))==1)
    f0 = f0(1);
  else
    error(['Sources must have same frequency.']);
  end
end
spatial_wavelength_f = vp_f/f0;
spatial_wavelength_s = vp_s/f0; % maybe will be useful

% Load stations.
STATIONS = importdata([OFD,'STATIONS']);
STATPOS = STATIONS.data(:,1:2);
nstat = size(STATIONS.data,1);
station_ids = 1:nstat;

% Distinguish type of coupling based on source altitude.
if(zs>zinterface)
  % source is above ground, coupling is thus fluid to solid
  f2s0_or_s2f1 = 0;
else
  if(zs==zinterface)
    error(['source must not be on interface']);
  end
  % source is above ground, coupling is thus solid to fluid
  f2s0_or_s2f1 = 1;
end

% Set title, RT coeffs, and theoretical values for conversion.
if(f2s0_or_s2f1)
  fig_title = ['S2F'];
  [R, T] = ReflexionTransmissionCoefs(Z_s, Z_f);
  R = abs(R); % work in magnitude only
  % Rationale:
  % pt = T*pi, pr = R*pi
  % pt = Zs*vi => pt = T*Zs*vi => vi/pt = 1/(T*Zs)
  % (vi+vr) = (1+R)*vi => (vi+vr)/pt = T*Zs/(1+R)
  Vi_over_Pt_th = 1/(T*Z_s); % v incident over p transmitted
  ViVr_over_Pt = (1+R)/(T*Z_s); % If too close to interface, v transmitted over (p incident + p reflected)
  disp(['CONVERSION COEFFICIENT NOT SURE']);
else
  fig_title = ['F2S'];
  [R, T] = ReflexionTransmissionCoefs(Z_f, Z_s);
  % Rationale:
  % pt = Tpi, pr=Rpi
  % pt=vtZs => vtZs=Tpi => vt/pi=T/Zs
  % (pi+pr) = (1+R)pi => vt/(pi+pr)=T/(Zs*(1+R))
  Vt_over_Pi_th = -T/Z_s; % v transmitted over p incident
  Vt_over_PiPr_th = -T/(Z_s*(1+R)); % if too close to interface, v transmitted over (p incident + p reflected)
end
disp([fig_title, ' coupling, Z_s = ',num2str(Z_s), ', Z_f = ',num2str(Z_f),', T = ',num2str(T), ', R = ',num2str(R), '.']);

% Sort everything by increasing/decreasing z.
[~, isorted]=sort(STATPOS(:,2)); % s2f: increasing z
if(f2s0_or_s2f1==0)
  isorted = flip(isorted); % f2s: decreasing z
end
station_ids = station_ids(isorted);
STATPOS = STATPOS(isorted, :);

% Test if positions of stations are ok.
% They should be linked couples of stations.
if(not(all((STATPOS(1:nstat/2,:)+flipud(STATPOS(nstat/2+1:end,:)))==0,'all')))
  error('stations must be placed as couples of stations at the same z-distance from the interface, and at the same x position');
end
% Build list of altitudes to be checked.
zlist = abs(STATPOS(1:nstat/2, 2));

% Prepare synthetics loading.
subsample=0;wanted_dt = -1; % do not subsample synthetics
channel = 'BXZ';
type_display = readExampleFiles_extractParam(parfile,'seismotype','int');
[extension, ylabel_unknown] = getUnknowns(type_display, channel);
% Load synthetics.
for istat = 1:nstat
  istatglob=station_ids(istat);
  [data, nsamples] = readAndSubsampleSynth(OFD, istatglob, channel, extension, subsample, wanted_dt, istatglob);
  Ztime(istat, 1:nsamples) = data(:, 1)';
  Zamp(istat, 1:nsamples) = data(:, 2)';
end

% If plot asked, plot.
if(plot_timeseries)
  plotOneByOne(Ztime, Zamp, station_ids, STATPOS(:,1), STATPOS(:,2), normalise_ylims, fig_title, ylabel_unknown);
end

format_positions = ['%',num2str(floor(log10(max(abs(STATPOS(:,2)))))+3)','.1f'];
format_values = ['%11.4e'];

% Loop over couples of stations and compute ratios.
for i=1:numel(zlist)
% for i=1
  ztotest=zlist(i);
  correspondingIDs = find(abs(STATPOS(:,2))==ztotest);
  % with sorting above, first of correspondingIDs is incoming and last is transmitted
  incoming_reflected = Zamp(correspondingIDs(1),:);
  incoming_reflected_t = Ztime(correspondingIDs(1),:);
  transmitted = Zamp(correspondingIDs(2),:);
  transmitted_t = Ztime(correspondingIDs(1),:);
  
  incoming_reflected = abs(incoming_reflected); % work in magnitude
  transmitted = abs(transmitted); % work in magnitude
  
  inc_peak = incoming_reflected(incoming_reflected_t==findFirstPeak(incoming_reflected_t,incoming_reflected));
  transm_peak = transmitted(transmitted_t==findFirstPeak(transmitted_t,transmitted));
  
  if(f2s0_or_s2f1)
    % s2f
    spatial_wavelength = spatial_wavelength_s;
    voverpname = 'Vi/Pt     ';
    roveriname = 'Vr/Vi     ';
    V_over_P_ex = inc_peak/transm_peak; % transmitted is pressure, incoming is velocity
    localV_over_P_th = Vi_over_Pt_th;
  else
    % f2s
    spatial_wavelength = spatial_wavelength_f;
    voverpname = 'Vt/Pi     ';
    roveriname = 'Pr/Pi     ';
    V_over_P_ex = transm_peak/inc_peak; % transmitted is velocity, incoming is pressure
    localV_over_P_th = Vt_over_Pi_th;
  end

  threshold_1_peak = 0.1; % (* spatial_wavelength)
  threshold_2_peaks = 0.5; % (* spatial_wavelength)

  % Boundary cases.
  if(ztotest<threshold_1_peak*spatial_wavelength)
    % close enough to boundary, consider double peak (incoming + reflected)
    if(f2s0_or_s2f1)
      % s2f
      localV_over_P_th = ViVr_over_Pt;
      voverpname = '(Vi+Vr)/Pt';
    else
      % f2s
      localV_over_P_th = Vt_over_PiPr_th;
      voverpname = 'Vt/(Pi+Pr)';
    end
  elseif(ztotest<threshold_2_peaks*spatial_wavelength)
    % far from boundary, but not far enough to separate the two peaks
    disp(['@z=+/-',sprintf(format_positions,ztotest),': Skipping this station, because we can neither separate 2 peaks nor consider they form only 1. Incoming spatial wavelength = ',num2str(spatial_wavelength),' m.']);
    continue;
  else
    % okay, change nothing (use first detected peak)
  end

  % Display.
  disp(['@z=+/-',sprintf(format_positions,ztotest),', TRANSMISSION: ',voverpname,' = ',sprintf(format_values,V_over_P_ex),' m/s/Pa, theoretical value = ',sprintf(format_values,localV_over_P_th),' m/s/Pa, relative error = ',sprintf('%6.3f',100*abs(V_over_P_ex-localV_over_P_th)/abs(localV_over_P_th)),'%.'])

  % Reflection, only if >0.5*spatial_wavelength
  if(ztotest>threshold_2_peaks*spatial_wavelength)
    % Find peak. Method below assumes there is exactly two peaks in the signal (tries to find 'first' peak of flipped array).
    refl_peak = incoming_reflected(incoming_reflected_t==findFirstPeak(fliplr(incoming_reflected_t),fliplr(incoming_reflected)));

    RoverI = refl_peak/inc_peak;

    % Display.
    disp(['@z=+/-',sprintf(format_positions,ztotest),', REFLEXION:    ',roveriname,' = ',sprintf(format_values,RoverI),'       , theoretical value = ',sprintf(format_values,R),'       , relative error = ',sprintf('%6.3f',100*abs(RoverI-R)/abs(R)),'%.'])
  end

  % TODO: save agreements.
end

% TODO: plot of saved agreements.

diary off; % End logging.

% Function for reflexion and transmission coefficients computation.
function [R, T] = ReflexionTransmissionCoefs(Z_incoming, Z_outgoing)
  R = (Z_outgoing - Z_incoming) / (Z_incoming + Z_outgoing);
  T = 2*Z_outgoing / (Z_incoming + Z_outgoing);
end