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
OFD='OUTPUT_FILES_LNS_F2S_1.0dx_1.0dt';
% OFD='OUTPUT_FILES_LNS_F2S_1.0dx_0.5dt';
% OFD='OUTPUT_FILES_LNS_F2S_0.5dx_0.5dt';

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
f0 = readExampleFiles_extractParam(sourcefile,'f0','float');
spatial_wavelength_f = vp_f/f0;
spatial_wavelength_s = vp_f/f0; % maybe will be useful

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
  f2s0_or_s2f1 = 0;
end

% Set title, RT coeffs, and theoretical values for conversion.
if(f2s0_or_s2f1)
  fig_title = ['S2F'];
  [R, T] = ReflexionTransmissionCoefs(Z_s, Z_f);
  Vi_over_Pt_th = -T*Z_s; % v incident over p transmitted
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
  correspondingStations = find(abs(STATPOS(:,2))==ztotest);
  if(f2s0_or_s2f1)
    error(['not implementd']);
  else
    % with sorting above, first of correspondingStations is incoming and last is transmitted
    incoming_reflected = Zamp(correspondingStations(1),:);
    incoming_reflected_t = Ztime(correspondingStations(1),:);
    transmitted = Zamp(correspondingStations(2),:);
    transmitted_t = Ztime(correspondingStations(1),:);
    
%     t_inc_peak = findFirstPeak(incoming_reflected_t,incoming_reflected);
%     findFirstPeak(incoming_reflected_t,incoming_reflected,0.5,2,10,1,1); % verbose
    inc_peak = incoming_reflected(incoming_reflected_t==findFirstPeak(incoming_reflected_t,incoming_reflected));
    
%     findFirstPeak(incoming_reflected_t,incoming_reflected,0.5,2,10,1,1); % verbose
    transm_peak = transmitted(transmitted_t==findFirstPeak(transmitted_t,abs(transmitted)));
    
    V_over_P_ex = transm_peak/inc_peak;
    localV_over_P_th = Vt_over_Pi_th;
    
    if(f2s0_or_s2f1)
      voverpname = 'Vi/Pt     ';
      roveriname = 'Vr/Vi     ';
      spatial_wavelength = spatial_wavelength_s;
    else
      voverpname = 'Vt/Pi     ';
      roveriname = 'Pr/Pi     ';
      spatial_wavelength = spatial_wavelength_f;
    end
    
    % Boundary cases.
    if(ztotest<0.5*spatial_wavelength)
      % If distance is under the wavelength, 'incident' peak felt is in fact some combination of incident+reflected
      if(ztotest>0.15*spatial_wavelength)
        % if somewhat too far from interface, skip it because the combination is composite
        disp(['@z=+/-',sprintf(format_positions,ztotest),': Skipping this station, because neither can separate the 2 peaks, nor can consider they form 1.']);
        continue;
      else
        if(f2s0_or_s2f1)
          % s2f
          error(['not implemented'])
          voverpname = '(Vi+Vr)/Pt';
        else
          % f2s
          localV_over_P_th = Vt_over_PiPr_th;
          voverpname = 'Vt/(Pi+Pr)';
        end
      end
    end
    
    % Display.
    disp(['@z=+/-',sprintf(format_positions,ztotest),', TRANSMISSION: ',voverpname,' = ',sprintf(format_values,V_over_P_ex),' m/s/Pa, theoretical value = ',sprintf(format_values,localV_over_P_th),' m/s/Pa, relative error = ',sprintf('%6.3f',100*abs(V_over_P_ex-localV_over_P_th)/abs(localV_over_P_th)),'%.'])
    
    % Reflection, only if >0.5*spatial_wavelength
    if(ztotest>0.5*spatial_wavelength)
      % Find peak. Method below assumes there is exactly two peaks in the signal (tries to find 'first' peak of flipped array).
      refl_peak = incoming_reflected(incoming_reflected_t==findFirstPeak(fliplr(incoming_reflected_t),fliplr(incoming_reflected)));
      
      RoverI = refl_peak/inc_peak;
      
      % Display.
      disp(['@z=+/-',sprintf(format_positions,ztotest),', REFLEXION:    ',roveriname,' = ',sprintf(format_values,RoverI),'       , theoretical value = ',sprintf(format_values,R),'       , relative error = ',sprintf('%6.3f',100*abs(RoverI-R)/abs(R)),'%.'])
    end
    
    % TODO: save agreements.
  end
end

% TODO: plot of saved agreements.

diary off; % End logging.

% Function for reflexion and transmission coefficients computation.
function [R, T] = ReflexionTransmissionCoefs(Z_incoming, Z_outgoing)
  R = (Z_outgoing - Z_incoming) / (Z_incoming + Z_outgoing);
  T = 2*Z_outgoing / (Z_incoming + Z_outgoing);
end