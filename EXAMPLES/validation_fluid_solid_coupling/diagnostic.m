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
debug_fig = 0;
summary_fig = 1;
normalise_ylims = 0;

classical1_zhang0 = 0; % theoretical values choice

mindelayincrefl = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format_positions = ['%',num2str(floor(log10(max(abs(STATPOS(:,2)))))+6)','.1f'];
vfmt = ['%11.4e'];
tfmt = ['%.4f'];
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


% Loop over couples of stations and compute ratios.
disp([' ']);
disp(['> Computing ratios.']);

if(summary_fig)
  figure('units','normalized','outerposition',[0,0,1,1]);
  tightAxes = tight_subplot(3, 2, [0.03,0.06], [0.11,0.07], [0.08, 0.01]);
  mapta = [1,3,5,2,4,6];
  facmag=1e6; unit='nm/s';
end

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
  i2 = snells(vp_1, vp_2, incident_angle); % deduce i2
  j2 = snells(vp_1, vs_2, incident_angle); % deduce j2
  
  % Compute experimental ratios.
  if(s2f1_or_f2s0)
%     error('not impl');
    incoming_reflected_vx = squeeze(amp(1, correspondingIDs(ID_inc),:)); % solid waves
    incoming_reflected_vz = squeeze(amp(2, correspondingIDs(ID_inc),:)); % solid waves
    incoming_reflected_t = time(correspondingIDs(ID_inc),:);
    transmitted_v = squeeze(amp(1,correspondingIDs(ID_out),:)); % pressure waves
    transmitted_t = time(correspondingIDs(ID_out),:);
    if(classical1_zhang0)
      incoming_reflected_v = (incoming_reflected_vx.^2+incoming_reflected_vz.^2).^0.5; % version not allowing to distinguish P from S
    else
      incoming_vp = incoming_reflected_vz; % incoming P-wave is along z axis
      [reflected_vp, reflected_vs] = vxz2vps(incoming_reflected_vx, incoming_reflected_vz, i2, j2);
    end
  else
    % FLUID-TO-SOLID
    incoming_reflected_v = squeeze(amp(1, correspondingIDs(ID_inc),:)); % pressure waves
    incoming_reflected_t = time(correspondingIDs(ID_inc),:);
    transmitted_vx = squeeze(amp(1,correspondingIDs(ID_out),:)); % solid waves
    transmitted_vz = squeeze(amp(2,correspondingIDs(ID_out),:)); % solid waves
    transmitted_t = time(correspondingIDs(ID_out),:);
    if(classical1_zhang0)
      transmitted_v = (transmitted_vx.^2+transmitted_vz.^2).^0.5; % version not allowing to distinguish P from S
    else
      [transmitted_vp, transmitted_vs] = vxz2vps(transmitted_vx, transmitted_vz, incident_angle, i2, j2);
      
      if(debug_fig); do_debug_figure;end

      if(summary_fig)
        axes(tightAxes(mapta(i)));
        [vang, vmag] = cart2pol(transmitted_vx, transmitted_vz);
        scatter(vang*180/pi, vmag*facmag, 50, transmitted_t, 'filled', 'displayname', ['ground velocity [',unit,']']); hold on;
        colormaps_fromPython('hsv', 1);
        plot((incident_angle*180/pi-90)*[1,1],ylim, 'displayname', '$n_2$');
        plot((incident_angle*180/pi-90-i2*180/pi)*[1,1],ylim, ':', 'displayname', '$i_2$');
        plot((incident_angle*180/pi-90-(j2*180/pi-90))*[1,1],ylim, ':', 'displayname', '$\left(j_2-90^\circ\right)$');
        if(ismember(mapta(i),[1,2]))
          title(['Incidence Angle = $',num2str(incident_angle*180/pi),'^\circ$']);
        end
        if(ismember(mapta(i),[5,6]))
          xlabel(['solid angle anti-clockwise from horizontal [$^\circ$]']);
        end
      end
    end
  end
  
  % Obtain theoretical values.
  if(s2f1_or_f2s0)
    if(classical1_zhang0)
      [R_th, T_th] = ReflexionTransmissionCoefs(Z_2P, Z_1, incident_angle, snells(vp_2, vp_1, incident_angle));
    else
      [R_th, T_th] = ReflexionTransmissionCoefsZhang(s2f1_or_f2s0, vp_1, Z_1, vp_2, vs_2, Z_2P, Z_2S, incident_angle);
      T_th = sum(T_th); % assume we are measuring both the P and the S
    end
    R_th = abs(R_th); % work in magnitude only
    % Rationale:
    % pt = T*pi, pr = R*pi
    % pt = Zs*vi => pt = T*Zs*vi => vi/pt = 1/(T*Zs)
    % (vi+vr) = (1+R)*vi => (vi+vr)/pt = T*Zs/(1+R)
%     Vi_over_Pt_th = 1/(T_s2f*Z_1); % v incident over p transmitted
    Pt_over_Vi_th = T_th*Z_1; % p transmitted over v incident
  %   ViVr_over_Pt = (1+R_s2f)/(T_s2f*Z_s); % If too close to interface, v transmitted over (p incident + p reflected)
%     disp(['CONVERSION COEFFICIENT NOT SURE']);
%     disp([fig_title, ' coupling, Z_s = ',num2str(Z_s), ', Z_f = ',num2str(Z_f),', T = ',num2str(T_s2f), ', R = ',num2str(R_s2f), '.']);
%     spatial_wavelength = spatial_wavelength_s;
%     localT_th = Vi_over_Pt_th; localT_th_name = 'Vi/Pt';
    T_th = Pt_over_Vi_th;
    R_th = R_th;
    localT_th_name = 'Pt/Vi';
    localR_th_name = 'Vr/Vi';
  else
    % FLUID-TO-SOLID
    if(classical1_zhang0)
      [R_th, T_th] = ReflexionTransmissionCoefs(Z_1, Z_2P, incident_angle, snells(vp_1, vp_2, incident_angle));
      % Rationale:
      % pt = Tpi, pr=Rpi
      % pt=vtZs => vtZs=Tpi => vt/pi=T/Zs
      % (pi+pr) = (1+R)pi => vt/(pi+pr)=T/(Zs*(1+R))
      Vt_over_Pi_th = T_th/Z_1; % v transmitted over p incident
    %   Vt_over_PiPr_th = T_f2s/(Z_s*(1+R_f2s)); % if too close to interface, v transmitted over (p incident + p reflected)
  %     disp([fig_title, ' coupling, Z_s = ',num2str(Z_s), ', Z_f = ',num2str(Z_f),', T = ',num2str(T_f2s), ', R = ',num2str(R_f2s), '.']);
  %     spatial_wavelength = spatial_wavelength_f;
      T_th = Vt_over_Pi_th;
      R_th = R_th;
      localT_th_name = 'Vt/Pi';
      localR_th_name = 'Pr/Pi';
    else
      [R_th, T_th] = ReflexionTransmissionCoefsZhang(s2f1_or_f2s0, vp_1, Z_1, vp_2, vs_2, Z_2P, Z_2S, incident_angle);
      VPSoP_th = T_th/Z_1;
      VPoP_th = VPSoP_th(1); localTP_th_name = 'VPt/Pi';
      VSoP_th = VPSoP_th(2); localTS_th_name = 'VSt/Pi';
      localR_th_name = 'Pr/Pi ';
    end
  end
  
  % Obtain experimental values and display.
  if(classical1_zhang0)
    [I_peak_t, I_peak_v] = findFirstPeak(incoming_reflected_t, incoming_reflected_v); % Find incoming wave.
    [T_peak_t, T_peak_v] = findFirstPeak(transmitted_t, transmitted_v); % Find transmitted wave.
    [R_peak_t, R_peak_v] = findFirstPeak(incoming_reflected_t(incoming_reflected_t>I_peak_t+mindelayincrefl),incoming_reflected_v(incoming_reflected_t>I_peak_t+mindelayincrefl)); % Find reflected wave.
    if(R_peak_t<I_peak_t); error(['good job, you found the reflection before the incident wave']); end
    T_exp = T_peak_v/I_peak_v;
    R_exp = R_peak_v/I_peak_v;
    %
  %     disp(['  > Couple at (x_1, x_2, z_1, z_2)=(',sprintf(format_positions,STATPOS(correspondingIDs,:)),'):']);
    disp(['  > Incident angle = ',sprintf('%3.0f',incident_angle*180/pi),'°. (Incident wave found at t=',sprintf(tfmt,I_peak_t),'s.)']);
    disp(['  > (t=',sprintf(tfmt,T_peak_t),'s) T^{',fig_title,'}_th = ',sprintf(vfmt,T_th),' m/s/Pa, T',fig_title,' = ',localT_th_name,' = ',sprintf(vfmt,T_exp),' m/s/Pa, relative error = ',sprintf('%6.3f',100*abs(T_exp-T_th)/abs(T_th)),'%.']);
    disp(['  > (t=',sprintf(tfmt,R_peak_t),'s) R^{',fig_title,'}_th = ',sprintf(vfmt,R_th),'       , R',fig_title,' = ',localR_th_name,' = ',sprintf(vfmt,R_exp),'       , relative error = ',sprintf('%6.3f',100*abs(R_exp-R_th)/abs(R_th)),'%.']);
  else
    if(s2f1_or_f2s0)
      %error()
    else
      [I_peak_t, I_peak_v] = findFirstPeak(incoming_reflected_t, incoming_reflected_v); % Find incoming pressure wave.
      [~, R_peak_v] = findFirstPeak(incoming_reflected_t(incoming_reflected_t>I_peak_t+mindelayincrefl),incoming_reflected_v(incoming_reflected_t>I_peak_t+mindelayincrefl)); % Find reflected wave.
      R_exp = R_peak_v/I_peak_v;
      
%       [TP_peak_t, TP_peak_v] = findFirstPeak(transmitted_t, transmitted_vp); % Find transmitted P-wave.
%       try [TS_peak_t, TS_peak_v] = findFirstPeak(transmitted_t, transmitted_vs); catch TS_peak_v=max(abs(transmitted_vs)); TS_peak_t=0; end % Find transmitted S-wave.
      TP_peak_v = range(transmitted_vp); % transmitted P-wave amplitude
      TS_peak_v = range(transmitted_vs); % transmitted S-wave amplitude
      VPoP_exp = TP_peak_v/I_peak_v;
      VSoP_exp = TS_peak_v/I_peak_v;
      
      disp(['  > R_{f}     = ',sprintf(vfmt,R_th),' | ',localR_th_name,' = ',sprintf(vfmt,R_exp),' | rel. err. = ',sprintf('%6.3f',100*abs(R_exp-R_th)/abs(R_th)),'%.']);
      disp(['  > T_{f->Pw} = ',sprintf(vfmt,T_th(1)),' | ',localTP_th_name,' = ',sprintf(vfmt,VPoP_exp),' | rel. err. = ',sprintf('%6.3f',100*abs(VPoP_exp-VPoP_th)/abs(VPoP_th)),'%.']);
      disp(['  > T_{f->Sw} = ',sprintf(vfmt,T_th(2)),' | ',localTS_th_name,' = ',sprintf(vfmt,VSoP_exp),' | rel. err. = ',sprintf('%6.3f',100*abs(VSoP_exp-VSoP_th)/abs(VSoP_th)),'%.']);
    end
  end
  disp(' ');
  % TODO: save agreements.
end

% TODO: plot of saved agreements.

if(summary_fig)
  linkaxes(tightAxes,'xy');
  set(tightAxes, 'xlim', [-1,1]*180, 'xtick', -180:45:180);
  set(tightAxes([1:end-2]), 'xticklabel', {});
  legend('location', 'southoutside', 'NumColumns', 4);
  colorbar;
end

diary off; % End logging.

% Function for Snell's law

% Function for reflexion and transmission coefficients computation.
