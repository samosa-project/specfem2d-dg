clear all;
close all;
clc;

LWangles = 3; LSangles = ':'; colangles = get(0,'defaultAxesColorOrder');
mS_radialplot = 20;
TLAB = ['time [s]'];
CMAP = colormaps_fromPython('jet', 0);

% Naming. Make agree with macros in TeX.
name_Pd = '\mathrm{P}_1^{\mathrm{d}}'; % incident P from fluid side
name_Pu = '\mathrm{P}_2^{\mathrm{u}}'; % incident P from solid side
name_PuPu = ['$Z_1T^{\mathrm{SF}}_{\mathrm{PP}}',name_Pu,'$']; % STF PTP transmission
name_PuPd = ['$R^{\mathrm{S}}_{\mathrm{PP}}',name_Pu,'$']; % STF PRP reflection
name_PuSd = ['$R^{\mathrm{S}}_{\mathrm{PS}}',name_Pu,'$']; % STF PRS reflection
name_PdPd = ['$T^{\mathrm{FS}}_{\mathrm{PP}}',name_Pd,'$']; % FTS PTP transmission
name_PdSd = ['$T^{\mathrm{FS}}_{\mathrm{PS}}',name_Pd,'$']; % FTS PTS transmission
name_PdPu = ['$R^{\mathrm{F}}_{\mathrm{PP}}',name_Pd,'$']; % FTS PRP reflection

setup;

% station_ids = 1:6;
% station_ids = [2, 5]-1; % left one
station_ids = [2, 5]; % middle one
% station_ids = [2, 5]+1; % right one

effective_R = {};
effective_T = {};
saveR = {};
saveT = {};
err_R = {};
err_T = {};

for i = 1:numel(cases)
% for i = 1:2
% for i = 3
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get working folders and
  % load synthetics.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  folder = folderz{i};
  OFD = [folder, filesep, 'OUTPUT_FILES', filesep];
  parfile = [OFD, 'input_parfile'];
  [time, amp] = load_synthetics(OFD, parfile, station_ids); % amp has shape [2, nstat, ntimes]
  sel = (time(1, :) >= min(cases{i}.tlim)) & (time(1, :) <= max(cases{i}.tlim));
  newtime = time(:, sel);
  newamp = amp(:, :, sel);
  time = newtime;
  amp = newamp;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Prepare.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TIT = '';
  if(cases{i}.fts0_stf1)
    TIT = [TIT, stftag];
    rad_lim = 0.3e-3;
    if(cases{i}.ortho0_slant1)
      ylim_pre = [-0.05, 0.5];
      ylim_vel = [-0.3, 0.3]*1e-3;
    else
      ylim_pre = [-0.05, 0.5];
      ylim_vel = [-0.2, 0.3]*1e-3;
    end
  else
    TIT = [TIT, ftstag];
    ylim_pre = [-0.25, 3];
    rad_lim = 1.5e-6;
    if(cases{i}.ortho0_slant1)
      ylim_vel = [-2.5, 2]*1e-6;
    else
      ylim_vel = [-2.5, 0.5]*1e-6;
    end
  end
  TIT = [TIT, ', '];
  if(cases{i}.ortho0_slant1)
    TIT = [TIT, '$\theta_i\neq0$'];
    theta_i = ic_rad;
  else
    TIT = [TIT, '$\theta_i=0$'];
    theta_i = 0;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get angles and ratios.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %extract time series
  time = time(1,:);
  dp = squeeze(amp(1, 1, :));
  vx = squeeze(amp(1, 2, :));
  vz = squeeze(amp(2, 2, :));
  
  [rho__1, alpha__1, rho__2, alpha__2, beta__2] = get_models(parfile);
  [i1_i2_j2t_j2r] = get_predicted_angles_deg(ic_rad, alpha__1, alpha__2, beta__2)* pi/180; % cast to [rad]
%   [R, T] = ReflexionTransmissionCoefsZhang(cases{i}.fts0_stf1, alpha__1, rho__1, alpha__2, beta__2, rho__2, theta_i);
  [R, T] = RTCoefs(cases{i}.fts0_stf1, alpha__1, rho__1, alpha__2, beta__2, rho__2, theta_i);
  if(cases{i}.fts0_stf1)
    % STF
    [~, vmag] = cart2pol(vx, vz); % compute magnitude of displacement
%     [~, incident_P] = findFirstPeak(time, vmag); % get incident P seismic wave
    if(cases{i}.ortho0_slant1)
      selangles = [pi/2-theta_i, theta_i-pi/2, i1_i2_j2t_j2r(4)]; % incident P @(pi/2-ic), reflected P @(ic-pi/2), reflected S @()
      angdnam = {['$',name_Pu,'$'], name_PuPd, name_PuSd};
      [iP_rP_rS] = find_max_amplitude_along_direction(vx, vz, selangles); incident_P=iP_rP_rS(1); reflected_P=iP_rP_rS(2); reflected_S=iP_rP_rS(3);
    else
      selangles = [theta_i-pi/2] * [1, 1]; % incident P seismic wave and reflected P seimic wave
      angdnam = {['$',name_Pu,'$'], name_PuPd};
      [iP_rP] = find_max_amplitude_along_direction(vx, vz, selangles); incident_P=iP_rP(1); reflected_P=iP_rP(2); reflected_S=0;
    end
    [~, transmitted_P] = findFirstPeak(1:numel(dp), dp); transmitted_P=transmitted_P/(alpha__1*rho__1); % find transmitted P-wave in terms of velocity
    expected_reflected_amplitude = R * incident_P;
    expected_transmitted_amplitude = (alpha__1*rho__1*T) * incident_P; % in terms of pressure (pay attention to the multiplication by the fluid impedance)
    effective_R{i} = [reflected_P, reflected_S]/incident_P;
    effective_T{i} = transmitted_P/incident_P;
  else
    % FTS
    [~, incident_P] = findFirstPeak(time, dp); % incident P airwave
    expected_reflected_amplitude = R * incident_P; % expected reflected P airwave
    expected_transmitted_amplitude = (T/(alpha__1*rho__1)) * incident_P; % expected P and S seismic waves (pay attention to the conversion to velocity using the multiplication by the fluid impedance)
    incident_P = incident_P/(alpha__1*rho__1); % convert incident to velocity for later treatments
    if(cases{i}.ortho0_slant1)
      selangles = [i1_i2_j2t_j2r(2)-pi/2, i1_i2_j2t_j2r(3)]; % transmitted P seismic wave, transmitted S seismic wave
      angdnam = {name_PdPd, name_PdSd};
      [tP_tS] = find_max_amplitude_along_direction(vx, vz, selangles); transmitted_P=tP_tS(1); transmitted_S=tP_tS(2);
    else
      selangles = [theta_i - pi/2]; % transmitted P seismic wave
      %[~, transmitted_P] = findFirstPeak(time, abs(vz)); % get transmitted P seismic wave
      angdnam = {name_PdPd};
      [transmitted_P] = find_max_amplitude_along_direction(vx, vz, selangles); transmitted_S=0;
    end
    [~, reflected_P] = findFirstPeak(time, flip(dp)); reflected_P=reflected_P/(alpha__1*rho__1); % find reflected P-wave in terms of velocity
    
    effective_R{i} = reflected_P/incident_P;
    effective_T{i} = [transmitted_P, transmitted_S]/incident_P;
  end
  saveR{i} = abs(R);
  saveT{i} = abs(T);
  err_R{i} = abs(abs(effective_R{i}) - abs(R))./abs(R);
  err_T{i} = abs(abs(effective_T{i}) - abs(T))./abs(T);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Do the plot.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  THEFIGURE = figure('units','normalized','outerposition',[(i-1)*0.15, 0, 0.5, 1]);
  tightAxes = tight_subplot(3, 1, [0.16, 0.], [0.1, 0.06], [0.11, 0.09]); % gap, marg_h, marg_w
  
  % air time series
  axes(tightAxes(1));
  [prefix, factor] = prefix_factor_values({dp});
  h=[];
  h=[h, plot(time, factor*dp, 'displayname', 'synthetic $p''$')]; hold on;
%   scatter(time(istatloc, :), time(istatloc, :)*0+min(ylim), 100, time(istatloc, :), 'filled'); hold on;
  ttt = [min(time), max(time)];
  if(cases{i}.fts0_stf1)
    % STF: add expected transmitted airwave amplitude
    h = [h, plot(ttt, expected_transmitted_amplitude*[1,1], 'displayname', name_PuPu, 'linewidth', LWangles, 'linestyle', LSangles)]; hold on;
  else
    % FTS: add expected reflected airwave amplitude
    h = [h, plot(ttt, expected_reflected_amplitude(1)*[1,1], 'displayname', name_PdPu, 'linewidth', LWangles, 'linestyle', LSangles)]; hold on;
  end
  ylabel(['$p''$ [',prefix,'Pa]']);
  xticklabels({});
  title(TIT);
  ylim(factor*ylim_pre);
  legend('location', 'best');

  % solid time series
  axes(tightAxes(2))
  [prefix, factor] = prefix_factor_values({vx, vz});
  h=[];
  h=[h, plot(time, factor*vx, 'displayname', 'synthetic $v_x$')]; hold on;
  h=[h, plot(time, factor*vz, 'displayname', 'synthetic $v_z$')]; hold on;
  scatter(time, time*0+factor*min(ylim_vel), 80, time, 'filled'); hold on;
  legend(h, 'location', 'south', 'numcolumns', 2);
  xlabel(TLAB); ylabel(['$v_{x,z}$ [',prefix,'m/s]']);
  ylim(factor*ylim_vel);

  % solid radial plot
  axes(tightAxes(3))
  % plot angles
  h=[];
  for ia = 1:numel(selangles)
    angl = selangles(ia);
    [x, y] = pol2cart(angl*[1,1], factor*2*rad_lim*[-1,1]);
    if(cases{i}.fts0_stf1)
      % STF
      h = [h, plot(x, y, 'color', colangles(ia,:), 'linewidth', LWangles, 'linestyle', LSangles, 'displayname', angdnam{ia})]; hold on;
      if(ia>1 & expected_reflected_amplitude(ia-1)~=0) % ia==1 reserved for incident P wave in this case
        ha = draw_ampl_circle(factor*expected_reflected_amplitude(ia-1), angl, 0.05); set(ha, 'color', colangles(ia,:), 'linewidth', LWangles, 'linestyle', LSangles);
      end
    else
      % FTS
      h = [h, plot(x, y, 'color', colangles(ia+1,:), 'linewidth', LWangles, 'linestyle', LSangles, 'displayname', angdnam{ia})]; hold on;
      if(expected_transmitted_amplitude(ia)~=0)
        ha = draw_ampl_circle(factor*expected_transmitted_amplitude(ia), angl, 0.25); set(ha, 'color', colangles(ia+1,:), 'linewidth', LWangles, 'linestyle', LSangles);
      end
    end
  end
  % plot time series
  h = [h, scatter(factor*vx, factor*vz, mS_radialplot, time, 'filled', 'displayname', 'ground motion')]; hold on;
  legend(h, 'location', 'westoutside');
  hcb = colorbar; ylabel(hcb, TLAB, 'interpreter', 'latex', 'fontsize', 26);
  xlim([-1,1]*factor*rad_lim); ylim([-1,1]*factor*rad_lim);
  xlabel(['$v_{x}$ [',prefix,'m/s]']); ylabel(['$v_{z}$ [',prefix,'m/s]']);
  daspect([1,1,1]);
  
  % adjust
  set(tightAxes, 'colormap', CMAP, 'clim', cases{i}.tlim);
  set(tightAxes([1,2]), 'xlim', cases{i}.tlim);
  vshift = 0.035;
  moveup = (tightAxes(1).Position([2])-tightAxes(1).Position([4])-vshift) - tightAxes(2).Position(2);
  tightAxes(2).Position = tightAxes(2).Position + [0, moveup, 0, 0];
  tightAxes(3).Position = tightAxes(3).Position + [0, 0, 0, moveup+2*vshift-0.02];
  ll = add_labels_subplots(gcf, 1.8, 0, [0,0], 'roman'); set(ll(3),'Position',get(ll(3),'Position')+[0.012,0,0,0]);
  
  % save
  customSaveFig(THEFIGURE, [plotFolder,filesep,cases{i}.code], extToSave, 9999);
end

printTexTable(ic_deg, cases, saveR, saveT, effective_R, effective_T, err_R, err_T, stftag, ftstag);