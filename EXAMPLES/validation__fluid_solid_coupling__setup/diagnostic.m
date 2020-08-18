clear all;
close all;
clc;

LWangles = 3; LSangles = ':'; colangles = get(0,'defaultAxesColorOrder');
mS_radialplot = 20;
TLAB = ['time [s]'];
CMAP = colormaps_fromPython('jet', 0);

setup;

% station_ids = 1:6;
station_ids = [2, 5];

% for i = 1:numel(cases)
for i = 1:2
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
    TIT = [TIT, 'STF'];
  else
    TIT = [TIT, 'FTS'];
    ylim_pre = [-0.25, 3];
    if(cases{i}.ortho0_slant1)
      ylim_vel = [-3, 4]*1e-6;
    else
      ylim_vel = [-5, 1]*1e-6;
    end
  end
  TIT = [TIT, ', '];
  if(cases{i}.ortho0_slant1)
    TIT = [TIT, '$\theta_i\neq0$'];
    theta_i = ic;
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
  
  [rhof, c, rhos, vp, vs] = get_models(parfile);
  [i1_i2_j2t_j2r] = get_predicted_angles_deg(ic, c, vp, vs)* pi/180; % cast to [rad]
  if(cases{i}.fts0_stf1)
    % STF
    if(cases{i}.ortho0_slant1)
      selangles = [pi/2-theta_i, theta_i-pi/2, i1_i2_j2t_j2r(4)]; % incident P @(pi/2-ic), reflected P @(ic-pi/2), reflected S @()
      angdnam = {'incident P', 'reflected P', 'reflected S'};
    else
      selangles = [pi/2-theta_i]; % incident P and reflected P @(pi/2-ic=0)
      angdnam = {'incident \& reflected P'};
    end
  else
    % FTS
    [R, T] = ReflexionTransmissionCoefsZhang(cases{i}.fts0_stf1, c, c*rhof, vp, vs, vp*rhos, vs*rhos, theta_i)
    [~, incident_P] = findFirstPeak(time, dp); % get incident P pressure wave
    expected_reflected_amplitude = R * incident_P; % expected P- and S-wave amplitude
    expected_transmitted_amplitude = (T/(c*rhof)) * incident_P; % expected P- and S-wave amplitude
    if(cases{i}.ortho0_slant1)
      selangles = [i1_i2_j2t_j2r(2)-pi/2, i1_i2_j2t_j2r(3)]; % transmitted P-wave, transmitted S-wave
      angdnam = {'transmitted P', 'transmitted S'};
    else
      selangles = [theta_i - pi/2];
      %[~, transmitted_P] = findFirstPeak(time, abs(vz)); % get transmitted P seismic wave
      angdnam = {'transmitted P'};
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Do the plot.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure('units','normalized','outerposition',[(i-1)*0.15, 0, 0.5, 1]);
  tightAxes = tight_subplot(3, 1, [0.16, 0.], [0.1, 0.06], [0.1, 0.09]); % gap, marg_h, marg_w
  
  % air time series
  axes(tightAxes(1));
  [prefix, factor] = prefix_factor_values({dp});
  h=[];
  h=[h, plot(time, factor*dp, 'displayname', '$p''$')]; hold on;
%   scatter(time(istatloc, :), time(istatloc, :)*0+min(ylim), 100, time(istatloc, :), 'filled'); hold on;
  if(cases{i}.fts0_stf1==0)
    ttt = [min(time), max(time)];
    h=[h, plot(ttt, expected_reflected_amplitude*[1,1], 'displayname', 'reflected P', 'linewidth', LWangles, 'linestyle', LSangles)]; hold on;
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
  h=[h, plot(time, factor*vx, 'displayname', '$v_x$')]; hold on;
  h=[h, plot(time, factor*vz, 'displayname', '$v_z$')]; hold on;
  scatter(time, time*0+factor*min(ylim_vel), 80, time, 'filled'); hold on;
  legend(h);
  xlabel(TLAB);
  ylabel(['$v_{x,z}$ [',prefix,'m/s]']);
  ylim(factor*ylim_vel);

  % solid radial plot
  axes(tightAxes(3))
  % plot angles
  h=[];
  for ia = 1:numel(selangles)
    angl = selangles(ia);
    [x, y] = pol2cart(angl*[1,1], factor*1.5*max(amp(:))*[-1,1]);
    h = [h, plot(x, y, 'color', colangles(ia,:), 'linewidth', LWangles, 'linestyle', LSangles, 'displayname', angdnam{ia})]; hold on;
    if(expected_transmitted_amplitude(ia)~=0)
      ha = draw_ampl_circle(factor*expected_transmitted_amplitude(ia), angl); set(ha, 'color', colangles(ia,:), 'linewidth', LWangles, 'linestyle', LSangles);
    end
  end
  % plot time series
  h = [h, scatter(factor*vx, factor*vz, mS_radialplot, time, 'filled', 'displayname', 'ground motion')]; hold on;
  legend(h, 'location', 'westoutside');
  hcb = colorbar;
  ylabel(hcb, TLAB, 'interpreter', 'latex', 'fontsize', 26);
  [~, r] = cart2pol(vx*factor, vz*factor);
  xlim([-1,1]*max(r)*1.2);
  ylim([-1,1]*max(r)*1.2);
  xlabel(['$v_{x}$ [',prefix,'m/s]']);
  ylabel(['$v_{z}$ [',prefix,'m/s]']);
  daspect([1,1,1]);

  % adjust
  set(tightAxes, 'colormap', CMAP, 'clim', cases{i}.tlim);
  set(tightAxes([1,2]), 'xlim', cases{i}.tlim);
  vshift = 0.035;
  moveup = (tightAxes(1).Position([2])-tightAxes(1).Position([4])-vshift) - tightAxes(2).Position(2);
  tightAxes(2).Position = tightAxes(2).Position + [0, moveup, 0, 0];
  tightAxes(3).Position = tightAxes(3).Position + [0, 0, 0, moveup+2*vshift];
end