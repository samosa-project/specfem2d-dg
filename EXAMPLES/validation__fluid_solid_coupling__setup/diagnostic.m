clear all;
close all;
clc;

setup;

% station_ids = 1:6;
station_ids = [2, 5];

for i = 1:numel(cases)
  % Get working folders.
  folder = folderz{i};
  OFD = [folder, filesep, 'OUTPUT_FILES', filesep];
  parfile = [OFD, 'input_parfile'];
  [time, amp] = load_synthetics(OFD, parfile, station_ids); % amp has shape [2, nstat, ntimes]
  
  sel = (time(1, :) >= min(cases{i}.tlim)) & (time(1, :) <= max(cases{i}.tlim));
  newtime = time(:, sel);
  newamp = amp(:, :, sel);
  time = newtime;
  amp = newamp;

  % plot
  LWangles = 2; LSangles = '-';
  TLAB = ['time [s]'];
  CMAP = colormaps_fromPython('jet', 0);
  figure('units','normalized','outerposition',[(i-1)*0.15, 0, 0.5, 1]);
  tightAxes = tight_subplot(3, 1, [0.16, 0.], [0.1, 0.06], [0.125, 0.085]); % gap, marg_h, marg_w

  % air time series
  axes(tightAxes(1));
  istatloc = 1;
  dp = squeeze(amp(1, istatloc, :));
  [prefix, factor] = prefix_factor_values({dp});
  h=[];
  h=[h, plot(time(istatloc, :), factor*dp, 'displayname', '$p''$')]; hold on;
%   scatter(time(istatloc, :), time(istatloc, :)*0+min(ylim), 100, time(istatloc, :), 'filled'); hold on;
  ylabel(['$p''$ [',prefix,'Pa]']);
  xticklabels({});
  TIT = '';
  if(cases{i}.fts0_stf1)
    TIT = [TIT, 'STF'];
  else
    TIT = [TIT, 'FTS'];
  end
  TIT = [TIT, ', '];
  if(cases{i}.ortho0_slant1)
    TIT = [TIT, '$\theta_i\neq0$'];
    theta_i = ic;
  else
    TIT = [TIT, '$\theta_i=0$'];
    theta_i = 0;
  end
  title(TIT);

  % solid time series
  axes(tightAxes(2))
  istatloc = 2;
  vx = squeeze(amp(1, istatloc, :));
  vz = squeeze(amp(2, istatloc, :));
  [prefix, factor] = prefix_factor_values({vx, vz});
  h=[];
  h=[h, plot(time(istatloc, :), factor*vx, 'displayname', '$v_x$')]; hold on;
  h=[h, plot(time(istatloc, :), factor*vz, 'displayname', '$v_z$')]; hold on;
  scatter(time(istatloc, :), time(istatloc, :)*0+min(ylim), 100, time(istatloc, :), 'filled'); hold on;
  legend(h);
  xlabel(TLAB);
  ylabel(['$v_{x,z}$ [',prefix,'m/s]']);

  % solid radial plot
  axes(tightAxes(3))
  % plot angles
  [i1_i2_j2t_j2r] = get_predicted_angles(parfile, ic) * pi/180;
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
    if(cases{i}.ortho0_slant1)
      selangles = [i1_i2_j2t_j2r(2) - pi/2, i1_i2_j2t_j2r(3)];
      angdnam = {'transmitted P', 'transmitted S'};
    else
      selangles = [pi/2-theta_i];
      angdnam = {'transmitted P'};
    end
  end
  for ia = 1:numel(selangles)
    a = selangles(ia);
    [x, y] = pol2cart(a*[1,1], factor*1.5*max(amp(:))*[-1,1]);
    plot(x, y, 'linewidth', LWangles, 'linestyle', LSangles, 'displayname', angdnam{ia}); hold on;
  end
  % plot time series
  scatter(factor*vx, factor*vz, 40, time(istatloc, :), 'filled', 'displayname', 'ground motion'); hold on;
  legend('location', 'best');
  hcb = colorbar;
  ylabel(hcb, TLAB, 'interpreter', 'latex', 'fontsize', 26);
  [~, r] = cart2pol(vx*factor, vz*factor);
  xlim([-1,1]*max(r));
  ylim([-1,1]*max(r));
  xlabel(['$v_{x}$ [',prefix,'m/s]']);
  ylabel(['$v_{z}$ [',prefix,'m/s]']);

  % adjust
  set(tightAxes, 'colormap', CMAP, 'clim', cases{i}.tlim);
  set(tightAxes([1,2]), 'xlim', cases{i}.tlim);
  
  vshift = 0.03;
  moveup = (tightAxes(1).Position([2])-tightAxes(1).Position([4])-vshift) - tightAxes(2).Position(2);
  tightAxes(2).Position = tightAxes(2).Position + [0, moveup, 0, 0];
  tightAxes(3).Position = tightAxes(3).Position + [0, 0, 0, moveup+2*vshift];
end