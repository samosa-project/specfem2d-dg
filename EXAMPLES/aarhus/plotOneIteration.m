function [fig_field] = plotOneIteration(distMic_Spkr, OFD, IT, boxx, boxy, dx, dz, forceDGMesh, time, pre_t, fac_t, pre_p, fac_p, cax)
  % Read dumps.
%   disp(['[',mfilename,'] Reading dumps (may be long).']);
  [X, Y, V] = readDumpsUnique(OFD, IT, 0);

  % Select box.
%   disp(['[',mfilename,'] Selecting values within box.']);
  select = ((Y>=min(boxy)) & (Y<=max(boxy)) & (X>=min(boxx)) & (X<=max(boxx)));
  X = X(select);
  Y = Y(select);
  V.pre = V.pre(select);

  % Interpolate for plotting.
%   disp(['[',mfilename,'] Interpolating.']);
  [Xi, Yi, Vi] = interpDumps(X, Y, V, range(X)/dx, range(Y)/dz, forceDGMesh);

%   [pre_t, fac_t] = prefix_factor_values({time});
%   [pre_p, fac_p] = prefix_factor_values({Vi.pre(:)});

  % Plot.
  TIT_addendum = [' at $t=',sprintf('%.2f', fac_t*time),'$~',pre_t,'s'];
  TIT = ['Pressure Field ',TIT_addendum];
  CBYLAB_PRESPEC = ['$p''$ [',pre_p,'Pa]'];
  fig_field = figure('units','normalized','outerposition',[0,0,1,0.55]);
  tightAxes = tight_subplot(1, 1, [0,0], [0.17,0.08], [0.07, 0.06]);
  pcolor(Xi, Yi, fac_p*Vi.pre); hold on;
  scatter(0, 0, 400, 'marker', 'pentagram', 'markerfacecolor', [1,0.6,0]*0.8, 'markeredgecolor', 'k', 'linewidth', 2);
  scatter(distMic_Spkr, 0*distMic_Spkr, 200, 'marker', 'v', 'markerfacecolor', [0,1,0]*0.6, 'markeredgecolor', 'k', 'linewidth', 2);
  daspect([1,1,1]);
  hcb = colorbar(); ylabel(hcb, CBYLAB_PRESPEC, 'interpreter', 'latex');
%   caxis([-1,1]*fac_p*max(abs(Vi.pre(:))));
  caxis(fac_p*cax);
  % colormaps_fromPython('seismic', 1);
  thrsh = 0.5;
  lim = 0.05;
%   colormaps_custom([-1, -1+lim, -(1+thrsh)/2, -thrsh, 0, thrsh, (1+thrsh)/2, 1-lim, 1], [[0,0,1].*[0.25,0.75,1]';[0.9,0.9,1];[1,1,1];[1,0.9,0.9];[1,0,0].*[1,0.75,0.25]'], 1);
  colormaps_fromPython('seismic', 1);
  xlabel(['$x$ [m]']); ylabel(['$y$ [m]']);
  title(TIT);
end

