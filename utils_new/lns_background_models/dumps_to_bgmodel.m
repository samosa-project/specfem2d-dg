function [ROWS] = dumps_to_bgmodel(OFD, IT, uniform)
  % Artificial modifications, tweaks.
  tweaks = 1;
  factor_dp = 2e6; % increase the pressure perturbation by a factor
  
  % Get quantities' order, to put model in right order.
  [order, ~] = order_bg_model();
  nb_qty = size(order, 1);
  
  % Default values.
  rho0 = 1.4;
  vx = 0;
  vz = 0;
  grav = 0;
  gamma = 1005/717.3447537473;
  pre = 340^2*rho0/gamma;
  mu = 0;
  kap = 0;
  
  % Read dumps and integrate those.
  [X, Z, pre, ~] = readDumpsUnique(OFD, IT, 0);
  Z = Z-min(Z); % recall background models should be put in ASL format along z (i.e. min(z) should be 0)
  
  % choose to interpolate
  if(uniform.do)
    [Xi, Yi, pr_i] = interpDumps(X, Z, pre, uniform.nx, uniform.nz);
    ninterp = numel(Xi);
    X = reshape(Xi, ninterp, 1);
    Z = reshape(Yi, ninterp, 1);
    pre = reshape(pr_i, ninterp, 1);
  end
  
  % lexicographic order, for readability
  [~, isort] = sortrows([X, Z]);
  X = X(isort);
  Z = Z(isort);
  pre = pre(isort);
  
  % handmade tweaks for tests
  if(tweaks)
    pre = min(pre) + factor_dp*(pre-min(pre));
  end
  
  % Format under large unit table.
  ROWS = zeros(numel(X), nb_qty);
  curQty = 'xx'; id = find(all((order==curQty)')); ROWS(:, id) = X;
  curQty = 'zz'; id = find(all((order==curQty)')); ROWS(:, id) = Z;
  curQty = 'rh'; id = find(all((order==curQty)')); ROWS(:, id) = rho0;
  curQty = 'vx'; id = find(all((order==curQty)')); ROWS(:, id) = vx;
  curQty = 'vz'; id = find(all((order==curQty)')); ROWS(:, id) = vz;
  curQty = 'pr'; id = find(all((order==curQty)')); ROWS(:, id) = pre;
  curQty = 'gr'; id = find(all((order==curQty)')); ROWS(:, id) = grav;
  curQty = 'ga'; id = find(all((order==curQty)')); ROWS(:, id) = gamma;
  curQty = 'mu'; id = find(all((order==curQty)')); ROWS(:, id) = mu;
  curQty = 'ka'; id = find(all((order==curQty)')); ROWS(:, id) = kap;
  
  % Some visualisation.
%   figure(); spy(ROWS); daspect([1,numel(X)/nb_qty,1]); xticks(1:nb_qty); xticklabels(order);
%   toplot = VAL;
  toplot = sqrt(gamma*pre/rho0);
  [Xi, pr_i, VALi] = interpDumps(X, Z, toplot, 1000, 1000); fh = figure(); pcolor(Xi, pr_i, VALi); shading interp; colormaps_fromPython('seismic', 1); colorbar;
end

