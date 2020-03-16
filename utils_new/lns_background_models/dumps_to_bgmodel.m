% Author:        Léo Martire.
% Description:   Loads dumps from SPECFEM and converts them to a LNS
%                generalised background model array. Shape is dictated by
%                the script 'order_bg_model.m'.
% Notes:         Needs scripts:
%                  utils_new/lns_background_models/order_bg_model.m
%                  utils_new/readDumpsUnique.m
%                  utils_new/interpDumps.m
%
% Usage:
%   [ROWS] = dumps_to_bgmodel(OFD, IT, uniform)
% with:
%   TODO.

function [ROWS, info] = dumps_to_bgmodel(OFD, IT, uniform)
  if(not(exist('uniform', 'var')))
    uniform = struct();
    uniform.do = 0;
  end
  
  addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new'));
  
  % Safety.
  if(not(OFD(end) == filesep))
    OFD = [OFD, filesep];
  end
  
  % Build paths.
  parfile = [OFD, 'input_parfile'];
  atmfile = [OFD, 'input_atmospheric_model.dat'];
  
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
  
  % Save/display information.
  info_i = 1;
  info{info_i} = ['[',mfilename,'] Constitutive variables [rho, vx, vz, p]: loaded from SPECFEM dumps located in ''',OFD,''', at IT=',num2str(IT),'.'];
  disp(info{info_i});
  
  % Read dumps and integrate those.
  % TODO: read rho, vx, vz.
  [X, Z, VALUES] = readDumpsUnique(OFD, IT, 0);
  rho0 = VALUES.rho;
  vx = VALUES.vel(:,1);
  vz = VALUES.vel(:,2);
  pre = VALUES.pre;
  Z = Z-min(Z); % recall background models should be put in ASL format along z (i.e. min(z) should be 0)
  
  % Choose to interpolate.
  if(uniform.do)
    % TODO: do it for rho, vx, vz.
    [Xi, Yi, pr_i] = interpDumps(X, Z, VALUES, uniform.nx, uniform.nz);
    ninterp = numel(Xi);
    X = reshape(Xi, ninterp, 1);
    Z = reshape(Yi, ninterp, 1);
    VALUES = reshape(pr_i, ninterp, 1);
  end
  
  % Recover physical quantities (gravity, gamma, mu, kappa).
  MODEL = readExampleFiles_extractParam(parfile, 'MODEL', 'string');
  switch(MODEL)
    case 'default'
      mu = readExampleFiles_extractParam(parfile, 'dynamic_viscosity', 'float');
      kap = readExampleFiles_extractParam(parfile, 'thermal_conductivity', 'float');
      cp = readExampleFiles_extractParam(parfile, 'constant_p', 'float');
      cv = readExampleFiles_extractParam(parfile, 'constant_v', 'float');
      gamma = cp/cv;
      
      USE_ISOTHERMAL_MODEL = readExampleFiles_extractParam(parfile, 'USE_ISOTHERMAL_MODEL', 'bool');
      USE_LNS = readExampleFiles_extractParam(parfile, 'USE_LNS', 'bool');
      
      if(USE_ISOTHERMAL_MODEL)
        grav = readExampleFiles_extractParam(parfile, 'gravity', 'float');
      else
        if(USE_LNS)
          error(['[',mfilename,', ERROR] Not implemented yet.']);
          grav = 0; % Isobaric, needs to be zero, at least in FNS.
        else
          grav = 0; % Isobaric, needs to be zero, at least in FNS.
        end
      end
      
      info_i = info_i + 1;
      info{info_i} = ['[',mfilename,'] Physical quantities [gravity, gamma, mu, kappa]: MODEL=default, hence extracted all from parfile ''',parfile,'''.'];
      disp(info{info_i});
    
    case 'external_DG'
      % Load the external_DG atmospheric model. It is in ASL convention, no need to change vertical reference.
      [ZDump, ~, ~, ~, ~, ~, GRAVDump, ~, KAPPADump, MUDump, ~, ~, ~, ~, ~, ~, GAMMADump, ~, ~] = extract_atmos_model(atmfile, 3, 0, -1);
      interpmethod = 'spline';
      grav = interp1(ZDump, GRAVDump, Z, interpmethod);
      gamma = interp1(ZDump, GAMMADump, Z, interpmethod);
      mu = interp1(ZDump, MUDump, Z, interpmethod);
      kap = interp1(ZDump, KAPPADump, Z, interpmethod);
      
      info_i = info_i + 1;
      info{info_i} = ['[',mfilename,'] Physical quantities [gravity, gamma, mu, kappa]: interpolated from external_DG file ''',atmfile,''', using Matlab''s interp1 with ''',interpmethod,''' method.'];
      disp(info{info_i});
      
    otherwise
      error(['[',mfilename,', ERROR] This MODEL (',MODEL,') was not implemented yet.']);
  end
%   pause
  
  % Put in lexicographic order, for readability.
  [~, isort] = sortrows([X, Z]);
  X = X(isort);
  Z = Z(isort);
  rho0 = rho0(isort);
  vx = vx(isort);
  vz = vz(isort);
  pre = pre(isort);
  
  % handmade tweaks for tests
  if(tweaks)
    pre = min(pre) + factor_dp*(pre-min(pre));
  end
  
  q=X; disp(['[',mfilename,'] (min, max) of X:    (',num2str(min(q)),', ',num2str(max(q)),')']);
  q=Z; disp(['[',mfilename,'] (min, max) of Z:    (',num2str(min(q)),', ',num2str(max(q)),')']);
  q=rho0; disp(['[',mfilename,'] (min, max) of rho0: (',num2str(min(q)),', ',num2str(max(q)),')']);
  q=pre; disp(['[',mfilename,'] (min, max) of pre:  (',num2str(min(q)),', ',num2str(max(q)),')']);
  
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
%   toplot = sqrt(gamma*pre/rho0);
%   [Xi, pr_i, VALi] = interpDumps(X, Z, toplot, 1000, 1000); fh = figure(); pcolor(Xi, pr_i, VALi); shading interp; colormaps_fromPython('seismic', 1); colorbar;
%   pause;
end

