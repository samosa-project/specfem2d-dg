function [] = write_bg_model(OFD, IT)
  nb_significantdigits = 16;
  do_interp = 0;
  do_interp_NINTERPX = 4;
  do_interp_NINTERPZ = do_interp_NINTERPX;
  
  output_file = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_lns_custom_wavefield_fromdg/background_model.dat';
%   output_file = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_lns_custom_wavefield_4elements/background_model.dat';

  tweaks = 1;
  factor_dp = 2e6; % increase the pressure perturbation by a factor
  
  % get quantities to put in model
  [order, lab] = order_bg_model();
  nb_qty = size(order, 1);
  
  % prepare format
  spacing = nb_significantdigits+7;
  format_number = ['%',num2str(spacing),'.',num2str(nb_significantdigits),'e '];
  format_line = repmat(format_number, [1,nb_qty]);
  format_line = [format_line(1:end-1), '\n'];
  
  % default values to fill out blanks
  rho0 = 1.4;
  vx = 0;
  vz = 0;
  gra = 0;
  gam = 1005/717.3447537473;
  mu = 0;
  ka = 0;
  
  % read some pressure dumps and integrate those within a default background model
  [X, Z, VAL, imagetype_wavefield_dumps] = readDumpsUnique(OFD, IT, 0);
  Z = Z-min(Z); % recall background models should be put in ASL format along z (i.e. min(z) should be 0)
  
  % choose to interpolate
  if(do_interp)
    [Xi, Yi, Zi] = interpDumps(X, Z, VAL, do_interp_NINTERPX, do_interp_NINTERPZ);
    ninterp = numel(Xi);
    X = reshape(Xi,ninterp,1);
    Z = reshape(Yi,ninterp,1);
    VAL = reshape(Zi,ninterp,1);
  end
  
  % lexicographic order, for readability
  [~, isort] = sortrows([X, Z]);
  X = X(isort);
  Z = Z(isort);
  VAL = VAL(isort);
  
  % handmade tweaks for tests
  if(tweaks)
    VAL = min(VAL) + factor_dp*(VAL-min(VAL));
  end
  
%   X = X(1:10);
%   Y = Y(1:10);
%   V = V(1:10);
  
  ROWS = zeros(numel(X), nb_qty);
  
  curQty = 'xx'; id = find(all((order==curQty)')); ROWS(:, id) = X;
  curQty = 'zz'; id = find(all((order==curQty)')); ROWS(:, id) = Z;
  curQty = 'pr'; id = find(all((order==curQty)')); ROWS(:, id) = VAL;
  
  curQty = 'rh'; id = find(all((order==curQty)')); ROWS(:, id) = rho0;
  curQty = 'vx'; id = find(all((order==curQty)')); ROWS(:, id) = vx;
  curQty = 'vz'; id = find(all((order==curQty)')); ROWS(:, id) = vz;
  curQty = 'gr'; id = find(all((order==curQty)')); ROWS(:, id) = gra;
  curQty = 'ga'; id = find(all((order==curQty)')); ROWS(:, id) = gam;
  curQty = 'mu'; id = find(all((order==curQty)')); ROWS(:, id) = mu;
  curQty = 'ka'; id = find(all((order==curQty)')); ROWS(:, id) = ka;
  
%   % visualisation
%   figure(); spy(ROWS); daspect([1,numel(X)/nb_qty,1]); xticks(1:nb_qty); xticklabels(order);
%   toplot = VAL;
  toplot = sqrt(gam*VAL/rho0);
  [Xi, Zi, VALi] = interpDumps(X, Z, toplot, 1000, 1000); fh = figure(); pcolor(Xi, Zi, VALi); shading interp; colormaps_fromPython('seismic', 1); colorbar;
  pause
  
  % Open file and write in it.
  foutput = fopen(output_file, 'w');
  
  % header
  fprintf(foutput, 'custom/default LNS generalised background model\n');
  fprintf(foutput, ['using dumps from IT=',num2str(IT),' in ',OFD,'\n']);
  for q = 1:nb_qty
    fprintf(foutput, pad(lab{q}, spacing+1));
  end
  fprintf(foutput, '\n');
  fprintf(foutput, format_line, ROWS');
  
  fclose('all');
  disp(['[',mfilename,'] Finished writing background model to file (''',output_file,''').']);
  
end

