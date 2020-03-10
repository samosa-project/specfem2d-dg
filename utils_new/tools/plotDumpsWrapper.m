function [fh] = plotDumpsWrapper(OFD, IT, verbose, nx, ny)
  [X,Y,V, ~] = readDumpsUnique(OFD, IT, verbose);
  [Xi, Yi, Zi] = interpDumps(X, Y, V, nx, ny);
  
  fh = figure();
  pcolor(Xi, Yi, Zi);
  shading interp;
  colormaps_fromPython('seismic', 1);
  
  caxis([-1,1]*max(abs(Zi(:))));
  h = colorbar;
end