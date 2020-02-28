function [Xi, Yi, Zi] = interpDumps(X, Y, V, nx, ny)
  x_exact = linspace(min(X), max(X), nx);
  y_exact = linspace(min(Y), max(Y), ny);
  [Xi, Yi] = meshgrid(x_exact,y_exact);

  % Interpolate dump.
  F = scatteredInterpolant(X,Y,V);
  Zi = F(Xi, Yi);
end

