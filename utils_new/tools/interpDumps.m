function [Xi, Yi, Vi] = interpDumps(X, Y, V, nx, ny)
  tags = {'rho', 'vel', 'pre'};
  
  x_exact = linspace(min(X), max(X), nx);
  y_exact = linspace(min(Y), max(Y), ny);
  [Xi, Yi] = meshgrid(x_exact,y_exact);

  % Interpolate dump, on each non-empty tag.
  for t=1:numel(tags)
    if(not(isempty(V.(tags{t}))))
      F = scatteredInterpolant(X, Y, V.(tags{t}));
      Vi.(tags{t}) = F(Xi, Yi);
    else
      Vi.(tags{t}) = [];
    end
  end
end