function [Xi, Yi, Vi] = interpDumps(X, Y, V, nx, ny)
  tags = {'rho', 'vel', 'pre'};
  
  x_exact = linspace(min(X), max(X), nx);
  y_exact = linspace(min(Y), max(Y), ny);
  [Xi, Yi] = meshgrid(x_exact,y_exact);

  % Interpolate dump, on each non-empty tag.
  for t=1:numel(tags)
    if(not(isempty(V.(tags{t}))))
      % For vel, make sure to do column by column.
      if(strcmp(tags{t}, 'vel'))
        ndim = size(V.(tags{t}),2);
        switch(ndim)
          case 2
            axx = 'xz';
          case 3
            axx = 'xyz';
          otherwise
            error('kek');
        end
        for dd = 1:ndim
          F = scatteredInterpolant(X, Y, V.(tags{t})(:,dd));
          Vi.([tags{t},axx(dd)]) = F(Xi, Yi);
        end
      else
        F = scatteredInterpolant(X, Y, V.(tags{t}));
        Vi.(tags{t}) = F(Xi, Yi);
      end
    else
      Vi.(tags{t}) = [];
    end
  end
end