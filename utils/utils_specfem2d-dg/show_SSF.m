% Author:        Léo Martire.
% Description:   Plots a spread Source Spatial Function (SSF), as output by SPECFEM2D-DG if the flag SPREAD_SSF_SAVE was set to .true. in the parfile.
% Notes:         N. A.
%
% Usage:
%   [] = show_EBF()
% with:
%   N. A.
% yields:
%   N. A.

function show_SSF()
  % User inputs.
  disp(['[',mfilename,'] User inputs.']);
  folder = input([blanks(numel(mfilename)),'   Folder containing the SSF files? > '], 's');
  if(not(strcmp(folder(end), filesep)))
    folder = [folder, filesep];
  end
  listSSF = dir(strcat(folder, 'SSF*'));
  if(numel(listSSF)==0)
    error(['[',mfilename,', ERROR] No SSF files found in this folder.']);
  end
  listSSF.name
  nsources = input([blanks(numel(mfilename)),'   Number of sources? > ']);
  gs = input([blanks(numel(mfilename)),'   Plotting grid size? > ']);
  ar = input([blanks(numel(mfilename)),'   Axes aspect ratio (format [ar_x, ar_y, ar_z])? > ']);
  th = input([blanks(numel(mfilename)),'   Threshold to recenter plot (show only where SSF>threshold)? > ']);

  % Load the data.
  disp(['Loading.']);
  for s = 1:nsources
    a = [];
    for i = 1:length(listSSF)
      if (not(isempty(regexp(listSSF(i).name, strcat(num2str(s), '_')))))
        a = unique([a; importdata(strcat(folder, listSSF(i).name))], 'rows');
        if (mod(i, floor(length(listSSF) / 10)) == 0 || i == length(listSSF))
           disp(strcat('Loading data (', num2str(100 * i / length(listSSF)), '%).'));
        end
      end
    end
    x{s} = a(:, 1); z{s} = a(:, 2); d{s} = a(:, 3);
  end

  % Interpolate the point clouds and plot.
  disp(['Interpolating.']);
  for s = 1:nsources
    figure(1000 * s);
    xx = x{s}; zz = z{s}; dd = d{s};
    F = scatteredInterpolant(x{s}, z{s}, d{s});
    tx = min(xx(:)):gs:max(xx(:));
    tz = min(zz(:)):gs:max(zz(:));
    [X, Z] = meshgrid(tx, tz);
    D = F(X, Z);
    surf(X, Z, D, 'FaceColor', 'interp', 'EdgeColor', 'white', 'LineStyle', ':');
    set(gca, 'DataAspectRatio', ar);
    xlabel('$x$'); ylabel('$z$'); title(strcat("Source Spatial Function ", num2str(s)));
    xlim([min(xx(d{s} > th)), max(xx(d{s} > th))]); ylim([min(zz(d{s} > th)), max(zz(d{s} > th))]); zlim([th, max(d{s})]);
  end
end