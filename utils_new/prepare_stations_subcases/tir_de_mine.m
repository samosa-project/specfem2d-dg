% Author:        Léo Martire.
% Description:   Prepares stations for a specfic type of run.
% Notes:         TODO.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function [xminmax, zminmax, interface, Xsource, debfin, d, name] = tir_de_mine()
  addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools');
  Xsource = [0, -25];
  xminmax = [-440, 240];
  zminmax = [-190, 150];
  ISAEIS2_z = 65;
  ISAEIS3_z = 104;
  to_gla08_z = (195.-167.);
  to_gla08_r = 187.;
  to_gla08_x = (to_gla08_r^2.-to_gla08_z^2.)^0.5;
  xr0 = -to_gla08_x/to_gla08_r; zr0 = -to_gla08_z/to_gla08_r; zx0 = -to_gla08_z/to_gla08_x;
  interface = [-400, 10;400*zx0, -10*zx0];

  lid = 0;
  lid = lid+1; d(lid) = 0; r = 0;   name{lid} = ['above source (@r = 0)', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 100; name{lid} = ['GLA04 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 193; name{lid} = ['GLA05 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 150; name{lid} = ['GLA06 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 5;   name{lid} = ['GLA07 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['GLA08 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0; idGLA08 = lid;
  lid = lid+1; d(lid) = 0; r = 185; name{lid} = ['GLA09 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 50;  name{lid} = ['GLA11 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  debfin(:, 2, 1) = debfin(:, 2, 1)+1; % beginning above ground
  debfin(:, 2, 2) = debfin(:, 2, 2)-1; % end under ground
  d = d+2; % separate correctly
  lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['ISAE IS Sensor 2, above GLA08 (@r = ', num2str(r), '), z = +', num2str(ISAEIS2_z), ' (ISAEIS ground is GLA08 z = +1m)']; debfin(lid, :, :) = debfin(idGLA08, :, :); debfin(lid, 2, :) = debfin(lid, 2, :)+ISAEIS2_z;
  lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['ISAE IS Sensor 3, above GLA08 (@r = ', num2str(r), '), z = +', num2str(ISAEIS3_z), ' (ISAEIS ground is GLA08 z = +1m)']; debfin(lid, :, :) = debfin(idGLA08, :, :); debfin(lid, 2, :) = debfin(lid, 2, :)+ISAEIS3_z;
  lid = lid+1; d(lid) = 1; r = 0; name{lid} = ['Monitor Source']; debfin(lid, 1, :) = Xsource(1)*[1, 1]; debfin(lid, 2, :) = (Xsource(2)+d(lid))*[1, 1];
end