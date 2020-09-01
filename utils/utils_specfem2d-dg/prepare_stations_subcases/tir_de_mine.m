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

function [xminmax, zminmax, interface, Xsource, debfin, d, name] = tir_de_mine(simulationfolder)
  addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools');
%   Xsource = [0, -25];
%   xminmax = [-440, 240];
%   zminmax = [-190, 150];
  parfile        = [simulationfolder, 'parfile_input'];
  sourcefile     = [simulationfolder, 'source_input'];
%   interfacesfile = [simulationfolder, 'interfaces_input'];
  interfacesfile = [];
  Nodes_extMesh = [simulationfolder, 'EXTMSH/Nodes_extMesh'];
  [xminmax, zminmax, Xsource] = readExampleFiles(parfile, sourcefile, interfacesfile, Nodes_extMesh);
  
  ISAEIS2_z = 65; % test report
  ISAEIS3_z = 104; % test report
  
  IMU_S5_z = ISAEIS3_z + 1; % test report page 15 (positions of IMUs on balloon) + page 11 (get altitude of bottom of balloon)
  IMU_S4S3_z = ISAEIS3_z + 1 + 3.10/2; % test report page 15 (positions of IMUs on balloon) + page 11 (get altitude of bottom of balloon) + page 6 (use balloon specs to get nose/tail altitude)
  IMU_S4S3_dx = 7/2; % put worst case scenario for dx (or dr): ballon is aligned with direction vector to shot, and thus S4/S3 are furthest/closest or close/furthest of shot, use page 6 (use balloon specs to get nose/tail distance)
  
  to_gla08_z = (195.-167.); % test report
  to_gla08_r = 187.; % test report
  to_gla08_x = (to_gla08_r^2.-to_gla08_z^2.)^0.5; % deduce
  xr0 = -to_gla08_x/to_gla08_r; zr0 = -to_gla08_z/to_gla08_r; zx0 = -to_gla08_z/to_gla08_x;
  
  groundclearance = 1.5;
  isitok=-1;
  while(not(ismember(isitok,[0,1])))
    isitok=input(['[', mfilename, '] Ground clearance is planned to be ',num2str(groundclearance),' [m]. Is that ok (0 for no, 1 for yes)? > ']);
  end
  
  % interface for plots and checks
%   interface = [-400, 10;400*zx0, -10*zx0];
  interface = [-450,      -400,      10,      10,         25,         25,         250        ; ...
                -60.5759,  -60.5759,  1.5143,  1.5143-15,  1.5143-15,  1.5143-30,   1.5143-30];
  isitok=-1;
  while(not(ismember(isitok,[0,1])))
    disp(['[', mfilename, '] Interface for plots and checks is planned to be [ ',sprintf('%7.2f ',interface(1,:))]);
    isitok=input([blanks(numel(mfilename)),'                                                     ',sprintf('%7.2f ',interface(2,:)),']. Check with ''cd EXTMSH/; gmsh extMesh.geo; cd ..'' Is that ok (0 for no, 1 for yes)? > ']);
  end

  lid = 0;
  lid = lid+1; d(lid) = 0; r = 0;   name{lid} = ['above source (@r = 0)', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 100; name{lid} = ['GLA04 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 193; name{lid} = ['GLA05 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 150; name{lid} = ['GLA06 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 5;   name{lid} = ['GLA07 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['GLA08 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0; idGLA08 = lid;
  lid = lid+1; d(lid) = 0; r = 185; name{lid} = ['GLA09 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  lid = lid+1; d(lid) = 0; r = 50;  name{lid} = ['GLA11 @r = ', num2str(r)]; debfin(lid, 1, :) = [1, 1]*r*xr0; debfin(lid, 2, :) = [1, 1]*r*zr0;
  debfin(:, 2, 1) = debfin(:, 2, 1)+groundclearance; % beginning above ground
  debfin(:, 2, 2) = debfin(:, 2, 2)-groundclearance; % end under ground
  d = d+2*groundclearance; % separate correctly the ground stations
  
  lid = lid+1; d(lid) = 1; r = 0; name{lid} = ['Monitor Source']; debfin(lid, 1, :) = Xsource(1)*[1, 1]; debfin(lid, 2, :) = (Xsource(2)+d(lid))*[1, 1];
  
  lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['ISAEIS2, above GLA08 (@r = ', num2str(r), '), z = +', num2str(ISAEIS2_z), ' (ISAEIS ground is GLA08 z = +1m)']; debfin(lid, :, :) = debfin(idGLA08, :, :); debfin(lid, 2, :) = debfin(lid, 2, :)+ISAEIS2_z;
  lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['ISAEIS3, above GLA08 (@r = ', num2str(r), '), z = +', num2str(ISAEIS3_z), ' (ISAEIS ground is GLA08 z = +1m)']; debfin(lid, :, :) = debfin(idGLA08, :, :); debfin(lid, 2, :) = debfin(lid, 2, :)+ISAEIS3_z;
  
%   lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['IMU S5 on bottom of balloon (@r = ', num2str(r), '), z = +', num2str(IMU_S5_z), ' (which is bottom of balloon altitude, or ISAEIS3 + 1m)']; debfin(lid, :, :) = debfin(idGLA08, :, :); debfin(lid, 2, :) = debfin(lid, 2, :)+IMU_S5_z;
%   lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['IMU S3 and/or S4 on nose/tail of balloon (@r = ', num2str(r), '-',num2str(IMU_S4S3_dx),' whichever is furthest to shot), z = +', num2str(IMU_S4S3_z), ' (which is nose/tail of balloon altitude)']; debfin(lid, 1, :) = debfin(idGLA08, 1, :)-IMU_S4S3_dx; debfin(lid, 2, :) = debfin(lid, 2, :)+IMU_S4S3_z;
%   lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['IMU S3 and/or S4 on nose/tail of balloon (@r = ', num2str(r), '+',num2str(IMU_S4S3_dx),' whichever is closest to shot), z = +', num2str(IMU_S4S3_z), ' (which is nose/tail of balloon altitude)']; debfin(lid, 1, :) = debfin(idGLA08, 1, :)+IMU_S4S3_dx; debfin(lid, 2, :) = debfin(lid, 2, :)+IMU_S4S3_z;
  
  lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['IMU S5 on bottom of balloon (@X = X_ISAEIS3 + [0, ',num2str(IMU_S5_z-ISAEIS3_z),'])']; debfin(lid, :, :) = debfin(idGLA08, :, :); debfin(lid, 2, :) = debfin(lid, 2, :)+IMU_S5_z;
  lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['IMU S3 and/or S4 on nose/tail (whichever is furthest to shot) of balloon (@X = X_ISAEIS3 + [-',num2str(IMU_S4S3_dx),', ',num2str(IMU_S4S3_z-ISAEIS3_z),'])']; debfin(lid, 1, :) = debfin(idGLA08, 1, :)-IMU_S4S3_dx; debfin(lid, 2, :) = debfin(idGLA08, 2, :)+IMU_S4S3_z;
  lid = lid+1; d(lid) = 0; r = 187; name{lid} = ['IMU S3 and/or S4 on nose/tail (whichever is closest  to shot) of balloon (@X = X_ISAEIS3 + [+',num2str(IMU_S4S3_dx),', ',num2str(IMU_S4S3_z-ISAEIS3_z),'])']; debfin(lid, 1, :) = debfin(idGLA08, 1, :)+IMU_S4S3_dx; debfin(lid, 2, :) = debfin(idGLA08, 2, :)+IMU_S4S3_z;
end