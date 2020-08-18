clear all;
% close all;
clc;

addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new'));
thisFolder = [regexprep(mfilename('fullpath'),mfilename,'')];
setup;

for i = 1:numel(cases)
  % Get working folders.
  folder = folderz{i};
  
  % Copy common configuration files.
  src = [thisFolder,filesep,'common',filesep,'*'];
  system(['cp ', src, ' ', folder]);
  
  % Prepare source.
  parfile = [folder, 'parfile_input'];
  soufile = [folder, 'source_input'];
  intfile = [folder, 'interfaces_input'];
  xmin = readExampleFiles_extractParam(parfile, 'xmin', 'float');
  xmax = readExampleFiles_extractParam(parfile, 'xmax', 'float');
  [~, ~, zmin, zmax] = readExampleFiles_meshfem_mesh(intfile);
  
  endpoints_x = [xmin, xmax] + [1, -1]*0.25;
  if(cases{i}.fts0_stf1)
    % STF
    if(cases{i}.ortho0_slant1)
      endpoints_z = zmin*0.5 + [1, -1]*0.5*range(endpoints_x)*tan(ic);
      if(range(endpoints_z)>abs(zmin))
        error(['[',mfilename,', ERROR] Cannot fit a slanted wavefront over a ',num2str(range(endpoints_z)),' m vertical span in this mesh.']);
      end
    else
      endpoints_z = zmin*[1, 1]*0.5;
    end
  else
    % FTS
    if(cases{i}.ortho0_slant1)
      endpoints_z = zmax*0.5 + [-1, 1]*0.5*range(endpoints_x)*tan(ic);
      if(range(endpoints_z)>abs(zmax))
        error(['[',mfilename,', ERROR] Cannot fit a slanted wavefront over a ',num2str(range(endpoints_z)),' m vertical span in this mesh.']);
      end
    else
      endpoints_z = zmax*[1, 1]*0.5;
    end
  end
%   endpoints_z
  writeContinuousCGSource(soufile, xmin, xmax, zmin, zmax, f0, endpoints_x, endpoints_z)
end

% Check predicted angles.
[i1_i2_j2t_j2r] = get_predicted_angles(parfile, ic);
disp(['[',mfilename,'] Predicted angle i1 = ',sprintf('%.2f', i1_i2_j2t_j2r(1)),'° (SP 2 FP)']);
disp(['[',mfilename,'] Predicted angle i2 = ',sprintf('%.2f', i1_i2_j2t_j2r(2)),'° (FP 2 SP)']);
disp(['[',mfilename,'] Predicted angle j2 = ',sprintf('%.2f', i1_i2_j2t_j2r(3)),'° (FP 2 SS) or ',sprintf('%.2f', i1_i2_j2t_j2r(4)),'° (SP 2 SS)']);
 % [i1, i2, j2]






