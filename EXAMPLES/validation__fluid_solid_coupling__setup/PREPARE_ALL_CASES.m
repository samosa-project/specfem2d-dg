clear all;
% close all;
clc;

addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new'));

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
      endpoints_z = zmin*0.5 + [1, -1]*0.5*range(endpoints_x)*tan(ic_rad);
      if(range(endpoints_z)>abs(zmin))
        error(['[',mfilename,', ERROR] Cannot fit a slanted wavefront over a ',num2str(range(endpoints_z)),' m vertical span in this mesh.']);
      end
      if(range(endpoints_z)<abs(zmin)-30)
        disp(['[',mfilename,'] Source points span ',sprintf('%.0f', range(endpoints_z)),' m, while mesh is set up to span ',sprintf('%.0f', abs(zmin)),' m. It is advised to reduce zmin to -',sprintf('%.0f', range(endpoints_z)+30),' m for performance.']);
      end
    else
      endpoints_z = zmin*[1, 1]*0.5;
    end
  else
    % FTS
    if(cases{i}.ortho0_slant1)
      endpoints_z = zmax*0.5 + [-1, 1]*0.5*range(endpoints_x)*tan(ic_rad);
      if(range(endpoints_z)>abs(zmax))
        error(['[',mfilename,', ERROR] Cannot fit a slanted wavefront over a ',num2str(range(endpoints_z)),' m vertical span in this mesh.']);
      end
      if(range(endpoints_z)<zmax-30)
        disp(['[',mfilename,'] Source points span ',sprintf('%.0f', range(endpoints_z)),' m, while mesh is set up to span ',sprintf('%.0f', zmax),' m. It is advised to reduce zmax to  ',sprintf('%.0f', range(endpoints_z)+30),' m for performance.']);
      end
    else
      endpoints_z = zmax*[1, 1]*0.5;
    end
  end
%   range(endpoints_z)
  N = writeContinuousCGSource(soufile, xmin, xmax, zmin, zmax, f0, endpoints_x, endpoints_z);
  if(readExampleFiles_extractParam(parfile, 'NSOURCES', 'int')~=N)
    error(['[',mfilename,', ERROR] Number of sources in parfile does not match produced source file. Please set NSOURCES = ',num2str(N),' in the parfile.']);
  end
end

% Check predicted angles.
[rho__1, alpha__1, rho__2, alpha__2, beta__2] = get_models(parfile);
[i1_i2_j2t_j2r] = get_predicted_angles_deg(ic_rad, alpha__1, alpha__2, beta__2);
disp(['[',mfilename,'] Predicted angle i1 = ',sprintf('%.2f', i1_i2_j2t_j2r(1)),'° (STF P2P)']);
disp(['[',mfilename,'] Predicted angle i2 = ',sprintf('%.2f', i1_i2_j2t_j2r(2)),'° (FTS P2P)']);
disp(['[',mfilename,'] Predicted angle j2 = ',sprintf('%.2f', i1_i2_j2t_j2r(3)),'° (FTS P2S) or ',sprintf('%.2f', i1_i2_j2t_j2r(4)),'° (SRS P2S)']);
 % [i1, i2, j2]
 
 % Get critical angles for FTS transmission.
 icPS_deg = asin(alpha__1./[alpha__2, beta__2])*180/pi;
 disp(['[',mfilename,'] Critical angle for P-waves = ',sprintf('%.2f', icPS_deg(1)),'°. In the ray limit, no P-wave may be created for incoming waves above this angle.']);
 disp(['[',mfilename,'] Critical angle for S-waves = ',sprintf('%.2f', icPS_deg(2)),'°. In the ray limit, no S-wave may be created for incoming waves above this angle.']);
 
% Check predicted angles are below their critical angles for FTS transmission.
ic_deg = ic_rad*180/pi;
if(abs(ic_deg-icPS_deg(1))<1e-6)
  disp(['[',mfilename,', WARNING] Chosen incidence angle (',sprintf('%.2f', ic_deg),'°) is close or above critical angle for FTS-transmitted P-waves (',sprintf('%.2f', icPS_deg(1)),'°)!']);
  disp([' ',blanks(numel(mfilename)),'           Decrease v_p in solid, increase cs in fluid, or decrease incidence angle.']);
end
if(abs(ic_deg-icPS_deg(2))<1e-6)
  disp(['[',mfilename,', WARNING] Chosen incidence angle (',sprintf('%.2f', ic_deg),'°) is close or above critical angle for FTS-transmitted S-waves (',sprintf('%.2f', icPS_deg(2)),'°)!']);
  disp([' ',blanks(numel(mfilename)),'           Decrease v_s in solid, increase cs in fluid, or decrease incidence angle.']);
end





