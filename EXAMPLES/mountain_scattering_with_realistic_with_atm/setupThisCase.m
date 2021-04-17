clear all;
close all;
clc;

system(['cp ../mountain_scattering_with_realistic/parfile_input ./']);
system(['cp ../mountain_scattering_with_realistic/source_input ./']);
if(ismac())
  system(['sed -i "" "s/MODEL                           = default/MODEL                           = external_DG/g" parfile_input']);
else
  system(['sed -i "s/MODEL                           = default/MODEL                           = external_DG/g" parfile_input']);
end
mkdir('EXTMSH');
system(['cp ../mountain_scattering_with_realistic/EXTMSH/extMesh.msh ./EXTMSH/']);

load('../mountain_scattering_with_realistic/cross_section.mat');
az = atan2(range(yq),range(xq));

testing = 0;
if(testing)
  testM=[1:12]; testH=[0:24];
%   testM=[1:12]; testH=[12];
%   testM=[6:8]; testH=[0:3:24];
  figure();
  maxv = -Inf;
  maxx = [0, 0];
  for i=1:numel(testM)
    for j=1:numel(testH)
      mm = msishwm([mean(xq), mean(yq)], [0, 16e3], 16e3/10, datetime([2000, testM(i), 15, testH(j), 0, 0],'inputformat','yyyyMMdd@HHmmss'));
      mm.wproj = cos(az)*mm.u+sin(az)*mm.v;
      DZW = gradient(mm.wproj, mm.z);
% %       plot(model.wproj, model.z, 'displayname', [num2str(testM(i)), ' ', num2str(testH(j))]); hold on;
%       plot(DZW, model.z, 'displayname', [num2str(testM(i)), ' ', num2str(testH(j))]); hold on;
%       legend();
      if(max(DZW)>maxv)
        maxv = max(DZW);
        maxx = [testM(i), testH(j)];
      end
    end
  end
  maxv
  
else
  time = datetime([2000, 7, 15, 23, 0, 0],'inputformat','yyyyMMdd@HHmmss');
  mm = msishwm([mean(xq), mean(yq)], [0, 16e3], 16e3/10, time);
  
  mm.wproj = cos(az)*mm.u+sin(az)*mm.v;
%   apo = .5*(1+erf((mm.z/1e3)/.4 - 7));
%   apo = .5*cos(pi*mm.z/2000 - pi)+.5; apo(mm.z>=2000)=1;
  apo = .5*cos(pi*(mm.z-4000)/(4000-2600))+.5; apo(mm.z<=2600)=0; apo(mm.z>=4000)=1;
%   pause
  mm.wproj = apo .* mm.wproj;

  fig = plotModel(mm, 15e3);
  
  customSaveFig(fig, ['.',filesep,'atmospheric_model'], {'jpg', 'eps'}, 9999);
  
  matrix = [mm.z, mm.d, mm.t, mm.c, mm.p, mm.h, mm.g, mm.nsqrd, mm.k, mm.mu, mm.muvol, mm.v, mm.u, mm.wproj, mm.cp, mm.cv, mm.gamma];
  
  output_file = 'atmospheric_model.dat';
  fid = fopen(output_file, 'w');
  fprintf(fid,['MSIS20 HWM14 model\n'],matrix);
  fprintf(fid,['LON = ',sprintf('%.2f', mean(xq)),', LAT = ',sprintf('%.2f', mean(yq)),', ',datestr(time,'YYYY/mm/DD @ HH:MM:SS UTC'),'\n']);
  fprintf(fid,['\n'],matrix);
  fprintf(fid,[repmat('%15.8e ',1,size(matrix,2)),'\n'],matrix');
  fclose(fid);
end