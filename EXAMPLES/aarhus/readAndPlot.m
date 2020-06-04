clear all;
close all;
clc;

addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools'));

thisFolder = [regexprep(mfilename('fullpath'),mfilename,'')];

OFD = [thisFolder,filesep,'OUTPUT_FILES_401255_2090'];
% ITs = 34000; do_save = 0;
ITs = [5, [134:2:160]*1e3]; do_save = 1;

boxx = [-1,4];
boxy = [-1,1]*0.5;
dx = 35e-3/5; dz = dx; forceDGMesh = 0;
pre_t='m'; fac_t=1e3;
pre_p=['$\mu$']; fac_p=1e6;
cax = [-1,1]*3*1e-6;

parfile = [OFD,filesep,'input_parfile'];
dt = readExampleFiles_extractParam(parfile, 'DT', 'float');

for IT = ITs
  time = IT*dt;
  
  [fig_field] = plotOneIteration(OFD, IT, boxx, boxy, dx, dz, forceDGMesh, time, pre_t, fac_t, pre_p, fac_p, cax);

  name_fig = ['field_',sprintf('%.09d', IT)];

  if(do_save)
    if(numel(ITs)==1)
      customSaveFig(fig_field, [OFD, filesep, name_fig], {'fig', 'eps', 'png', 'tex'}, 9999);
    else
      customSaveFig(fig_field, [OFD, filesep, name_fig], {'jpg'}, 9999);
      close(fig_field);
      pause(1);
    end
  end
end

if(numel(ITs)>1)
  % make movie
  command = ['convert -delay 25 -resize 1000x1000 -loop 1 ',OFD,filesep,'field_*.jpg ',OFD,filesep,'field_movie.gif'];
  system(command);
end