clear all;
close all;
clc;

OF_rangeDep = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/LNS_GWB/OUTPUT_FILES_4236576';
OF_rangeIndep = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/LNS_GWB_rangeIndep/OUTPUT_FILES_4236578';
OFs = {OF_rangeDep, OF_rangeIndep};
% IT = 18000;
IT = 25000;
NX = 100*5;
NZ = 500*5;
maxpercenterr = 10;

% Start.
DT = 0.005;
C = 371;
D_wavefront = C*IT*DT;
XLAB = '$x$ [km]';
YLAB = '$z$ [km]';

% Read.
Vi = {};
for io = 1:numel(OFs)
  OF = OFs{io};
  [X, Y, V] = readDumpsUnique(OF, IT, 0);
  [Xi, Yi, Vi{io}] = interpDumps(X, Y, V, NX, NZ, 0);
end

% Plot.
for io = 1:numel(OFs)
  if(io==1)
    maxabsval = max(abs(Vi{io}.pre(Yi>=D_wavefront-3e3 & Yi<=(D_wavefront+1e3)))); % grab amplitude of wavefront
  end

  figure();
  pcolor(Xi/1e3, Yi/1e3, Vi{io}.pre); hold on;
  shading interp; daspect([1,1,1]);
  colormaps_fromPython('seismic', 1);
  colorbar;
  caxis([-1,1]*maxabsval);
  xlabel(XLAB); ylabel(XLAB); title(OF);
end

difference_percent = 100*((Vi{1}.pre-Vi{2}.pre)./maxabsval); % reference to rangedep (io=1)
figure();
pcolor(Xi/1e3, Yi/1e3, difference_percent); hold on;
shading interp; daspect([1,1,1]);
colormaps_fromPython('bwr', 1);
colorbar;
caxis([-1,1]*maxpercenterr);
xlabel(XLAB); ylabel(XLAB); title('difference [\%]');