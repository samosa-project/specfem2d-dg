clear all;
close all;
clc;

with1_or_without0 = 0;

box_absxmax = [8e3];
box_yminmax = [0, 2e3];

forceDGMesh = 0; % Force the use of the DG mesh for interpolation? Might be heavy, but is less wrong.
dx = 10; dz = dx; % If forceDGMesh=0, wanted dx and dz for interpolation.

if(with1_or_without0)
  OFDroot = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_mountains';
else
  OFDroot = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_without_mountains';
end

OFD = [OFDroot, filesep, 'OUTPUT_FILES_4213746_test']; IT = 5000;

[X, Y, V] = readDumpsUnique(OFD, IT, 0);

% Select box.
select = ((Y>=min(box_yminmax)) & (Y<=max(box_yminmax)) & (abs(X)<=box_absxmax));
X=X(select);
Y=Y(select);
V.pre=V.pre(select);

% Remove background atmosphere.
parfile = [OFD, filesep, 'parfile_input'];
MODEL = readExampleFiles_extractParam(parfile, 'MODEL', 'string');
switch(MODEL)
  case {'default'}
    mu = readExampleFiles_extractParam(parfile, 'dynamic_viscosity', 'float');
%     kap = readExampleFiles_extractParam(parfile, 'thermal_conductivity', 'float');
%     cp = readExampleFiles_extractParam(parfile, 'constant_p', 'float');
%     cv = readExampleFiles_extractParam(parfile, 'constant_v', 'float');
%     gamma = cp/cv;
    USE_ISOTHERMAL_MODEL = readExampleFiles_extractParam(parfile, 'USE_ISOTHERMAL_MODEL', 'bool');
    if(USE_ISOTHERMAL_MODEL)
      grav = readExampleFiles_extractParam(parfile, 'gravity', 'float');
      rho0 = readExampleFiles_extractParam(parfile, 'surface_density', 'float');
      H = readExampleFiles_extractParam(parfile, 'SCALE_HEIGHT', 'float');
      bg_rho = rho0 * exp(-Y/H);
      bg_pre = bg_rho * grav * H;
    else
      error('not implemented');
    end
  otherwise
    error('not implemented');
end
V.pre = V.pre - bgpre; % Note: this probably fuccs up below ground, but we do not care.

% Interpolate for plotting.
[Xi, Yi, Vi] = interpDumps(X, Y, V, range(X)/dx, range(Y)/dz, forceDGMesh);

% Figure.
figure();
pcolor(Xi/1e3, Yi/1e3, Vi.pre);
colorbar;
daspect([1,1,1]);
caxis([-1,1]*20);
colormaps_fromPython('seismic', 1);
xlabel(['$x$ [km]']);
ylabel(['$z$ [km]']);