clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% with1_or_without0 = 1; OF = 'OUTPUT_FILES_4213746_test'; IT = 45000; boxabsx = 18e3; boxy = [1500, 7.5e3];
with1_or_without0 = 1; OF = 'OUTPUT_FILES_4213746_test'; IT = 45000; boxabsx = 28e3; boxy = [1501, 25e3];

forceDGMesh = 0; % Force the use of the DG mesh for interpolation? Might be heavy, but is less wrong.
dx = 10; dz = dx; % If forceDGMesh=0, wanted dx and dz for interpolation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(with1_or_without0)
  OFDroot = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_mountains';
  TIT_addendum = 'with topography';
else
  OFDroot = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_without_mountains';
  TIT_addendum = 'without topography';
end

OFD = [OFDroot, filesep, OF];

% Read dumps.
disp(['[',mfilename,'] Reading dumps.']);
[X, Y, V] = readDumpsUnique(OFD, IT, 0);

% Select box.
disp(['[',mfilename,'] Selecting values within box.']);
select = ((Y>=min(boxy)) & (Y<=max(boxy)) & (abs(X)<=boxabsx));
X=X(select);
Y=Y(select);
V.pre=V.pre(select);

% Remove background atmosphere.
disp(['[',mfilename,'] Removing background atmophere.']);
parfile = [OFD, filesep, 'input_parfile'];
MODEL = readExampleFiles_extractParam(parfile, 'MODEL', 'string');
switch(MODEL)
  case {'default'}
    mu = readExampleFiles_extractParam(parfile, 'dynamic_viscosity', 'float');
%     kap = readExampleFiles_extractParam(parfile, 'thermal_conductivity', 'float');
    cp = readExampleFiles_extractParam(parfile, 'constant_p', 'float');
    cv = readExampleFiles_extractParam(parfile, 'constant_v', 'float');
    gamma = cp/cv;
    USE_ISOTHERMAL_MODEL = readExampleFiles_extractParam(parfile, 'USE_ISOTHERMAL_MODEL', 'bool');
    if(USE_ISOTHERMAL_MODEL)
      grav = readExampleFiles_extractParam(parfile, 'gravity', 'float');
      rho0 = readExampleFiles_extractParam(parfile, 'surface_density', 'float');
      H = readExampleFiles_extractParam(parfile, 'SCALE_HEIGHT', 'float');
      bg_rho = rho0 * exp(-Y/H);
      bg_pre = bg_rho * grav * H;
      soundspeed = sqrt(gamma*grav*H);
    else
      error('not implemented');
    end
  otherwise
    error('not implemented');
end
V.pre = V.pre - bg_pre; % Note: this probably fuccs up below ground, but we do not care.

% Interpolate for plotting.
disp(['[',mfilename,'] Interpolating.']);
[Xi, Yi, Vi] = interpDumps(X, Y, V, range(X)/dx, range(Y)/dz, forceDGMesh);

% Figure.
figure('units','normalized','outerposition',[0,0,1,1]);
pcolor(Xi/1e3, Yi/1e3, Vi.pre);
colorbar;
daspect([1,1,1]);
caxis([-1,1]*20);
colormaps_fromPython('seismic', 1);
xlabel(['$x$ [km]']);
ylabel(['$z$ [km]']);

% Compute FFT.
x = unique(Xi); z = unique(Yi);
Nx = numel(x); Nz = numel(z);
dkx = 0.5/range(x); kx_2side = (-floor(Nx/2):floor(Nx/2)) * dkx; kx = (1:floor(Nx/2)) * dkx; % [1/m]
dkz = 0.5/range(z); kz_2side = (-floor(Nz/2):floor(Nz/2)) * dkz; kz = (1:floor(Nz/2)) * dkz; % [1/m]
PRE_fft2s = fftshift(fft2(Vi.pre));
PRE_fft2s = PRE_fft2s/max(abs(PRE_fft2s(:))); % normalise
PRE_fft1s = PRE_fft2s(ceil(Nz/2+1):end, ceil(Nx/2+1):end);

% Get theoretical wavelengths.
sourcefile = [OFD, filesep, 'input_source'];
f0 = readExampleFiles_extractParam(sourcefile, 'f0', 'float');
th_kxkz = f0/soundspeed;

% Plot FFT.
COL_thkxkz = [0,1,0]*0.66; LW_thkxkz = 4; dnam_thkxkz = ['$\lambda_{\mathrm{th}}=c/f_0=',sprintf('%.0f', 1/th_kxkz),'$~m'];
XLAB = '$k_x$ [1/m]';
YLAB = '$k_z$ [1/m]';
TIT = ['Pressure Field Spectrum ',TIT_addendum];
toPlot=log10(abs(PRE_fft1s)); CBYLAB = '$\log_{10}\left(\left|\widehat{P}_{\mathrm{norm}}\left(k_x,k_z\right)\right|\right)$';
figure('units','normalized','outerposition',[0,0,0.6,1]);
pcolor(kx, kz, toPlot); hold on;

h_kxkzth = plot([min(kx), [1,1]*th_kxkz], [[1,1]*th_kxkz, min(kz)], 'color', COL_thkxkz, 'linewidth', LW_thkxkz, 'displayname', dnam_thkxkz);

daspect([1,1,1]);
set(gca, 'xscale', 'log', 'yscale', 'log');
hcb = colorbar(); ylabel(hcb, CBYLAB, 'interpreter', 'latex');
colormaps_fromPython('inferno', 1); caxis([-3,0]);
xlabel(XLAB);
ylabel(YLAB);
title(TIT);
legend(h_kxkzth, 'location', 'northeast');