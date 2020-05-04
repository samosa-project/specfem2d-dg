clear all;
close all;
clc;

addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% without_0__with033L_1__with300L_2 = 1; OF = 'OUTPUT_FILES_4213746_test'; IT = 45000; boxabsx = 18e3; boxy = [1500, 7.5e3];

IT = 100000; boxabsx = 28e3; boxy = [1501, 15e3];
LTOPO_over_L0 = Inf; OF = 'OUTPUT_FILES_4214695_fns';
% LTOPO_over_L0 = 0.331; OF = 'OUTPUT_FILES_4214696_fns';
% LTOPO_over_L0 = 0.332; OF = 'OUTPUT_FILES_4219868_fns';
% LTOPO_over_L0 = 3.00; OF = 'OUTPUT_FILES_4214697_fns';

forceDGMesh = 0; % Force the use of the DG mesh for interpolation? Might be heavy, but is less wrong.
dx = 20; dz = dx; % If forceDGMesh=0, wanted dx and dz for interpolation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thisFolder = [regexprep(mfilename('fullpath'),mfilename,'')];
output_figures_folder = [thisFolder, filesep, 'FIGURES'];
figfieldpath = [output_figures_folder, filesep, '',regexprep(num2str(LTOPO_over_L0),'\.','p'),'_field'];
figspectpath = [output_figures_folder, filesep, '',regexprep(num2str(LTOPO_over_L0),'\.','p'),'_spect'];
figspectdiffpath = [output_figures_folder, filesep, '',regexprep(num2str(LTOPO_over_L0),'\.','p'),'_spectdiff'];
extToSave = {'eps'};

switch(LTOPO_over_L0)
  case Inf
    OFDroot = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_without_mountains';
    TIT_addendum = 'without topography';
  case 0.331
    OFDroot = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_mountains_0.33L0';
    TIT_addendum = 'with $L_0/3$ high topography';
  case 0.332
    OFDroot = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_mountains_0.33L0_lower';
    TIT_addendum = 'with $L_0/3$ low topography';
  case 3
    OFDroot = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_mountains_3.00L0';
    TIT_addendum = 'with $3L_0$ topography';
  otherwise
    error('not implemented');
end

OFD = [OFDroot, filesep, OF];

% Read dumps.
disp(['[',mfilename,'] Reading dumps (may be long).']);
[X, Y, V] = readDumpsUnique(OFD, IT, 0);

% Select box.
disp(['[',mfilename,'] Selecting values within box.']);
select = ((Y>=min(boxy)) & (Y<=max(boxy)) & (abs(X)<=boxabsx));
X = X(select);
Y = Y(select);
V.pre = V.pre(select);

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
% V.pre = V.pre - bg_pre; % Note: this probably fuccs up below ground, but we do not care.

% Interpolate for plotting.
disp(['[',mfilename,'] Interpolating.']);
[Xi, Yi, Vi] = interpDumps(X, Y, V, range(X)/dx, range(Y)/dz, forceDGMesh);

% Figure.
TIT = ['Pressure Field ',TIT_addendum];
CBYLAB_PRESPEC = ['$\Delta P$ [Pa]'];
fig_field = figure('units','normalized','outerposition',[0,0,1,0.65]);
tightAxes = tight_subplot(1, 1, [0,0], [0.18,0.12], [0.06, 0.03]);
pcolor(Xi/1e3, Yi/1e3, Vi.pre); hold on;
daspect([1,1,1]);
hcb = colorbar(); ylabel(hcb, CBYLAB_PRESPEC, 'interpreter', 'latex');
caxis([-1,1]*max(abs(Vi.pre(:))));
colormaps_fromPython('seismic', 1);
xlabel(['$x$ [km]']); ylabel(['$z$ [km]']);
title(TIT);
customSaveFig(fig_field, [figfieldpath], extToSave, 9999);

% Compute FFT.
x = unique(Xi); z = unique(Yi);
Nx = numel(x); Nz = numel(z);
dkx = 0.5/range(x); kx_2side = (-floor(Nx/2):floor(Nx/2)) * dkx; kx = (1:floor(Nx/2)) * dkx; % [1/m]
dkz = 0.5/range(z); kz_2side = (-floor(Nz/2):floor(Nz/2)) * dkz; kz = (1:floor(Nz/2)) * dkz; % [1/m]
PRE_fft2s = fftshift(fft2(Vi.pre));
% PRE_fft2s = PRE_fft2s/max(abs(PRE_fft2s(:))); % normalise
PRE_fft1s = PRE_fft2s(ceil(Nz/2+1):end, ceil(Nx/2+1):end);
save([output_figures_folder,filesep,num2str(LTOPO_over_L0),'_PRESSURE_SPECTRUM'], 'PRE_fft1s');

% Get theoretical wavelengths.
sourcefile = [OFD, filesep, 'input_source'];
f0 = readExampleFiles_extractParam(sourcefile, 'f0', 'float');
th_kxkz = f0/soundspeed;

% Plot FFT.
COL_thkxkz = [0,1,0]*0.66; LW_thkxkz = 4; dnam_thkxkz = ['$\lambda_{\mathrm{th}}=c/f_0=',sprintf('%.0f', 1/th_kxkz),'$~m'];
XLAB = '$k_x$ [1/m]';
YLAB = '$k_z$ [1/m]';
TIT = {['Pressure Field Spectrum'],[TIT_addendum]};
toPlot=log10(abs(PRE_fft1s)); CBYLAB_PRESPEC = '$\log_{10}\left(\left|\widehat{P}\left(k_x,k_z\right)\right|\right)$';
fig_spect=figure('units','normalized','outerposition',[0,0,0.6,1]);
pcolor(kx, kz, toPlot); hold on;
h_kxkzth = plot([min(kx), [1,1]*th_kxkz], [[1,1]*th_kxkz, min(kz)], 'color', COL_thkxkz, 'linewidth', LW_thkxkz, 'displayname', dnam_thkxkz);
daspect([1,1,1]);
set(gca, 'xscale', 'log', 'yscale', 'log');
hcb = colorbar(); ylabel(hcb, CBYLAB_PRESPEC, 'interpreter', 'latex');
colormaps_fromPython('inferno', 1); caxis(max(toPlot(:))+[-2,0]);
xlabel(XLAB); ylabel(YLAB); title(TIT);
legend(h_kxkzth, 'location', 'northeast');
customSaveFig(fig_spect, [figspectpath], extToSave, 9999);

% Compare.
if(not(LTOPO_over_L0==Inf))
  basefile = [output_figures_folder,filesep,'Inf_PRESSURE_SPECTRUM.mat'];
  if(exist(basefile,'file'))
    PRE_fft1s_BASE = load(basefile); PRE_fft1s_BASE = PRE_fft1s_BASE.PRE_fft1s;
    TIT = {['Pressure Field Spectrum Difference'],[TIT_addendum, ' \textit{vs.} without topogaphy']};
    fig_spectdiff=figure('units','normalized','outerposition',[0,0,0.6,1]);
    toPlot = abs(PRE_fft1s)-abs(PRE_fft1s_BASE); % signed difference
    sgn = sign(toPlot);
    toPlot = log10(toPlot./sgn).*sgn;
    CBYLAB = ['signed difference of ',CBYLAB_PRESPEC];
    pcolor(kx, kz, toPlot); hold on;
    set(gca, 'xscale', 'log', 'yscale', 'log');
    daspect([1,1,1]);
    hcb = colorbar(); ylabel(hcb, CBYLAB, 'interpreter', 'latex');
    % adjust colormap
    maxneg = abs(min(toPlot(:))); maxpos = max(toPlot(:)); thresh = max([maxneg, maxpos])-1;
    colormaps_custom([-maxneg, [-1,-0.8,0,0.8,1]*thresh, maxpos], [[0,0,1].*[0.25,1]';[0.8,0.8,1];[1,1,1];[1,0.8,0.8];[1,0,0].*[1,0.25]'], 1);
    caxis([-maxneg,maxpos]);
    xlabel(XLAB); ylabel(YLAB); title(TIT);
    customSaveFig(fig_spectdiff, [figspectdiffpath], extToSave, 9999);
  else
    error('base file does not exist');
  end
end