do_load = 1; if(do_load); clear all; do_load = 1; end
close all;
clc;

dosave = 0;

addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils'));

IT = 120000;
OFD = ['/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_with_realistic/OUTPUT_FILES_411107_lns'];

sel_boxy = [0,12.5]*1e3;
sel_boxabsx = 50e3;
interp_dx = 15;
interp_dz = interp_dx;
fft_pcolor1_contour0 = 0;
radpatxlim = [-1,1]*30; radpatylim = [0.1, 2];
fft_selFreBand = [0.5, 2]; fft_selFreBandN = 25;
fft_selAngBand = (radpatxlim+[-1,1]*1)*pi/180; fft_selAngBandN = 150;
do_analytic_wavefront = 0;
CMAP_fft = colormaps_fromPython('bone', 0);
CMAP_fft = flipud(CMAP_fft);

fft_colourISBand=[0,1,0]*0.85;
thrsh_pfield = 0.15; blk_pfield=0.97;
COL_h_wave = [0,1,0]*0.66;
pfield_clim = [-1,1]*400;
XLIM = [-1,1]*sel_boxabsx;
YLIM = [sel_boxy];
XLAB_fft = 'horizontal wavelength $(2\pi/k_x)$ [km]';
YLAB_fft = {['vertical'], ['wavelength'],['$(2\pi/k_z)$ [km]']};
fft_clim = [-1,0];
XLIMfft = [100, 100e3]; YLIMfft = [80, 1e3];
fft_xtick = [0.1,1,10]; fft_ytick = [0.1, 1];
fftxtl = {'0.1', '1', '10'}; fftytl = {'0.1', '1'};

parfile = [OFD, filesep, 'input_parfile'];
dt = readExampleFiles_extractParam(parfile, 'DT', 'float');
if(do_load)
  wavefield01 = [OFD,filesep,'wavefield',sprintf('%07d',IT),'_01_000.txt'];
  matfile = [OFD,filesep,'wavefield',sprintf('%09d',IT),'.mat'];
  chkDateMat = dir(matfile);
  chkDateWVF = dir(wavefield01);
  if(exist(matfile,'file') & datenum(chkDateMat.date)>datenum(chkDateWVF.date))
    disp(['[',mfilename,'] Matfile exists and is more recent than .txt file, loading it instead']);
    load(matfile, 'X', 'Y', 'V');
  else
    disp(['[',mfilename,'] No matfile found, or is older than .txt file, loading .txt and saving to mat']);
    % Read dumps.
    disp(['[',mfilename,'] Reading dumps (may be long).']);
    [X, Y, V] = readDumpsUnique(OFD, IT, 0);
    save(matfile, 'X', 'Y', 'V'); % save to mat for later
  end
  % Select box.
  disp(['[',mfilename,'] Selecting values within box.']);
  select = ((Y>=min(sel_boxy)) & (Y<=max(sel_boxy)) & (abs(X)<=sel_boxabsx));
  X = X(select);
  Y = Y(select);
  V.pre = V.pre(select);
  % Interpolate for plotting.
  disp(['[',mfilename,'] Interpolating.']);
  [Xi, Yi, Vi] = interpDumps(X, Y, V, range(X)/interp_dx, range(Y)/interp_dz, 0);
end

% CBYLAB_PRESPEC = ['field of pressure perturbations $p''$ at $t=',num2str(IT*dt),'$~s [Pa]'];
CBYLAB_PRESPEC = ['$p''\left(t=',num2str(IT*dt),'\right)$ [Pa]'];

% Get speed of sound.
cp = readExampleFiles_extractParam(parfile, 'constant_p', 'float'); cv = readExampleFiles_extractParam(parfile, 'constant_v', 'float'); gamma = cp/cv;
grav = readExampleFiles_extractParam(parfile, 'gravity', 'float');
H = readExampleFiles_extractParam(parfile, 'SCALE_HEIGHT', 'float');
soundspeed = sqrt(gamma*grav*H);

% Load topo cross-section.
CS = load('cross_section.mat');
[topof, topos, ~, ~] = fft_wrap(CS.rq, CS.elvq/1e3);

if(do_analytic_wavefront)
  % Get analytic wavefront.
  z0 = 5e3;
  rho = 2730.14; vp = 6078.06; vs = 3537.23;
  [E_nu] = conversion('r', rho, 'p', vp, 's', vs, 'E', 'n'); nu = E_nu(2);
  vr = vs / ((1+nu)/(0.862+1.14*nu));
  dtair = IT*dt - z0/vp;
  x_wavefront=(CS.rq-mean(CS.rq));
  h_wavefront = CS.elvq + soundspeed*(dtair-abs(x_wavefront)/vr);
  offset = 1500;
end

% FFT.
selz = (Yi(:,1)>max(CS.elvq));
[kx, kz, PRE_fft1s, kx2s, kz2s, PRE_fft2s] = fft2_wrap(Xi(1,:), Yi(selz,1), Vi.pre(selz, :));

PRE_fft2s = PRE_fft2s(kz2s>0, kx2s~=0);
kx2s = kx2s(kx2s~=0);
kz2s = kz2s(kz2s>0); % select only upward propagating
[KX, KZ] = meshgrid(kx2s, kz2s);
[KXf, KZf, PRE_fft2s_folded] = fftfoldlr(KX, KZ, PRE_fft2s);
[LaunchAngles, Frequencies] = KxKz2LaFr(1./KX, 1./KZ, soundspeed);
[V_probed_ang, V_probed_fre, ~, V_probed_band] = LaFr2RadPattern(LaunchAngles, Frequencies, abs(PRE_fft2s), min(fft_selFreBand), max(fft_selFreBand), fft_selFreBandN, min(fft_selAngBand), max(fft_selAngBand), fft_selAngBandN);

% Plot.
fig_summary = figure('units','normalized','outerposition',[0,0,1,1]);
toppanel_z = 0.3; marg_left = 0.1; marg_right = 0.035; leftpanel_h=0.3; gaph = 0.12; gapz=0.13; gapzradpat = 0.01; margbot=0.1;
tightAxesPField = tight_subplot(1, 1, [0, 0], [1-toppanel_z, 0.018], [marg_left, marg_right]);
tightAxesSpectraDefault = tight_subplot(2, 1, [gapzradpat, 0], [margbot, toppanel_z+gapz], [marg_left, 1-leftpanel_h]);
tightAxesRadPat = tight_subplot(2, 1, [gapzradpat, 0], [margbot, toppanel_z+gapz], [leftpanel_h+gaph, marg_right]);

axes(tightAxesPField);
pcolor(Xi/1e3, Yi/1e3, Vi.pre); hold on
plot(polyshape(([CS.rq, max(CS.rq), min(CS.rq)]-mean(CS.rq))/1e3,[CS.elvq,0*[1,1]]/1e3), 'facealpha', 1, 'linestyle', 'none'); % add silhouette on top
if(do_analytic_wavefront); plot(x_wavefront/1e3, (h_wavefront-offset)/1e3, 'color', COL_h_wave, 'linestyle', ':', 'linewidth', 4); end
ylabel(['$z$ [km]']); xlabel(['$\leftarrow$ towards South $|$ $x$ [km] $|$ towards North $\rightarrow$']);
hcb1 = colorbar(); ylabel(hcb1, CBYLAB_PRESPEC, 'interpreter', 'latex');
CMAP_field = colormaps_custom([-1,-blk_pfield, -(1+thrsh_pfield)/2, -thrsh_pfield, 0, thrsh_pfield, (1+thrsh_pfield)/2, blk_pfield,1], [[0,0,1].*[0.1,0.75,1]';[0.9,0.9,1];[1,1,1];[1,0.9,0.9];[1,0,0].*[1,0.75,0.1]'], 0);
set(tightAxesPField, 'CLim', pfield_clim, 'Colormap', CMAP_field, 'XLim', XLIM/1e3, 'YLim', YLIM/1e3, 'ytick',0:2:20, 'xtick', -60:5:60);
set(tightAxesPField, 'dataaspectratio', [1,0.7,1]);

axes(tightAxesSpectraDefault(1));
pow = -1; fac = 1e-3;
absfftname = ['\left|\widehat{p''}\right|'];
CBYLAB_fft = ['radiated infrasound, $\log_{10}\left(',absfftname,'\right)$ [Pa]'];
% toPlot = log10(abs(PRE_fft1s));
% toPlot = toPlot(kz>0, kx>0);
% contourf(fac*kx(kx>0).^pow, fac*kz(kz>0).^pow, toPlot, [min(fft_clim):(range(fft_clim)/25):max(fft_clim)], 'edgecolor', 'none'); hold on;
contourf(fac*KXf.^pow, fac*KZf.^pow, log10(abs(PRE_fft2s_folded)), [min(fft_clim):(range(fft_clim)/25):max(fft_clim)], 'edgecolor', 'none'); hold on;
ylabel(YLAB_fft);
axes(tightAxesSpectraDefault(2));
loglog(fac*topof.^pow, abs(topos));
set(tightAxesSpectraDefault, 'xscale', 'log', 'yscale', 'log', 'xlim', XLIMfft*fac, 'xtick', fft_xtick);
set(tightAxesSpectraDefault(1), 'ylim', YLIMfft*fac, 'ytick', fft_ytick, 'yticklabel', fftytl, 'xticklabel', {});
set(tightAxesSpectraDefault(2), 'xticklabel', fftxtl, 'ylim', [0.1, 1e3]*fac, 'ytick', logspace(-1,3,5)*fac, 'yticklabel', {'$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '0.1', '1'});
xlabel(XLAB_fft); ylabel({'topography', 'spectrum', '$\widehat{z_\mathrm{int}}$ [km]'});

axes(tightAxesRadPat(1));
contourf(LaunchAngles*180/pi, Frequencies, log10(abs(PRE_fft2s)), [min(fft_clim):(range(fft_clim)/25):max(fft_clim)], 'edgecolor', 'none'); hold on;
for j=1:numel(fft_selFreBand); plot(radpatxlim, [1,1]*fft_selFreBand(j), 'color', fft_colourISBand); end
ylabel(['frequency [Hz]']);
hcb = colorbar(); ylabel(hcb, CBYLAB_fft, 'interpreter', 'latex');
axes(tightAxesRadPat(2));
semilogy(V_probed_ang*180/pi, V_probed_band, 'color', fft_colourISBand); hold on;
ylabel({['radiated'],'infrasound', ['$',absfftname,'$ [Pa]']});
xlabel(['$\leftarrow$ towards South $|$ launch angle [deg] $|$ towards North $\rightarrow$']);

set(tightAxesRadPat, 'yscale', 'log', 'xlim', radpatxlim);
set(tightAxesRadPat(1), 'xticklabel', {}, 'ylim', radpatylim);
set(tightAxesRadPat(2), 'ylim', 10.^fft_clim);
linkaxes(tightAxesRadPat, 'x');
set([tightAxesSpectraDefault(1); tightAxesRadPat(1)], 'colormap', CMAP_fft, 'clim', fft_clim);

hcb.Position([1, 2, 4]) = [hcb1.Position(1), tightAxesRadPat(2).Position(2), tightAxesRadPat(2).Position(4)+gapzradpat+tightAxesRadPat(1).Position(4)];
for i=1:numel(tightAxesRadPat); tightAxesRadPat(i).Position(3) = hcb1.Position(1)-0.015-tightAxesRadPat(i).Position(1); end
tightAxesRadPat(1).Position([2,4]) = tightAxesRadPat(1).Position([2,4])+[-1,1]*0.05;
tightAxesRadPat(2).Position(4) = tightAxesRadPat(2).Position(4)-0.075;
tightAxesSpectraDefault(1).Position([2,4]) = tightAxesRadPat(1).Position([2,4]);
tightAxesSpectraDefault(2).Position(4) = tightAxesRadPat(2).Position(4);

ll = add_labels_subplots(fig_summary, 0.9, 2);
set(ll(1), 'position', get(ll(1), 'position')+[0.03,0,0,0]);

set(findall(gcf,'type','text'), 'fontsize', 24);

% Save
extToSave = {'eps'};
fig_summary_path = [OFD, filesep, 'summary_realistic'];
if(dosave); customSaveFig(fig_summary, [fig_summary_path], extToSave, 9999);end;

% % move produced figures to thesis
% disp(['[] Starting to move Figures to thesis folder.']);
% system(['cp ', fig_summary_path, '.* /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/im_7_topo/realistic/']);
