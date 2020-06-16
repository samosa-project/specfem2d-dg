do_load = 0;

if(do_load)
  clear all;
  do_load = 1;
end
close all;
clc;

addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new'));

rtdwprfx = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_';
thisFolder = [regexprep(mfilename('fullpath'),mfilename,'')];
output_figures_folder = [thisFolder, filesep, 'FIGURES'];
baseline_fft_path = [output_figures_folder,filesep,'BASELINE_PRESSURE_SPECTRUM'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% without_0__with033L_1__with300L_2 = 1; OF = 'OUTPUT_FILES_4213746_test'; IT = 45000; boxabsx = 18e3; boxy = [1500, 7.5e3];
% IT = 140000;
% LTOPO_over_L0 = Inf; OF = 'OUTPUT_FILES_4214695_lns'; LTOPO_over_L0 = 0.331; OF = 'OUTPUT_FILES_4214696_lns'; LTOPO_over_L0 = 0.332; OF = 'OUTPUT_FILES_4219868_lns'; LTOPO_over_L0 = 3.00; OF = 'OUTPUT_FILES_4214697_lns';
% LTOPO_over_L0 = Inf; OF = 'OUTPUT_FILES_398709_fns'; LTOPO_over_L0 = 0.331; OF = 'OUTPUT_FILES_4214696_lns'; LTOPO_over_L0 = 0.332; OF = 'OUTPUT_FILES_4219868_lns'; LTOPO_over_L0 = 3.00; OF = 'OUTPUT_FILES_398703_fns';

wavenumber0_wavelength1 = 1; % plot in terms of wavenumber (0) or wavelength (1)
boxabsx = 28e3; boxy = [1501, 15e3];
forceDGMesh = 0; dx = 10; dz = dx; % dx=10 ok, dx<10 chugs hard

climpfield = [-1, 1]*175; thrsh_pfield = 0.25; blk_pfield=0.98;
climfft = [4.5, 6.25];
climcomp = [-1, 1]*6; maxneg_fftcomp = abs(min(climcomp)); maxpos_fftcomp = abs(max(climcomp)); thresh_fftcomp = max([maxneg_fftcomp, maxpos_fftcomp])-0.5;
extToSave = {'eps'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figfieldpath = [output_figures_folder, filesep, '',regexprep(num2str(LTOPO_over_L0),'\.','p'),'_field'];
% figspectpath = [output_figures_folder, filesep, '',regexprep(num2str(LTOPO_over_L0),'\.','p'),'_spect'];
% figspectdiffpath = [output_figures_folder, filesep, '',regexprep(num2str(LTOPO_over_L0),'\.','p'),'_spectdiff'];

% switch(LTOPO_over_L0)
%   case Inf
%     OFDroot = [rtdwprfx,'without_mountains'];
%     TIT_addendum = 'without topography';
%   case 0.331
%     OFDroot = [rtdwprfx,'with_mountains_0.33L0'];
%     TIT_addendum = 'with $L_0/3$ high topography';
%   case 0.332
%     OFDroot = [rtdwprfx,'with_mountains_0.33L0_lower'];
%     TIT_addendum = 'with $L_0/3$ low topography';
%   case 3
%     OFDroot = [rtdwprfx,'with_mountains_3.00L0'];
%     TIT_addendum = 'with $3L_0$ topography';
%   otherwise
%     error('not implemented');
% end

% OFDroots = {[rtdwprfx,'without_mountains']};
% OFs{}

% tg_bsln = '4214695_lns'; tg_03hig = '4214696_lns'; tg_03low = '4219868_lns'; tg_3 = '4214697_lns'; IT = 100000;
% tg_bsln = '398709_fns'; tg_03hig = '4214696_lns'; tg_03low = '4219868_lns'; tg_3 = '398703_fns'; IT = 140000;
tg_bsln = '408827_lns'; tg_03hig = '4214696_lns'; tg_03low = '4219868_lns'; tg_3 = '408826_lns'; IT = 100000;

zoombox_x = [-1,1]*14e3; zoombox_z = [8,12]*1e3;

OFDs = {[rtdwprfx, 'without_mountains', filesep, 'OUTPUT_FILES_', tg_bsln], ...
        [rtdwprfx, 'with_mountains_0.33L0', filesep, 'OUTPUT_FILES_', tg_03hig], ...
        [rtdwprfx, 'with_mountains_0.33L0_lower', filesep, 'OUTPUT_FILES_', tg_03low], ...
        [rtdwprfx, 'with_mountains_3.00L0', filesep, 'OUTPUT_FILES_', tg_3]};
IDs_to_process = 1:4;
figfieldpath = [output_figures_folder, filesep, 'all_fields'];
figspectpath = [output_figures_folder, filesep, 'all_spects'];
figspectdiffpath = [output_figures_folder, filesep,'all_spectdiffs'];

do_pfield = 1;
do_fft = 1;
do_comparefft = 1;

% Load.
if(do_load)
  Xi = {};
  Yi = {};
  Vi = {};
  th_kxkz = {};
  soundspeed = {};
  dt=[];
  for i=IDs_to_process
    OFD = OFDs{i};%[OFDroots{i}, filesep, OF{i}];

    % Read dumps.
    disp(['[',mfilename,'] Reading dumps (may be long).']);
    [X, Y, V] = readDumpsUnique(OFD, IT, 0);

    % Select box.
    disp(['[',mfilename,'] Selecting values within box.']);
    select = ((Y>=min(boxy)) & (Y<=max(boxy)) & (abs(X)<=boxabsx));
    X = X(select);
    Y = Y(select);
    V.pre = V.pre(select);

    % Eventually remove background atmosphere (if FNS).
    parfile = [OFD, filesep, 'input_parfile'];
    MODEL = readExampleFiles_extractParam(parfile, 'MODEL', 'string');
    dt(i) = readExampleFiles_extractParam(parfile, 'DT', 'float');
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
          soundspeed{i} = sqrt(gamma*grav*H);
        else
          error('not implemented');
        end
      otherwise
        error('not implemented');
    end
    % Get theoretical wavelengths.
    sourcefile = [OFD, filesep, 'input_source'];
    f0 = readExampleFiles_extractParam(sourcefile, 'f0', 'float');
    th_kxkz{i} = f0/soundspeed{i};
      
    USE_LNS = readExampleFiles_extractParam(parfile, 'USE_LNS', 'bool');
    if(not(USE_LNS))
      disp(['[',mfilename,'] Removing background atmophere.']);
      V.pre = V.pre - bg_pre; % Note: this probably fuccs up below ground, but we do not care.
    end

    % Interpolate for plotting.
    disp(['[',mfilename,'] Interpolating.']);
    [Xi{i}, Yi{i}, Vi{i}] = interpDumps(X, Y, V, range(X)/dx, range(Y)/dz, forceDGMesh);
  end
end

if(not(numel(unique(dt*IT))==1))
  error('dt varies from one simulation to the other');
end

% Figure.
if(do_pfield)
%   TIT = ['Pressure Field ',TIT_addendum];
  CBYLAB_PRESPEC = ['field of pressure perturbations $p''$ at $t=',num2str(IT*dt(1)),'$~s [Pa]'];
  fig_field = figure('units','normalized','outerposition',[0,0,1,1]);
  tightAxes = tight_subplot(4, 1, [0.03,0.], [0.11,0.05], [0.06, 0.06]);
  
  for i=IDs_to_process
    axes(tightAxes(i));
    selx = (Xi{i}(1,:)>=min(zoombox_x)*1.1) & (Xi{i}(1,:)<=max(zoombox_x)*1.1);
    selz = (Yi{i}(:,1)>=min(zoombox_z)*0.9) & (Yi{i}(:,1)<=max(zoombox_z)*1.1);
    pcolor(Xi{i}(selz,selx)/1e3, Yi{i}(selz,selx)/1e3, Vi{i}.pre(selz,selx)); hold on;
    ylabel(['$z$ [km]']);
  end
  xlabel(['$x$ [km]']);
  
  hcb = colorbar(); ylabel(hcb, CBYLAB_PRESPEC, 'interpreter', 'latex');
  XTCK_pfield = (min(zoombox_x):2e3:max(zoombox_x))/1e3;
  ZTCK_pfield = (min(zoombox_z):2e3:max(zoombox_z))/1e3;
%   caxis([-1,1]*max(abs(Vi.pre(:))));
  
  CMAP = colormaps_custom([-1,-blk_pfield, -(1+thrsh_pfield)/2, -thrsh_pfield, 0, thrsh_pfield, (1+thrsh_pfield)/2, blk_pfield,1], [[0,0,1].*[0.2,0.75,1]';[0.9,0.9,1];[1,1,1];[1,0.9,0.9];[1,0,0].*[1,0.75,0.2]'], 0);
  
%   set(tightAxes, 'DataAspectRatio', [1,1,1]);
  set(tightAxes, 'CLim', climpfield, 'Colormap', CMAP, 'XLim', zoombox_x/1e3, 'YLim', zoombox_z/1e3, 'XTick', XTCK_pfield, 'YTick', ZTCK_pfield);
  for i=1:numel(tightAxes)
    tightAxes(i).Position = [tightAxes(end).Position(1), tightAxes(i).Position(2), tightAxes(4).Position(3:4)];
  end
  set(tightAxes(1:end-1), 'XTickLabel', {});
  
  height_spanning_all_axes = sum(tightAxes(1).Position([2,4]))-tightAxes(end).Position(2);
  set(hcb, 'Position', [hcb.Position(1) + 0.005, tightAxes(end).Position(2), hcb.Position(3), height_spanning_all_axes],'fontsize', 26);
  ll = add_labels_subplots(fig_field,0.9);
  
%   title(TIT);
  customSaveFig(fig_field, [figfieldpath], extToSave, 9999);
end

% Compute FFT.
PRE_fft1s = {};
if(any([do_fft, do_comparefft]))
  for i=IDs_to_process
    x = unique(Xi{i}); z = unique(Yi{i});
    Nx = numel(x); Nz = numel(z);
    dkx = 0.5/range(x); kx_2side = (-floor(Nx/2):floor(Nx/2)) * dkx; kx = (1:floor(Nx/2)) * dkx; % [1/m]
    dkz = 0.5/range(z); kz_2side = (-floor(Nz/2):floor(Nz/2)) * dkz; kz = (1:floor(Nz/2)) * dkz; % [1/m]
    PRE_fft2s = fftshift(fft2(Vi{i}.pre));
    % PRE_fft2s = PRE_fft2s/max(abs(PRE_fft2s(:))); % normalise
    PRE_fft1s{i} = PRE_fft2s(ceil(Nz/2+1):end, ceil(Nx/2+1):end);
    
    if(i==1)
      baseline_fft = PRE_fft1s{i};
%       save([output_figures_folder,filesep,num2str(LTOPO_over_L0),'_PRESSURE_SPECTRUM'], 'PRE_fft1s');
      save([baseline_fft_path], 'baseline_fft');
      disp(['saving baseline']);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots.
if(wavenumber0_wavelength1)
  fac = 1e-3; XLAB = 'horizontal wavelength $(1/k_x)$ [km]'; YLAB = 'vertical wavelength $(1/k_z)$ [km]';
  pow = -1;
  minx = 100*fac; minz = minx;
else
  XLAB = '$k_x$ [1/m]';
  YLAB = '$k_z$ [1/m]'; pow=1;
end
margz = [0.16, 0.04]; margh = [0.06, 0.085]; gap=[0, 0.01];

% Plot FFT.
if(do_fft)
  COL_thkxkz = [0,1,0]*0.66; LW_thkxkz = 4; dnam_thkxkz = ['$\lambda_{\mathrm{th}}=c/f_0=',sprintf('%.0f', 1/th_kxkz{i}),'$~m'];
%   TIT = {['Pressure Field Spectrum'],[TIT_addendum]};
  
  CBYLAB_PRESPEC = '$\log_{10}\left(\left|\widehat{P}\left(k_x,k_z\right)\right|\right)$';
  fig_spect = figure('units', 'normalized', 'outerposition', [0,0,1,0.7]);
  tightAxes_spec = tight_subplot(1, 4, gap, margz, margh);
  
  for i=IDs_to_process
%   for i=1:4
    axes(tightAxes_spec(i));
    toPlot = log10(abs(PRE_fft1s{i}));
    pcolor(fac*kx.^pow, fac*kz.^pow, toPlot); hold on;
    h_kxkzth = plot(fac*[min(kx), [1,1]*th_kxkz{i}].^pow, fac*[[1,1]*th_kxkz{i}, min(kz)].^pow, 'color', COL_thkxkz, 'linewidth', LW_thkxkz, 'displayname', dnam_thkxkz);
    if(i==2)
      xlabel(XLAB);
    end
    if(i==1)
      ylabel(YLAB);
    end
  end
  
  % daspect([1,1,1]);
%   set(gca, 'xscale', 'log', 'yscale', 'log');
  hcb1 = colorbar(); ylabel(hcb1, CBYLAB_PRESPEC, 'interpreter', 'latex');
  % colormaps_fromPython('inferno', 1);
  % caxis(max(toPlot(:))+[-2,0]);
%   caxis(climfft);
  
%   title(TIT);
  
%   CMAP = colormaps_fromPython('magma', 0, floor(0.25*512):floor(1*512));
  CMAP = colormaps_fromPython('bone', 0);
  CMAP = flipud(CMAP);
  set(tightAxes_spec, 'xscale', 'log', 'yscale', 'log', 'colormap', CMAP, 'clim', climfft);
  set(tightAxes_spec(2:end), 'YTickLabel', {});
  
  if(wavenumber0_wavelength1)
%     legend(h_kxkzth, 'location', 'northwest');
    set(tightAxes_spec, 'xlim', [minx, max(xlim)], 'ylim', [minz, max(ylim)], 'xtick', logspace(log10(minx), ceil(log10(max(xlim))),ceil(log10(max(xlim)))-log10(minx)+1));
    xtl = split(sprintf('%.2g|',xticks),'|'); %xticklabels(xtl(1:end-1));
    ytl = split(sprintf('%.2g|',yticks),'|'); %yticklabels(ytl(1:end-1));
    set(tightAxes_spec(1), 'yticklabels', ytl(1:end-1));
    set(tightAxes_spec, 'xticklabels', xtl(1:end-1));
  else
    legend(h_kxkzth, 'location', 'northeast');
    error('kek');
  end
  
  hshift_cb = 0.04;
  hcb1.Position([2,4]) = tightAxes_spec(end).Position([2,4]);
  hcb1.Position([1]) = sum(tightAxes_spec(end).Position([1,3]))+hshift_cb;
  
  ll = add_labels_subplots(fig_spect,1);
  
  customSaveFig(fig_spect, [figspectpath], extToSave, 9999);
end

% Compare.
if(do_comparefft)
%   if(not(LTOPO_over_L0==Inf))
%     basefile = [output_figures_folder,filesep,'Inf_PRESSURE_SPECTRUM.mat'];
  basefile = [baseline_fft_path,'.mat'];
  if(exist(basefile,'file'))
    PRE_fft1s_BASE = load(basefile);
%       PRE_fft1s_BASE = PRE_fft1s_BASE.PRE_fft1s;
    PRE_fft1s_BASE = PRE_fft1s_BASE.baseline_fft;
%     TIT = {['Pressure Field Spectrum Difference'],[TIT_addendum, ' \textit{vs.} without topogaphy']};
    fig_spectdiff = figure('units','normalized','outerposition',[0,0,1,0.7]);
    margh = margh + [0, -0.01];
    tightAxes_diff = tight_subplot(1, 3, gap, margz, margh);

    CBYLAB = ['signed difference of ',CBYLAB_PRESPEC];

    for i=IDs_to_process(2:end)
%       for i=2:4
      axes(tightAxes_diff(i-1));
      toPlot = abs(PRE_fft1s{i})-abs(PRE_fft1s_BASE); % signed difference
      sgn = sign(toPlot);
      toPlot = log10(toPlot./sgn).*sgn;
      pcolor(fac*kx.^pow, fac*kz.^pow, toPlot); hold on;
      if(i==3)
        xlabel(XLAB);
      end
      if(i==2)
        ylabel(YLAB);
      end
    end

%     daspect([1,1,1]);
    hcb = colorbar(); ylabel(hcb, CBYLAB, 'interpreter', 'latex');
    % adjust colormap
%     maxneg = abs(min(toPlot(:))); maxpos = max(toPlot(:));

    CMAP = colormaps_custom([-maxneg_fftcomp, [-1,-0.85,0,0.85,1]*thresh_fftcomp, maxpos_fftcomp], [[0,0,1].*[0.25,1]';[0.9,0.9,1];[1,1,1];[1,0.9,0.9];[1,0,0].*[1,0.25]'], 0);
%     caxis([-maxneg,maxpos]);
%       caxis(climcomp);
    set(tightAxes_diff, 'xscale', 'log', 'yscale', 'log', 'colormap', CMAP, 'clim', climcomp);
    set(tightAxes_diff(2:end), 'YTickLabel', {});

    if(wavenumber0_wavelength1)
%         legend(h_kxkzth, 'location', 'northwest');
%         xlim([minx, max(xlim)]); ylim([minz, max(ylim)]);
      set(tightAxes_diff, 'xlim', [minx, max(xlim)], 'ylim', [minz, max(ylim)], 'xtick', logspace(log10(minx), ceil(log10(max(xlim))),ceil(log10(max(xlim)))-log10(minx)+1));
      xtl = split(sprintf('%.2g|',xticks),'|'); %xticklabels(xtl(1:end-1));
      ytl = split(sprintf('%.2g|',yticks),'|'); %yticklabels(ytl(1:end-1));
      set(tightAxes_diff(1), 'yticklabels', ytl(1:end-1));
      set(tightAxes_diff, 'xticklabels', xtl(1:end-1));
    else
%         legend(h_kxkzth, 'location', 'northeast');
      error('kek');
    end
%       xlabel(XLAB); ylabel(YLAB);
%       title(TIT);
    hcb.Position([2,4]) = tightAxes_diff(end).Position([2,4]);
    hcb.Position([1]) = sum(tightAxes_diff(end).Position([1,3]))+hshift_cb;

    ll = add_labels_subplots(fig_spectdiff,1);

    customSaveFig(fig_spectdiff, [figspectdiffpath], extToSave, 9999);
  else
    error('base file does not exist');
  end
%   end
end

% move produced figures to thesis
disp(['[] Starting to move Figures to thesis folder.']);
system(['cp ', figfieldpath, '.* /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/im_7_topo']);
system(['cp ', figspectpath, '.* /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/im_7_topo']);
system(['cp ', figspectdiffpath, '.* /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/im_7_topo']);

for i=1:numel(OFDs)
  system(['cp ', OFDs{i},filesep,'image0000005.jpg /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/im_7_topo']);
end
