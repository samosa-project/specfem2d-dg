do_load = 0;
if(do_load)
  clear all;
  do_load = 1;
end
close all;
clc;

addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils'));

rtdwprfx = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/mountain_scattering_';
thisFolder = [regexprep(mfilename('fullpath'),mfilename,'')];
output_figures_folder = [thisFolder, filesep, 'FIGURES'];
baseline_fft_path = [output_figures_folder,filesep,'BASELINE_PRESSURE_SPECTRUM'];
figfieldpath = [output_figures_folder, filesep, 'all_fields'];
figspectpath = [output_figures_folder, filesep, 'all_spects'];
figspectdiffpath = [output_figures_folder, filesep,'all_spectdiffs'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tg_bsln = '4214695_lns'; tg_03hig = '4214696_lns'; tg_03low = '4219868_lns'; tg_3 = '4214697_lns'; IT = 100000;
% tg_bsln = '398709_fns'; tg_03hig = '4214696_lns'; tg_03low = '4219868_lns'; tg_3 = '398703_fns'; IT = 140000;
tg_bsln = '408827_lns'; tg_03hig = '408932_lns'; tg_03low = '409807_lns'; tg_3 = '408826_lns'; IT = [100, 240, 100, 100]*1e3;

sel_boxabsx = 28e3;
sel_boxy = [1501, 15e3];

interp_forceDGMesh = 0;
interp_dx = 20;
interp_dz = interp_dx; % dx=10 ok, dx<10 chugs hard

do_pfield = 1;
do_fft = 0;
do_comparefft = 0;

pfield_zoombox_x = [-1,1]*14e3;
pfield_zoombox_z = [7,11.5]*1e3;

fft_wavenumber0_wavelength1 = 1; % plot in terms of wavenumber (0) or wavelength (1)
fft_pcolor1_contour0 = 0; % contour prettier, but slower

pfield_clim = [-1, 1]*175; thrsh_pfield = 0.15; blk_pfield=0.95;
fft_clim = [-2, 0];
fft_colourISBand=[0,1,0]*0.85;
fftcomp_clim = [-1, 1]*1;
fft_xlim = [0.1, 10]; fft_xtick = [0.1,1,10]; fft_xticklabel = {'0.1', '1', '10'};
fft_ylim = fft_xlim; fft_ytick = fft_xtick;
fft_xlim_ang = [0, 90]; fft_xtick_ang = 0:15:90;
fft_ylim_fre = [1e-1, 3]; fft_ytick_fre = logspace(-2,0,3);
fft_selFreBand = [0.5, 2]; fft_selFreBandN = 10;
extToSave = {'eps'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OFDs = {[rtdwprfx, 'without_mountains', filesep, 'OUTPUT_FILES_', tg_bsln], ...
        [rtdwprfx, 'with_mountains_0.33L0', filesep, 'OUTPUT_FILES_', tg_03hig], ...
        [rtdwprfx, 'with_mountains_0.33L0_lower', filesep, 'OUTPUT_FILES_', tg_03low], ...
        [rtdwprfx, 'with_mountains_3.00L0', filesep, 'OUTPUT_FILES_', tg_3]};
TAGS = {'0p000', '0p33h', '0p33l', '3p00h'};
IDs_to_process = 1:4;

[S] = zint(max(abs(pfield_zoombox_x)));
fft_radPatAngles = [0, atan(pi*[S{2}.H, S{3}.H, S{4}.H]./[S{2}.L, S{3}.L, S{4}.L])]*180/pi;

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
    
    wavefield01 = [OFD,filesep,'wavefield',sprintf('%07d',IT(i)),'_01_000.txt'];
    matfile = [OFD,filesep,'wavefield',sprintf('%09d',IT(i)),'.mat'];
    chkDateMat = dir(matfile);
    chkDateWVF = dir(wavefield01);
    if(exist(matfile,'file') & datenum(chkDateMat.date)>datenum(chkDateWVF.date))
      disp(['[',mfilename,'] Matfile exists and is more recent than .txt file, loading it instead']);
      load(matfile, 'X', 'Y', 'V');
    else
      disp(['[',mfilename,'] No matfile found, or is older than .txt file, loading .txt and saving to mat']);
      % Read dumps.
      disp(['[',mfilename,'] Reading dumps (may be long).']);
      [X, Y, V] = readDumpsUnique(OFD, IT(i), 0);
      save(matfile, 'X', 'Y', 'V'); % save to mat for later
    end

    % Select box.
    disp(['[',mfilename,'] Selecting values within box.']);
    select = ((Y>=min(sel_boxy)) & (Y<=max(sel_boxy)) & (abs(X)<=sel_boxabsx));
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
%     % Get theoretical wavelengths.
%     sourcefile = [OFD, filesep, 'input_source'];
%     f0 = readExampleFiles_extractParam(sourcefile, 'f0', 'float');
%     th_kxkz{i} = f0/soundspeed{i};
      
    USE_LNS = readExampleFiles_extractParam(parfile, 'USE_LNS', 'bool');
    if(not(USE_LNS))
      disp(['[',mfilename,'] Removing background atmophere.']);
      V.pre = V.pre - bg_pre; % Note: this probably fuccs up below ground, but we do not care.
    end

    % Interpolate for plotting.
    disp(['[',mfilename,'] Interpolating.']);
    [Xi{i}, Yi{i}, Vi{i}] = interpDumps(X, Y, V, range(X)/interp_dx, range(Y)/interp_dz, interp_forceDGMesh);
  end
end

% Figure.
if(do_pfield)
  if(max(abs(diff(dt.*IT))) > 1e-9)
    error('dt*IT varies from one simulation to the other');
  end
  CBYLAB_PRESPEC = ['field of pressure perturbations $p''$ [Pa], at $t=',num2str(IT(1)*dt(1)),'$~s'];
  fig_field = figure('units','normalized','outerposition',[0,0,0.8,1]);
  tightAxes = tight_subplot(8, 1, [0.03,0.], [0.09,0.01], [0.06, 0.12]);
  for i=IDs_to_process
    axes(tightAxes(i*2-1));
    selx = (Xi{i}(1,:)>=min(pfield_zoombox_x)*1.1) & (Xi{i}(1,:)<=max(pfield_zoombox_x)*1.1);
    selz = (Yi{i}(:,1)>=min(pfield_zoombox_z)*0.9) & (Yi{i}(:,1)<=max(pfield_zoombox_z)*1.1);
    pcolor(Xi{i}(selz,selx)/1e3, Yi{i}(selz,selx)/1e3, Vi{i}.pre(selz,selx)); hold on;
    if(i==3); ylabel(['altitude $z$ [km]']); end
    if(i==4); hcb = colorbar(); ylabel(hcb, CBYLAB_PRESPEC, 'interpreter', 'latex'); end
    
    axes(tightAxes(i*2));
    plot(S{i}.x/1e3, S{i}.z/1e3);
%     ylabel(['$z$ [km]']);
  end
  xlabel(['$x$ [km]']);
  XTCK_pfield = (min(pfield_zoombox_x):2e3:max(pfield_zoombox_x))/1e3;
  ZTCK_pfield = (min(pfield_zoombox_z):2e3:max(pfield_zoombox_z))/1e3;
  linkaxes(tightAxes, 'x');
  CMAP = colormaps_custom([-1,-blk_pfield, -(1+thrsh_pfield)/2, -thrsh_pfield, 0, thrsh_pfield, (1+thrsh_pfield)/2, blk_pfield,1], [[0,0,1].*[0.25,0.75,1]';[0.9,0.9,1];[1,1,1];[1,0.9,0.9];[1,0,0].*[1,0.75,0.25]'], 0);
  set(tightAxes, 'tickdir', 'out');
  set(tightAxes((1:4)*2-1), 'CLim', pfield_clim, 'Colormap', CMAP, 'XLim', pfield_zoombox_x/1e3, 'YLim', pfield_zoombox_z/1e3, 'XTick', XTCK_pfield, 'YTick', ZTCK_pfield);
  set(tightAxes((1:4)*2), 'xlim', pfield_zoombox_x/1e3, 'ylim', [-100, 1700]/1e3, 'XTick', XTCK_pfield, 'YTick', [0,1500]/1e3);
  for i=1:numel(tightAxes); tightAxes(i).Position = [tightAxes(end).Position(1), tightAxes(i).Position(2), tightAxes(4).Position(3:4)]; end
  vreduce = 0.05;
  for i=((1:4)*2); tightAxes(i).Position(4) = tightAxes(i).Position(4)-vreduce; end
  for i=((1:4)*2-1); tightAxes(i).Position([2,4]) = tightAxes(i).Position([2,4]) + [-1,1]*(vreduce+0.01); end
  set(tightAxes(1:end-1), 'XTickLabel', {});
  for i=1:numel(tightAxes); set(tightAxes(i).YAxis, 'TickLength', [0.001, 0.002]*4); set(tightAxes(i).XAxis, 'TickLength', [0.001, 0.002]*4); end
  height_spanning_all_axes = sum(tightAxes(1).Position([2,4]))-tightAxes(end).Position(2);
  set(hcb, 'Position', [hcb.Position(1) + 0.005, tightAxes(end).Position(2), hcb.Position(3), height_spanning_all_axes],'fontsize', 26);
  ll = add_labels_subplots_local(fig_field, tightAxes((1:4)*2-1), 0.9);
  customSaveFig(fig_field, [figfieldpath], extToSave, 9999);
end

% Compute FFT.
PRE_fft1s = {};
LaunchAngles = {}; Frequencies = {};
V_probed_ang = {}; V_probed_fre = {}; V_probed_band = {};
if(any([do_fft, do_comparefft]))
  for i=IDs_to_process
    x = unique(Xi{i}); z = unique(Yi{i});
    [kx, kz, PRE_fft1s{i}, ~, ~, ~] = fft2_wrap(x, z, Vi{i}.pre);
    [KX, KZ] = meshgrid(kx(kx>0), kz(kz>0)); 
    [LaunchAngles{i}, Frequencies{i}] = KxKz2LaFr(1./KX, 1./KZ, soundspeed{i});
    [V_probed_ang{i}, V_probed_fre{i}, ~, V_probed_band{i}] = LaFr2RadPattern(LaunchAngles{i}, Frequencies{i}, abs(PRE_fft1s{i}(kz>0, kx>0)), min(fft_selFreBand), max(fft_selFreBand), fft_selFreBandN);
    if(i==1)
      baseline_fft = PRE_fft1s{i};
      save([baseline_fft_path], 'baseline_fft');
      disp(['saving baseline']);
    end
  end
  
  % Plot parameters.
  if(fft_wavenumber0_wavelength1)
    pow = -1;
    fac = 1e-3;
    XLAB1 = 'horizontal wavelength $(2\pi/k_x)$ [km]';
    YLAB1 = {'vertical wavelength','$(2\pi/k_z)$ [km]'};
    XLAB2 = 'launch angle [deg]';
    YLAB2 = 'frequency [Hz]';
    YLAB3 = {'radiated','infrasound'};
  end
  margz = [0.09, 0.025]; margh = [0.075, 0.09]; gap = [0.13, 0.016];
  hshift_cb = 0.01;
  absfftname = ['\left|\widehat{P}\left(k_x,k_z\right)\right|'];
  CBYLAB_PRESPEC = ['radiated infrasound, $\log_{10}\left(',absfftname,'\right)$'];
  CBYLAB_diff = ['$',absfftname,'$ difference'];
  
  % Build theoretical curves.
  LW_thkxkz = 3;
  c = soundspeed{1}; vp = 6078.06; vs = 3537.23; rho = 2730.14; f0 = 2;
  [E_nu] = conversion('r', rho, 'p', vp, 's', vs, 'E', 'n'); nu = E_nu(2);
  vrayleigh = vs / ((1+nu)/(0.862+1.14*nu));
  thcurvs = {}; i=1;
  thcurvs{i}.Lx = [0.1,10]; thcurvs{i}.Lz = fac*[1,1]*c/f0; thcurvs{i}.col = [1,0,0]*0.66; thcurvs{i}.ls = '-'; i=i+1;
  thcurvs{i}.Lx = [0.1,10]; thcurvs{i}.Lz = fac*[1,1]*c*vp/(f0*vrayleigh); thcurvs{i}.col = [1,0,0]*0.66; thcurvs{i}.ls = '--'; i=i+1;
end

% Plot FFT.
if(do_fft)
  fig_spect = figure('units', 'normalized', 'outerposition', [0,0,1,1]);
  tightAxes_spec = tight_subplot(2, 4, gap, margz, margh);
  tightAxes_radPat = tight_subplot(1, 4, gap, [margz(1), 0.78], margh);
  for i=IDs_to_process
    % start with wavenumbers
    axes(tightAxes_spec(i));
    toPlot = log10(abs(PRE_fft1s{i}));
    if(fft_pcolor1_contour0)
      pcolor(fac*kx(kx>0).^pow, fac*kz(kz>0).^pow, toPlot(kz>0, kx>0)); hold on;
    else
      contourf(fac*kx(kx>0).^pow, fac*kz(kz>0).^pow, toPlot(kz>0, kx>0), [min(fft_clim):0.25:max(fft_clim)], 'edgecolor', 'none'); hold on;
    end
    for j=1:numel(thcurvs)
%       h_kxkzth = plot(fac*[min(kx(kx>0)), [1,1]*th_kxkz{i}].^pow, fac*[[1,1]*th_kxkz{i}, min(kz(kz>0))].^pow, 'color', COL_thkxkz, 'linewidth', LW_thkxkz, 'displayname', dnam_thkxkz);
      plot(thcurvs{j}.Lx, thcurvs{j}.Lz, 'color', thcurvs{j}.col, 'linewidth', LW_thkxkz, 'linestyle', thcurvs{j}.ls);
    end
    if(i==2); xlabel(XLAB1); end
    if(i==1); ylabel(YLAB1); end
    
    % now plot launch angles and frequencies
    axes(tightAxes_spec(i + 4));
    if(fft_pcolor1_contour0)
      pcolor(LaunchAngles{i}*180/pi, Frequencies{i}, toPlot(kz>0, kx>0)); hold on;
    else
      contourf(LaunchAngles{i}*180/pi, Frequencies{i}, toPlot(kz>0, kx>0), [min(fft_clim):0.25:max(fft_clim)], 'edgecolor', 'none'); hold on;
    end
    for j=1:numel(fft_selFreBand)
      plot(fft_xlim_ang, [1,1]*fft_selFreBand(j), 'color', fft_colourISBand);
    end
    if(i==1); ylabel(YLAB2); end
    if(i==4); pos=get(gca,'Position'); hcb1 = colorbar(); ylabel(hcb1, CBYLAB_PRESPEC, 'interpreter', 'latex'); set(gca,'Position',pos); end
    
    % now plot radiation pattern in infrasound band
    axes(tightAxes_radPat(i))
    plot(fft_radPatAngles(i)*[1,1], 10.^fft_clim, ':', 'color', 'k'); hold on;
    semilogy(V_probed_ang{i}*180/pi, V_probed_band{i}, 'color', fft_colourISBand); hold on;
    if(i==1); ylabel(YLAB3); end
    if(i==2); xlabel(XLAB2); end
  end
  
  CMAP = colormaps_fromPython('bone', 0);
  CMAP = flipud(CMAP);
  set(tightAxes_spec, 'colormap', CMAP, 'clim', fft_clim);
  set(tightAxes_spec(1:4), 'xscale', 'log', 'yscale', 'log', 'xlim', fft_xlim, 'ylim', fft_ylim, 'ytick', fft_xtick, 'xtick', fft_xtick, 'yticklabel', fft_xticklabel);
  set(tightAxes_spec(1:3), 'xticklabel', [fft_xticklabel(1:end-1), ' ']); set(tightAxes_spec(4), 'xticklabel', fft_xticklabel);
  set(tightAxes_spec(5:8), 'xscale', 'lin', 'yscale', 'log', 'xlim', fft_xlim_ang, 'ylim', fft_ylim_fre);
  set(tightAxes_spec(5:8), 'xtick', fft_xtick_ang, 'ytick', fft_ytick_fre);
  set(tightAxes_radPat, 'xscale', 'lin', 'yscale', 'log', 'xtick', fft_xtick_ang, 'xlim', fft_xlim_ang, 'ylim', 10.^fft_clim);
  set(tightAxes_radPat([2:4]), 'YTickLabel', {});
  set(tightAxes_spec([2:4, 6:8]), 'YTickLabel', {});
  set(tightAxes_spec([5:8]), 'XTickLabel', {});
  hcb1.Position(1) = sum(tightAxes_spec(8).Position([1,3]))+hshift_cb;
  hcb1.Position(2) = tightAxes_spec(8).Position(2);
  hcb1.Position(4) = sum(tightAxes_spec(4).Position([2,4])) - tightAxes_spec(8).Position(2);
  for i=5:8; tightAxes_spec(i).Position([2,4]) = tightAxes_spec(i).Position([2,4]) + [1,-1]*(0.03 + sum(tightAxes_radPat(1).Position([2,4]))-tightAxes_spec(i).Position(2)); end
  for i=1:4
    vshift = 0.075;
    tightAxes_spec(i).Position([2,4]) = tightAxes_spec(i).Position([2,4]) + [1,-1]*vshift;
    tightAxes_spec(i+4).Position([4]) = tightAxes_spec(i+4).Position([4]) + vshift;
  end
  ll1 = add_labels_subplots_local(fig_spect, tightAxes_spec(1:4), 1);
%   ll2 = add_labels_subplots_local(fig_spect, tightAxes_spec(5:8), 1);
  customSaveFig(fig_spect, [figspectpath], extToSave, 9999);
end

% Compare.
if(do_comparefft)
  basefile = [baseline_fft_path,'.mat'];
  if(exist(basefile,'file'))
    PRE_fft1s_BASE = load(basefile);
    PRE_fft1s_BASE = PRE_fft1s_BASE.baseline_fft;
    fig_spectdiff = figure('units','normalized','outerposition',[0,0,1,0.7]);
    margh = margh + [0, +0.005];
    tightAxes_diff = tight_subplot(1, 3, gap, margz, margh);
    for i=IDs_to_process(2:end)
      axes(tightAxes_diff(i-1));
      toPlot = abs(PRE_fft1s{i})-abs(PRE_fft1s_BASE); % signed difference
      toPlot = toPlot(kz>0, kx>0);
      if(fft_pcolor1_contour0)
        pcolor(fac*kx(kx>0).^pow, fac*kz(kz>0).^pow, toPlot); hold on;
      else
        step = 0.25;
        contourf(fac*kx(kx>0).^pow, fac*kz(kz>0).^pow, toPlot, [min(fftcomp_clim-step):step:max(fftcomp_clim+step)], 'edgecolor', 'none');
      end
      if(i==3)
        xlabel(XLAB1);
      end
      if(i==2)
        ylabel(YLAB1);
      end
    end
    hcb = colorbar(); ylabel(hcb, CBYLAB_diff, 'interpreter', 'latex');
%     tickcmap = [-maxneg_fftcomp, [-1,-0.85,0,0.85,1]*thresh_fftcomp, maxpos_fftcomp];
%     CMAP = colormaps_custom(tickcmap, [[0,0,1].*[0.25,1]';[0.9,0.9,1];[1,1,1];[1,0.9,0.9];[1,0,0].*[1,0.25]'], 0);
    CMAP = colormaps_fromPython('seismic', 0);
    set(tightAxes_diff, 'xscale', 'log', 'yscale', 'log', 'xlim', [minx, maxx], 'ylim', [minz, max(ylim)], 'xtick', fft_xtick, 'ytick', fft_ytick, 'colormap', CMAP, 'clim', fftcomp_clim);
    set(tightAxes_diff(2:end), 'YTickLabel', {});
    if(fft_wavenumber0_wavelength1)
      xtl = split(sprintf('%.3g|',xticks),'|');
      ytl = split(sprintf('%.3g|',yticks),'|');
      set(tightAxes_diff(1), 'yticklabels', ytl(1:end-1));
      set(tightAxes_diff, 'xticklabels', xtl(1:end-1));
    else
      error('kek');
    end
    hcb.Position([2,4]) = tightAxes_diff(end).Position([2,4]);
    hcb.Position([1]) = sum(tightAxes_diff(end).Position([1,3]))+hshift_cb;
%     set(hcb, 'ticks', round(tickcmap,1));
    ll = add_labels_subplots_local(fig_spectdiff,tightAxes_diff,1,1);
    customSaveFig(fig_spectdiff, [figspectdiffpath], extToSave, 9999);
  else
    error('base file does not exist');
  end
end

% % move produced figures to thesis
% disp(['[] Starting to move Figures to thesis folder.']);
% system(['cp ', figfieldpath, '.* /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/im_7_topo']);
% system(['cp ', figspectpath, '.* /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/im_7_topo']);
% system(['cp ', figspectdiffpath, '.* /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/im_7_topo']);
% if(do_mosaic)
%   system(['cp ', figpath_mosaic, '.* /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/im_7_topo/']);
% else
%   for i=1:numel(OFDs)
%     system(['cp ', OFDs{i},filesep,'image0000005.jpg /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/im_7_topo/snap_',TAGS{i},'.jpg']);
%   end
% end
% 



