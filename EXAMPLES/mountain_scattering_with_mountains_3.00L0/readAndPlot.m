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

do_pfield = 0;
do_fft = 1;
do_comparefft = 0;

pfield_zoombox_x = [-1,1]*14e3;
pfield_zoombox_z = [7,11.5]*1e3;

fft_wavenumber0_wavelength1 = 1; % plot in terms of wavenumber (0) or wavelength (1)
fft_pcolor1_contour0 = 0; % contour prettier, but slower

pfield_clim = [-1, 1]*175; thrsh_pfield = 0.15; blk_pfield=0.95;
fft_clim = [-1.5, 0];
% fftcomp_clim = [-1, 1]*5.5; maxneg_fftcomp = abs(min(fftcomp_clim)); maxpos_fftcomp = abs(max(fftcomp_clim)); thresh_fftcomp = max([maxneg_fftcomp, maxpos_fftcomp])-0.5;
fftcomp_clim = [-1, 1]*1; %maxneg_fftcomp = abs(min(fftcomp_clim)); maxpos_fftcomp = abs(max(fftcomp_clim)); thresh_fftcomp = max([maxneg_fftcomp, maxpos_fftcomp])-0.5;
fft_xtick = [0.1,1,10]; fft_ytick = fft_xtick;
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

if(max(abs(diff(dt.*IT))) > 1e-9)
  error('dt*IT varies from one simulation to the other');
end

% Figure.
if(do_pfield)
  CBYLAB_PRESPEC = ['field of pressure perturbations $p''$ at $t=',num2str(IT(1)*dt(1)),'$~s [Pa]'];
  fig_field = figure('units','normalized','outerposition',[0,0,1,1]);
  tightAxes = tight_subplot(4, 1, [0.03,0.], [0.11,0.05], [0.06, 0.06]);
  for i=IDs_to_process
    axes(tightAxes(i));
    selx = (Xi{i}(1,:)>=min(pfield_zoombox_x)*1.1) & (Xi{i}(1,:)<=max(pfield_zoombox_x)*1.1);
    selz = (Yi{i}(:,1)>=min(pfield_zoombox_z)*0.9) & (Yi{i}(:,1)<=max(pfield_zoombox_z)*1.1);
    pcolor(Xi{i}(selz,selx)/1e3, Yi{i}(selz,selx)/1e3, Vi{i}.pre(selz,selx)); hold on;
    ylabel(['$z$ [km]']);
  end
  xlabel(['$x$ [km]']);
  hcb = colorbar(); ylabel(hcb, CBYLAB_PRESPEC, 'interpreter', 'latex');
  XTCK_pfield = (min(pfield_zoombox_x):2e3:max(pfield_zoombox_x))/1e3;
  ZTCK_pfield = (min(pfield_zoombox_z):2e3:max(pfield_zoombox_z))/1e3;
  CMAP = colormaps_custom([-1,-blk_pfield, -(1+thrsh_pfield)/2, -thrsh_pfield, 0, thrsh_pfield, (1+thrsh_pfield)/2, blk_pfield,1], [[0,0,1].*[0.25,0.75,1]';[0.9,0.9,1];[1,1,1];[1,0.9,0.9];[1,0,0].*[1,0.75,0.25]'], 0);
  set(tightAxes, 'CLim', pfield_clim, 'Colormap', CMAP, 'XLim', pfield_zoombox_x/1e3, 'YLim', pfield_zoombox_z/1e3, 'XTick', XTCK_pfield, 'YTick', ZTCK_pfield);
  for i=1:numel(tightAxes)
    tightAxes(i).Position = [tightAxes(end).Position(1), tightAxes(i).Position(2), tightAxes(4).Position(3:4)];
  end
  set(tightAxes(1:end-1), 'XTickLabel', {});
  height_spanning_all_axes = sum(tightAxes(1).Position([2,4]))-tightAxes(end).Position(2);
  set(hcb, 'Position', [hcb.Position(1) + 0.005, tightAxes(end).Position(2), hcb.Position(3), height_spanning_all_axes],'fontsize', 26);
  ll = add_labels_subplots_local(fig_field,tightAxes,0.9);
  customSaveFig(fig_field, [figfieldpath], extToSave, 9999);
end

% Compute FFT.
PRE_fft1s = {};
if(any([do_fft, do_comparefft]))
  for i=IDs_to_process
    x = unique(Xi{i}); z = unique(Yi{i});
    [kx, kz, PRE_fft1s{i}, ~, ~, ~] = fft2_wrap(x, z, Vi{i}.pre);
    if(i==1)
      baseline_fft = PRE_fft1s{i};
      save([baseline_fft_path], 'baseline_fft');
      disp(['saving baseline']);
    end
  end
  
  % Plot parameters.
  if(fft_wavenumber0_wavelength1)
    fac = 1e-3; XLAB = 'horizontal wavelength $(1/k_x)$ [km]'; YLAB = 'vertical wavelength $(1/k_z)$ [km]';
    pow = -1;
    minx = 50*fac; minz = minx;
    maxx = round(max(kx(kx>0).^pow)*fac, -1);
  else
    XLAB = '$k_x$ [1/m]';
    YLAB = '$k_z$ [1/m]'; pow=1;
  end
  margz = [0.16, 0.04]; margh = [0.06, 0.085]; gap=[0, 0.015];
  hshift_cb = 0.04;
  absfftname = ['\left|\widehat{P}\left(k_x,k_z\right)\right|'];
  CBYLAB_PRESPEC = ['$\log_{10}\left(',absfftname,'\right)$'];
%   CBYLAB_diff = ['signed difference of ',CBYLAB_PRESPEC];
  CBYLAB_diff = ['$',absfftname,'$ difference'];
end

% Build theoretical curves.
LW_thkxkz = 3;
c = 342;
vp = 6078.06;
vs = 3537.23;
rho = 2730.14;
f0 = 2;
[E_nu] = conversion('r', rho, 'p', vp, 's', vs, 'E', 'n');
nu = E_nu(2);
vrayleigh = vs / ((1+nu)/(0.862+1.14*nu));
thcurvs = {}; i=1;
thcurvs{i}.Lx = [0.1,10]; thcurvs{i}.Lz = fac*[1,1]*c/f0; thcurvs{i}.col = [1,0,0]*0.66; thcurvs{i}.ls = '-'; i=i+1;
thcurvs{i}.Lx = [0.1,10]; thcurvs{i}.Lz = fac*[1,1]*c*vp/(f0*vrayleigh); thcurvs{i}.col = [1,0,0]*0.66; thcurvs{i}.ls = '--'; i=i+1;

% Plot FFT.
if(do_fft)
  fig_spect = figure('units', 'normalized', 'outerposition', [0,0,1,0.7]);
  tightAxes_spec = tight_subplot(1, 4, gap, margz, margh);
  for i=IDs_to_process
    axes(tightAxes_spec(i));
    toPlot = log10(abs(PRE_fft1s{i}));
    toPlot = toPlot(kz>0, kx>0);
    if(fft_pcolor1_contour0)
      pcolor(fac*kx(kx>0).^pow, fac*kz(kz>0).^pow, toPlot); hold on;
    else
      contourf(fac*kx(kx>0).^pow, fac*kz(kz>0).^pow, toPlot, [min(fft_clim):0.25:max(fft_clim)], 'edgecolor', 'none'); hold on;
    end
    for j=1:numel(thcurvs)
%       h_kxkzth = plot(fac*[min(kx(kx>0)), [1,1]*th_kxkz{i}].^pow, fac*[[1,1]*th_kxkz{i}, min(kz(kz>0))].^pow, 'color', COL_thkxkz, 'linewidth', LW_thkxkz, 'displayname', dnam_thkxkz);
      plot(thcurvs{j}.Lx, thcurvs{j}.Lz, 'color', thcurvs{j}.col, 'linewidth', LW_thkxkz, 'linestyle', thcurvs{j}.ls);
    end
    if(i==2)
      xlabel(XLAB);
    end
    if(i==1)
      ylabel(YLAB);
    end
  end
  hcb1 = colorbar(); ylabel(hcb1, CBYLAB_PRESPEC, 'interpreter', 'latex');
  CMAP = colormaps_fromPython('bone', 0);
  CMAP = flipud(CMAP);
  set(tightAxes_spec, 'xscale', 'log', 'yscale', 'log', 'xlim', [minx, maxx], 'ylim', [minz, max(ylim)], 'xtick', fft_xtick, 'ytick', fft_xtick, 'colormap', CMAP, 'clim', fft_clim);
  set(tightAxes_spec(2:end), 'YTickLabel', {});
  if(fft_wavenumber0_wavelength1)
    xtl = split(sprintf('%.3g|',xticks),'|');
    ytl = split(sprintf('%.3g|',yticks),'|');
    set(tightAxes_spec(1), 'yticklabels', ytl(1:end-1));
    set(tightAxes_spec, 'xticklabels', xtl(1:end-1));
  else
    legend(h_kxkzth, 'location', 'northeast');
    error('kek');
  end
  hcb1.Position([2,4]) = tightAxes_spec(end).Position([2,4]);
  hcb1.Position([1]) = sum(tightAxes_spec(end).Position([1,3]))+hshift_cb;
  ll = add_labels_subplots_local(fig_spect, tightAxes_spec, 1);
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
%       sgn = sign(toPlot);
%       toPlot = log10(toPlot./sgn).*sgn;
      toPlot = toPlot(kz>0, kx>0);
      if(fft_pcolor1_contour0)
        pcolor(fac*kx(kx>0).^pow, fac*kz(kz>0).^pow, toPlot); hold on;
      else
        step = 0.25;
        contourf(fac*kx(kx>0).^pow, fac*kz(kz>0).^pow, toPlot, [min(fftcomp_clim-step):step:max(fftcomp_clim+step)], 'edgecolor', 'none');
      end
      if(i==3)
        xlabel(XLAB);
      end
      if(i==2)
        ylabel(YLAB);
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

% move produced figures to thesis
disp(['[] Starting to move Figures to thesis folder.']);
system(['cp ', figfieldpath, '.* /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/im_7_topo']);
system(['cp ', figspectpath, '.* /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/im_7_topo']);
system(['cp ', figspectdiffpath, '.* /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/im_7_topo']);

for i=1:numel(OFDs)
  system(['cp ', OFDs{i},filesep,'image0000005.jpg /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/im_7_topo/snap_',TAGS{i},'.jpg']);
end





function [labels] = add_labels_subplots_local(figure_handle, list_of_axes, scale_factor, offset, xzshift)
  if(not(exist('scale_factor', 'var')))
    scale_factor = 1;
  end
  if(not(exist('offset', 'var')))
    offset = 0;
  end
  if(not(exist('xzshift', 'var')))
    xzshift = [1,1]*0;
  end
  
%   labtype = 'roman'; % a), b), ...
  labtype = 'numeric'; % 1), 2), ...
  
%   children = figure_handle.Children;
%   
%   list_of_axes = [];
  AX_leftCoordinate = [];
  AX_topCoordinate = [];
%   c = 1;
  for c=1:numel(list_of_axes)
%     if(strcmp(class(children(i)),'matlab.graphics.axis.Axes'))
%       list_of_axes = [list_of_axes, children(i)];
      AX_leftCoordinate(c) = list_of_axes(c).Position(1); % left
  %     list_of_axes_y(c) = list_of_axes(c).Position(2); % bottom
      AX_topCoordinate(c) = sum(list_of_axes(c).Position([2,4])); % top
%       c = c+1;
%     end
  end
% %   list_of_axes'
% 
%   % Sort from top to bottom first.
%   [~, isort_y] = sort(AX_topCoordinate, 'descend');
%   list_of_axes = list_of_axes(isort_y);
%   AX_leftCoordinate = AX_leftCoordinate(isort_y);
%   AX_topCoordinate = AX_topCoordinate(isort_y);
% %   list_of_axes'
% 
%   % Find axes which have same y but different x, and sort them left to right.
%   isort_x = 1:numel(list_of_axes);
%   unique_top_coords = unique(AX_topCoordinate);
%   for tt=1:numel(unique_top_coords)
%     original_ids = find(AX_topCoordinate==unique_top_coords(tt));
%     x_to_sort = AX_leftCoordinate(AX_topCoordinate==unique_top_coords(tt));
%     [~, x_sorted_ids] = sort(x_to_sort);
%     isort_x(original_ids) = original_ids(x_sorted_ids);
%   end
%   list_of_axes = list_of_axes(isort_x);
%   AX_leftCoordinate = AX_leftCoordinate(isort_x);
%   AX_topCoordinate = AX_topCoordinate(isort_x);
% %   list_of_axes'
%   
  % Prepare size of annotations.
  set(figure_handle,'units','pixels');
  figsiz = get(figure_handle, 'Position');
  set(figure_handle,'units','normalized');
%   hsiz = 0.025; % proportion
  hsiz = 0.03*scale_factor; % proportion
%   hsiz = 0.035; % proportion
  vsiz = hsiz*figsiz(3)/figsiz(4);

  % Create annotations.
  labels = [];
  for i=1:numel(list_of_axes)
    switch(labtype)
      case 'roman'
        lab = [char(96+offset+i),')'];
      case 'numeric'
        lab = ['S',num2str(offset+i)];
      otherwise
        error([' label type unknown']);
    end

    labels(i) = annotation('textbox',[min(AX_leftCoordinate(i)+xzshift(1),1), min(AX_topCoordinate(i)-vsiz+xzshift(2),1), hsiz, vsiz], ...
                          'String',lab, ...
                          'FitBoxToText', 'off', ...
                          'color', 'k', ...
                          'backgroundcolor', 'w', ...
                          'interpreter', 'latex', ...
                          'horizontalalignment', 'center', ...
                          'linewidth', get(list_of_axes(i),'linewidth')*0.8);
  end
end
