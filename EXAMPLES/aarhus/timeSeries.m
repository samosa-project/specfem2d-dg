clear all;
close all;
clc;

addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools'));

cFldr = [regexprep(mfilename('fullpath'),mfilename,'')];
figfolder = [cFldr, filesep, 'FIGURES', filesep];
addpath(genpath(cFldr));

% plot parameters
distOverPTP = 0.75;
do_save_plot = 1;
normallise_ts = 1;
synth_on_nrg_plot = 0;

% data parameters
data_timestamp = '182608'; data_nrun = 397; data_pre_Pa = 20*100; data_freq = 2090;
data_filter_kind = 'bp';
% data_filter_fcut = [1150, 2923];
% data_filter_fcut = [1050, 3495];
data_filter_fcut = [958.7, 4085];
data_filter_order = 2;

% synthetics parameters
switch(data_freq)
  case 2090
    OFD = [cFldr,filesep,'OUTPUT_FILES_402141_2090_lns_allatt']; rsclFctr = 1.693495; tShft = (3.204-2504.925)*1e-3; TIT = 'LNS, with all attenuations'; addend='__lns_allatt';
%     OFD = [cFldr,filesep,'OUTPUT_FILES_402036_2090_fns_allatt']; rsclFctr = 1.693642; tShft = (3.204-2504.925)*1e-3; TIT = 'FNS, with all attenuations'; addend='__fns_allatt';
    
%     OFD = [cFldr,filesep,'OUTPUT_FILES_401255_2090_lns_acl']; rsclFctr = 36931.5321885834; tShft = (1.74760004-2504.2)*1e-3; TIT = 'LNS, w/ $\alpha_\mathrm{cl}\neq0$, $\alpha_\mathrm{rot}=\alpha_\mathrm{vib}=0$'; addend='__lns_acl';
%     OFD = [cFldr,filesep,'OUTPUT_FILES_401758_2090_lns_noatt']; rsclFctr = 1.309279e+00; tShft = (2.25799996-2504.2-5.068+5.309)*1e-3; TIT = 'LNS, w/ $\alpha_\mathrm{cl}=\alpha_\mathrm{rot}=\alpha_\mathrm{vib}=0$'; addend='__lns_noatt';
%     OFD = [cFldr,filesep,'OUTPUT_FILES_401736_2090_fns_acl_avib']; rsclFctr = 6.242541e+04; tShft = (2.25099991-2504.2)*1e-3; TIT = 'FNS, w/ $\alpha_\mathrm{cl}\neq0$, $\alpha_\mathrm{rot}=0$, $\alpha_\mathrm{vib}\neq0$'; addend='__fns_acl_avib';
%     OFD = [cFldr,filesep,'OUTPUT_FILES_401741_2090_fns_noatt']; rsclFctr = 5.191609e+04; tShft = (2.25099991-2504.2)*1e-3; TIT = 'FNS, w/ $\alpha_\mathrm{cl}=\alpha_\mathrm{rot}=\alpha_\mathrm{vib}=0$'; addend='__fns_noatt';
  otherwise
    error('kek');
end
istattab = 1:5; % should not change
typeDisplay = 2; % should not change
distChoice = 1; % should not change
doGeometricAttenuation = 0;
subsample = 1;
subsample_wantedDt = 1e-6;

% load synthetics
[Tsy, Psy, Dsy, Nsy, ~] = gji2020_loadSomeSynthetics(OFD, istattab, typeDisplay, 'BXZ', distChoice, doGeometricAttenuation, rsclFctr, subsample, subsample_wantedDt);
NSsy = size(Tsy,1);
Dsy = Dsy.vals';
Csy = repmat([1,0,0], NSsy, 1);

% load and stacc data
stacc_file = [figfolder,filesep,'data_stack.mat'];
if(not(exist(stacc_file, 'file')))
  warning off;
  [Tda, Pda, Dda] = Load_Data_ATN(data_timestamp, data_nrun, data_pre_Pa/100, data_freq);
  warning on;
  Pda = Pda';
  NSda = size(Pda, 1);
  Tda = repmat(Tda, NSda, 1);
  tstacctime = {};
  tstacctime{1} = [2.504205, 4.004375, 5.50361, 7.001905, 8.499265, 9.995695, 11.49119, 12.985755, 14.479395, 15.972095];
  tstacctime{2} = [2.50493, 4.0051, 5.504335, 7.00263, 8.499995, 9.996425, 11.49192, 12.986485, 14.48012, 15.972825];
  tstacctime{3} = [2.506775, 4.006945, 5.50618, 7.004475, 8.50184, 9.998265, 11.49376, 12.98833, 14.481965, 15.97467];
  tstacctime{4} = [2.51049, 4.01066, 5.509895, 7.008185, 8.505555, 10.00198, 11.49748, 12.992045, 14.485675, 15.97838];
  tstacctime{5} = [2.514205, 4.01437, 5.5136, 7.0119, 8.50926, 10.00569, 11.501185, 12.995745, 14.48939, 15.982085];
  [t_stack, p_stack_save, p_stack] = get_staccs(Tda, Pda, tstacctime, data_filter_kind, data_filter_fcut, data_filter_order);
  save(stacc_file, 'Dda', 't_stack', 'p_stack_save', 'p_stack');
else
  load(stacc_file);
  NSda = numel(t_stack);
end
ratio = range(p_stack{1})/range(Psy(1,:));
Cda = repmat([0,0,0], NSda, 1);
if(abs(1-ratio)>1e-6)
  disp(['[',mfilename,'] Ratio between data amplitude and synthetic amplitude (i.e. recommended scaling of synthetics): ',sprintf('%.6e', ratio),'.']);
  pause;
end

% Compute Acoustic NRG.
ANRGda = zeros(NSda, 1);
ANRGsy = zeros(NSsy, 1);
window = [0,16]*1e-3; % [ms]
for i = 1:NSda
  sel = (t_stack{i}+tShft>=min(window)) & (t_stack{i}+tShft<=max(window));
  ANRGda(i) = trapz(t_stack{i}(sel), p_stack{i}(sel).^2);
end
for i = 1:NSsy
  sel = (Tsy(i, :)>=min(window)) & (Tsy(i, :)<=max(window));
  ANRGsy(i) = trapz(Tsy(i, sel), Psy(i, sel).^2);
end
% compute amplitude of first peak
AFPda_t = zeros(NSda, 1);
AFPda = zeros(NSda, 1);
AFPsy_t = zeros(NSsy, 1);
AFPsy = zeros(NSsy, 1);
for i = 1:NSsy
  [AFPda_t(i), AFPda(i), ~] = findFirstPeak(t_stack{i}, max(p_stack{i},0), 0.35);
  [AFPsy_t(i), AFPsy(i), ~] = findFirstPeak(Tsy(i,:), max(Psy(i,:),0), 0.3);
end

% prepare a scale to print
idsy = 3;
tamp = 1e-3;
amp = range(Psy(idsy, :));
amp = round(amp,abs(ceil(log10(amp)))+3);

% plot
[pre_t, fac_t] = prefix_factor_values({Tsy(:)});
XLAB = ['time [',pre_t,'s]'];
XLIM = [0,16e-3]; XTIC = (0:1:16)*1e-3;
NAda = 'filtered data stack';
NAsy = 'synthetics';
if(normallise_ts)
  YLAB = ['distance to the source $r$ [m]'];
  distOverPTP = 0.5;
else
  YLAB = ['$r \mathrm{ [m]} + ',sprintf('%.2f',distOverPTP),' \mathrm{ [m/Pa]} \times p'' \mathrm{ [Pa]}$'];
end
fig_tvd = figure('units','normalized','outerposition',[0,0,1,1]);
tightAxes = tight_subplot(1, 1, [0, 0.05], [0.12,0.08], [0.07, 0.02]);
LS = '-';
handlesToLines = [];
for i = 1:NSda
  colour = Cda(i, :);
  if(normallise_ts)
    fac_y = 1/(range(p_stack{i}));
  end
  for ti=1:floor(size(p_stack_save{i}, 1)/2)
    plot(fac_t*(t_stack{i}+tShft), Dda(i) + distOverPTP * fac_y*p_stack_save{i}(ti, :), 'linewidth', 2, 'color', min((colour+[1,1,1])/2,1),'linestyle', LS); hold on;
  end
  h = plot(fac_t*(t_stack{i}+tShft), Dda(i) + distOverPTP * fac_y*p_stack{i}, 'color', colour, 'linestyle', LS, 'displayname', NAda);
%   plot(fac_t*(AFPda_t(i)+tShft)*[1,1], Dda(i)+[-1,1]*0.25, 'g');
  handlesToLines = [handlesToLines, h];
end
for i = 1:NSsy
  colour = Csy(i, :);
  if(normallise_ts)
    fac_y = 1/range(Psy(i, :));
  end
  h = plot(fac_t*Tsy(i, :), Dsy(i) + distOverPTP * fac_y*Psy(i, :), 'displayname', NAsy, 'color', colour,'linestyle',LS); hold on;
%   plot(fac_t*AFPsy_t(i)*[1,1], Dsy(i)+[-1,1]*0.25, 'g:');
  handlesToLines = [handlesToLines, h];
end
xlabel(XLAB); ylabel(YLAB);
xlim(fac_t*XLIM); xticks(fac_t*XTIC);
if(normallise_ts)
  legend([handlesToLines([1, NSda+1])], 'location', 'northwest');
else
  [pa, fa] = prefix_factor_values({amp});
  hscale = plot(fac_t*[1,1]*tamp, Dsy(idsy)+distOverPTP*[-1,1]*amp, 'linewidth', 6, 'color', [0,1,0]*0.5, 'displayname', [sprintf('%.2f',fa*amp),'~',pa,'Pa scale']);
  legend([handlesToLines([1, NSda+1]), hscale], 'location', 'northwest');
end
if(do_save_plot)
  name_fig_tvd = ['synth_v_data__f=',sprintf('%04d', data_freq), 'Hz_p=', sprintf('%04.0f', data_pre_Pa),'Pa'];
  customSaveFig(fig_tvd, [figfolder, name_fig_tvd, addend], {'fig', 'eps', 'png', 'tex'}, 9999);
end

% acoustic energy
alpha_th = 0.319137;
one_over_dist = 1./abs(1+Dda-min(Dda));
% do fits
fitt = 'c - 2*a*x - 2*b*log(x)'; fitfunc = fittype(fitt); op = fitoptions(fitfunc); op.Weights = one_over_dist.^2; fits={}; colfit={}; i = 1;
op.Lower = [0, 0, -Inf]; op.Upper = [Inf, 1, Inf]; fits{i} = fit(Dda', log(ANRGda), fitfunc, op); colfit{i}=[0,1,1]*0.5; i=i+1; % default
op.Lower = [0, 1, -Inf]; op.Upper = [Inf, 1, Inf]; fits{i} = fit(Dda', log(ANRGda), fitfunc, op); colfit{i}=[0,1,1]*0.1; i=i+1; % force spherical
op.Lower = [alpha_th, 1, -Inf]; op.Upper = [alpha_th, 1, Inf]; fits_nrg_th = fit(Dda', log(ANRGda), fitfunc, op); colfit_nrg_th=[1,0,0]; i=i+1; % theoretical
x = linspace(min(Dda), max(Dda), 40);
fitt_amp = 'c * exp(-a*x) / x.^b'; fitfunc_amp = fittype(fitt_amp); op_amp = fitoptions(fitfunc_amp); op_amp.Weights = one_over_dist; fits_amp={}; colfit_amp={}; i = 1;
op_amp.Lower = [0, 0, -Inf]; op_amp.Upper = [Inf, 1, Inf]; fits_amp{i} = fit(Dda', AFPda, fitfunc_amp, op_amp); colfit_amp{i}=[0,1,1]*0.5; i=i+1; % default
op_amp.Lower = [0, 1, -Inf]; op_amp.Upper = [Inf, 1, Inf]; fits_amp{i} = fit(Dda', AFPda, fitfunc_amp, op_amp); colfit_amp{i}=[0,1,1]*0.1; i=i+1; % force spherical
op.Lower = [alpha_th, 1, -Inf]; op.Upper = [alpha_th, 1, Inf]; fits_amp_th = fit(Dda', AFPda, fitfunc_amp, op); colfit_amp_th=[1,0,0]; i=i+1; % theoretical
x = linspace(min(Dda), max(Dda), 40);

% figure
fig_anrg = figure('units','normalized','outerposition',[0,0,1,1]);
tightAxes = tight_subplot(2, 1, [0.05, 0.05], [0.12,0.06], [0.09, 0.01]);
axes(tightAxes(1));
hfit=[];
for i=1:numel(fits)
  cfintcurfit = confint(fits{i});
  dnam = ['$\alpha=', sprintf('%.3f',fits{i}.a), '\pm',sprintf('%.3f',mean(abs(cfintcurfit(:,1)-fits{i}.a))),'$, $\beta=',sprintf('%.2f', fits{i}.b),'$'];
  semilogy(x, exp(fits{i}.c-2*max(cfintcurfit(1,1),0)*x-2*fits{i}.b*log(x)), 'displayname', dnam, 'linestyle', ':', 'color', colfit{i}); hold on;
  semilogy(x, exp(fits{i}.c-2*max(cfintcurfit(2,1),0)*x-2*fits{i}.b*log(x)), 'displayname', dnam, 'linestyle', ':', 'color', colfit{i}); hold on;
  hfit=[hfit,semilogy(x, exp(fits{i}(x)), 'displayname', dnam, 'color', colfit{i})]; hold on;
end
hfit=[hfit,semilogy(x, exp(fits_nrg_th(x)), 'displayname', ['$\alpha=\alpha_\mathrm{th}=', sprintf('%.6f',alpha_th), '$, $\beta=1.00$'], 'color', colfit_nrg_th)]; hold on;
XLIM = [0.25, 3.5]; ms = 250;
hpts = [];
hpts = [hpts, scatter(Dda, ANRGda, ms, 'o', 'markerfacecolor', Cda(1, :), 'markeredgecolor', 'none', 'displayname', NAda)]; hold on;
if(synth_on_nrg_plot)
  dnam_nrg_sy = [NAsy];
  hpts = [hpts, scatter(Dsy, ANRGsy, ms, 'o', 'markerfacecolor', Csy(1, :), 'markeredgecolor', 'none', 'displayname', dnam_nrg_sy)];
  hpts = [hpts, scatter(Dsy, ANRGsy.*(0.05./Dsy'.^2), ms*1.5, 'x', 'markeredgecolor', Csy(1, :), 'linewidth', 4, 'displayname', [NAsy, ', scaled $\times0.05/r^2$'])];
end
legend([hpts, hfit],'location', 'eastoutside');
ylabel('acoustic energy [Pa$^2\cdot$s]');
axes(tightAxes(2));
hpts = [];
hpts = [hpts, scatter(Dsy, AFPda, ms, '^', 'markerfacecolor', Cda(1, :), 'markeredgecolor', 'none', 'displayname', NAda)]; hold on;
hfit=[];
for i=1:numel(fits_amp)
  cfintcurfit = confint(fits_amp{i});
  dnam = ['$\alpha=', sprintf('%.3f',fits_amp{i}.a), '\pm',sprintf('%.3f',mean(abs(cfintcurfit(:,1)-fits_amp{i}.a))),'$, $\beta=',sprintf('%.2f', fits_amp{i}.b),'$'];
  semilogy(x, fits_amp{i}.c*exp(-max(cfintcurfit(1,1),0)*x)./x.^(fits_amp{i}.b), 'displayname', dnam, 'linestyle', ':', 'color', colfit{i}); hold on;
  semilogy(x, fits_amp{i}.c*exp(-max(cfintcurfit(2,1),0)*x)./x.^(fits_amp{i}.b), 'displayname', dnam, 'linestyle', ':', 'color', colfit{i}); hold on;
  hfit = [hfit, semilogy(x, fits_amp{i}(x), 'displayname', dnam, 'color', colfit{i})]; hold on;
end
hfit=[hfit,semilogy(x, fits_amp_th(x), 'displayname', ['$\alpha=\alpha_\mathrm{th}=', sprintf('%.6f',alpha_th), '$, $\beta=1.00$'], 'color', colfit_amp_th)]; hold on;
ylabel(['amplitude [Pa]']);
xlabel(['distance to the source $r$ [m]']);
legend([hpts, hfit],'location', 'eastoutside');
set(tightAxes, 'yscale', 'log');
set(tightAxes(1), 'xticklabels', {});
set(tightAxes(1), 'ylim', [1e-8, 1e-4]);
set(tightAxes(2), 'ylim', [5e-3, 2e-1]);
linkaxes(tightAxes, 'x');
xlim(XLIM);
linkprop(tightAxes,'xtick');
tightAxes(1).Position([3]) = tightAxes(2).Position([3]);

if(do_save_plot)
  name_fig_nrg = ['synth_v_data__f=',sprintf('%04d', data_freq), 'Hz_p=', sprintf('%04.0f', data_pre_Pa),'Pa__nrg'];
%   customSaveFig(fig_tvd, [figfolder, name_fig, addend], {'png'}, 9999);
  customSaveFig(fig_anrg, [figfolder, name_fig_nrg, addend], {'fig', 'eps', 'png', 'tex'}, 9999);
end
