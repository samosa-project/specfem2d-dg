clear all;
close all;
clc;

addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools'));

cFldr = [regexprep(mfilename('fullpath'),mfilename,'')];
figfolder = [cFldr, filesep, 'FIGURES', filesep];
addpath(genpath(cFldr));

% plot parameters
distOverPTP = 2;
do_save_plot = 1;

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
%     OFD = [cFldr,filesep,'OUTPUT_FILES_401255_2090_lns_acl']; rsclFctr = 36931.5321885834; tShft = (1.74760004-2504.2)*1e-3; TIT = 'LNS, w/ $\alpha_\mathrm{cl}\neq0$, $\alpha_\mathrm{rot}=\alpha_\mathrm{vib}=0$'; addend='__lns_acl';
    OFD = [cFldr,filesep,'OUTPUT_FILES_401758_2090_lns_noatt']; rsclFctr = 1.309279e+00; tShft = (2.25799996-2504.2-5.068+5.309)*1e-3; TIT = 'LNS, w/ $\alpha_\mathrm{cl}=\alpha_\mathrm{rot}=\alpha_\mathrm{vib}=0$'; addend='__lns_noatt';
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
for i = 1:NSda
  ANRGda(i) = trapz(t_stack{i}, p_stack{i}.^2);
end
for i = 1:NSsy
  ANRGsy(i) = trapz(Tsy(i, :), Psy(i, :).^2);
end
ANRGda = ANRGda.^0.5;
ANRGsy = ANRGsy.^0.5;
% scale acoustic NRG with distance
if(max(abs(Dsy-Dda))>1e-6)
  error('kek');
end
fitte = fit(Dsy', (ANRGsy./ANRGda), 'a*x+b');
NRG_scale_to_data = (Dsy*fitte.a+fitte.b)';
% Psy = Psy./NRG_scale_to_data; ANRGsy = ANRGsy./NRG_scale_to_data;
% disp(['[',mfilename,'] WARNING: SCALED SYNTHETICS BY ',sprintf('%.3f', fitte.a),' * R + ',sprintf('%.3f', fitte.b),'.']);

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
YLAB = ['$r \mathrm{ [m]} + ',sprintf('%.2f',distOverPTP),' \mathrm{ [m/Pa]} \times p'' \mathrm{ [Pa]}$'];
fig_tvd = figure('units','normalized','outerposition',[0,0,1,1]);
tightAxes = tight_subplot(1, 2, [0, 0.05], [0.12,0.08], [0.07, 0.02]);
LS = '-';
handlesToLines = [];
axes(tightAxes(1));
for i = 1:NSda
%   NAsy = 'data';
  colour = Cda(i, :);
  for ti=1:size(p_stack_save{i}, 1)
    plot(fac_t*(t_stack{i}+tShft), Dda(i) + distOverPTP * p_stack_save{i}(ti, :), 'linewidth', 2, 'color', min((colour+[1,1,1])/2,1),'linestyle', LS); hold on;
  end
  h = plot(fac_t*(t_stack{i}+tShft), Dda(i) + distOverPTP * p_stack{i}, 'color', colour, 'linestyle', LS, 'displayname', NAda);
  handlesToLines = [handlesToLines, h];
end
for i = 1:NSsy
  colour = Csy(i, :);
  h = plot(fac_t*Tsy(i, :), Dsy(i) + distOverPTP * Psy(i, :), 'displayname', NAsy, 'color', colour,'linestyle',LS); hold on;
  handlesToLines = [handlesToLines, h];
end
[pa, fa] = prefix_factor_values({amp});
hscale = plot(fac_t*[1,1]*tamp, Dsy(idsy)+distOverPTP*[-1,1]*amp, 'linewidth', 6, 'color', [0,1,0]*0.5, 'displayname', [sprintf('%.2f',fa*amp),'~',pa,'Pa scale']);
xlabel(XLAB);
ylabel(YLAB);
xlim(fac_t*XLIM);
xticks(fac_t*XTIC);
title(TIT);
legend([handlesToLines([1, NSda+1]), hscale], 'location', 'northwest');

if(do_save_plot)
  name_fig = ['synth_v_data__f=',sprintf('%04d', data_freq), 'Hz_p=', sprintf('%04.0f', data_pre_Pa),'Pa'];
  customSaveFig(fig_tvd, [figfolder, name_fig, addend], {'png'}, 9999);
end

% acoustic energy
XLIM = [1e-4, 1e-2];
axes(tightAxes(2));
semilogx(ANRGda, Dda, '.', 'color', Cda(1, :), 'markersize', 30, 'displayname', NAda); hold on;
semilogx(ANRGsy, Dsy, '.', 'color', Csy(1, :), 'markersize', 30, 'displayname', NAsy);
semilogx(ANRGsy./NRG_scale_to_data, Dsy, 'x', 'color', Csy(1, :), 'markersize', 15, 'displayname', [NAsy, ' $\times\left(',sprintf('%.1f', fitte.a),'r',sprintf('%+.1f', fitte.b),'\right)^{-1}$']);
yticklabels({});
legend('location', 'southwest');
xlim(XLIM);
xlabel('Acoustic Energy [Pa$\cdot$s]');

linkaxes(tightAxes, 'y');
shift = 0.1;
tightAxes(1).Position([3]) = tightAxes(1).Position([3]) + shift;
tightAxes(2).Position([1, 3]) = tightAxes(2).Position([1, 3]) + [1,-1]*shift;

