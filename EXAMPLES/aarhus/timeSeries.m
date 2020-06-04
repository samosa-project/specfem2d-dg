clear all;
close all;
clc;

thisFolder = [regexprep(mfilename('fullpath'),mfilename,'')];
figfolder = [thisFolder, filesep, 'FIGURES', filesep];
addpath(genpath(thisFolder));

% plot parameters
distOverPTP = 0.75;
do_save_plot = 0;

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
%     OFD = [thisFolder,filesep,'OUTPUT_FILES_401255_2090_old']; rescaleFactor = 36931.5321885834; tshift = (1.74760004-2504.2)*1e-3; TIT = 'LNS, w/ $\alpha_\mathrm{cl}\neq0$, $\alpha_\mathrm{rot}=\alpha_\mathrm{vib}=0$';
%     OFD = [thisFolder,filesep,'OUTPUT_FILES_401736_2090']; rescaleFactor = 6.242541e+04; tshift = (2.25099991-2504.2)*1e-3; TIT = 'FNS, w/ $\alpha_\mathrm{cl}\neq0$, $\alpha_\mathrm{rot}=0$, $\alpha_\mathrm{vib}\neq0$';
%     OFD = [thisFolder,filesep,'OUTPUT_FILES_401741_2090_noatt']; rescaleFactor = 5.191609e+04; tshift = (2.25099991-2504.2)*1e-3; TIT = 'FNS, w/ $\alpha_\mathrm{cl}=\alpha_\mathrm{rot}=\alpha_\mathrm{vib}=0$';
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
[Tsy, Psy, Dsy, Nsy, ~] = gji2020_loadSomeSynthetics(OFD, istattab, typeDisplay, 'BXZ', distChoice, doGeometricAttenuation, rescaleFactor, subsample, subsample_wantedDt);
nstat_synth = size(Tsy,1);
Dsy = Dsy.vals;
Csy = repmat([1,0,0], nstat_synth, 1);

% load data
warning off;
[Tda, Pda, Dda] = Load_Data_ATN(data_timestamp, data_nrun, data_pre_Pa/100, data_freq);
warning on;
Pda = Pda';
NSda = size(Pda, 1);
Tda = repmat(Tda, NSda, 1);
Cda = repmat([0,0,0], NSda, 1);

% stacc data
tstacctime = {};
tstacctime{1} = [2.504205, 4.004375, 5.50361, 7.001905, 8.499265, 9.995695, 11.49119, 12.985755, 14.479395, 15.972095];
tstacctime{2} = [2.50493, 4.0051, 5.504335, 7.00263, 8.499995, 9.996425, 11.49192, 12.986485, 14.48012, 15.972825];
tstacctime{3} = [2.506775, 4.006945, 5.50618, 7.004475, 8.50184, 9.998265, 11.49376, 12.98833, 14.481965, 15.97467];
tstacctime{4} = [2.51049, 4.01066, 5.509895, 7.008185, 8.505555, 10.00198, 11.49748, 12.992045, 14.485675, 15.97838];
tstacctime{5} = [2.514205, 4.01437, 5.5136, 7.0119, 8.50926, 10.00569, 11.501185, 12.995745, 14.48939, 15.982085];
[tstack, pselsave, pstacc] = get_staccs(Tda, Pda, tstacctime, data_filter_kind, data_filter_fcut, data_filter_order);
ratio = range(pstacc{1})/range(Psy(1,:));
if(abs(1-ratio)>1e-6)
  disp(['[',mfilename,'] Ratio between data amplitude and synthetic amplitude (i.e. recommended scaling of synthetics): ',sprintf('%.6e', ratio),'.']);
  pause;
end

% Adjust timing w.r.t. synthetics.
% tshift = 0;

% prepare a scale to print
idsy = 4;
tamp = 1e-3;
amp = range(Psy(idsy, :));
amp = round(amp,abs(ceil(log10(amp)))+3);

% plot
[pre_t, fac_t] = prefix_factor_values({Tsy(:)});
XLAB = ['time [',pre_t,'s]'];
XLIM = [0,16e-3]; XTIC = (0:1:16)*1e-3;
YLAB = ['$r \mathrm{ [m]} + ',sprintf('%.2f',distOverPTP),' \mathrm{ [m/Pa]} \times p'' \mathrm{ [Pa]}$'];
fig_tvd = figure('units','normalized','outerposition',[0,0,1,1]);
tightAxes = tight_subplot(1, 1, [0,0], [0.12,0.08], [0.07, 0.02]);
LS = '-';
handlesToLines = [];
for i = 1:NSda
  dname = 'data';
  colour = Cda(i, :);
  for ti=1:size(pselsave{i}, 1)
    plot(fac_t*(tstack{i}+tshift), Dda(i) + distOverPTP * pselsave{i}(ti, :), 'linewidth', 2, 'color', min((colour+[1,1,1])/2,1),'linestyle', LS); hold on;
  end
  dname = 'filtered data stack';
  h = plot(fac_t*(tstack{i}+tshift), Dda(i) + distOverPTP * pstacc{i}, 'color', colour, 'linestyle', LS, 'displayname', dname);
  handlesToLines = [handlesToLines, h];
end
for i = 1:nstat_synth
  dname = 'non-filtered synthetics';
  colour = Csy(i, :);
  h = plot(fac_t*Tsy(i, :), Dsy(i) + distOverPTP * Psy(i, :), 'displayname', dname, 'color', colour,'linestyle',LS); hold on;
  handlesToLines = [handlesToLines, h];
end
[pa, fa] = prefix_factor_values({amp});
hscale = plot(fac_t*[1,1]*tamp, Dsy(idsy)+distOverPTP*[-1,1]*amp, 'linewidth', 6, 'color', [0,1,0]*0.5, 'displayname', [sprintf('%.2f',fa*amp),'~',pa,'Pa scale']);
xlabel(XLAB);
ylabel(YLAB);
xlim(fac_t*XLIM);
xticks(fac_t*XTIC);
title(TIT);
legend([handlesToLines([1, NSda+1]), hscale], 'location', 'best');

if(do_save_plot)
  name_fig = ['synth_v_data__f=',sprintf('%04d', data_freq), 'Hz_p=', sprintf('%04.0f', data_pre_Pa),'Pa'];
  customSaveFig(fig_tvd, [figfolder, name_fig], {'fig', 'eps', 'png', 'tex'}, 9999);
end
