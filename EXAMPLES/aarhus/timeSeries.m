clear all;
close all;
clc;

thisFolder = [regexprep(mfilename('fullpath'),mfilename,'')];
figfolder = [thisFolder, filesep, 'FIGURES', filesep];
addpath(genpath(thisFolder));

do_save = 0;

timestamp = '182608'; nrun = 397; pre_Pa = 20*100; freq = 2090;

switch(freq)
  case 2090
    OFD = [thisFolder,filesep,'OUTPUT_FILES_401255_2090']; rescaleFactor = 36931.5408386775;
  otherwise
    error('kek');
end

istattab = 1:5;
typeDisplay = 2;
distChoice = 1; % x, |x|, z, d
doGeometricAttenuation = 0;
subsample = 1;
subsample_wantedDt = 1e-6;
distOverPTP = 0.75;

% synthetics
[Tsy, Psy, Dsy, Nsy, ~] = gji2020_loadSomeSynthetics(OFD, istattab, typeDisplay, 'BXZ', distChoice, doGeometricAttenuation, rescaleFactor, subsample, subsample_wantedDt);
nstat_synth = size(Tsy,1);
Dsy = Dsy.vals;
Csy = repmat([1,0,0], nstat_synth, 1);

% data
% disp(['PRODUCING SIMULATED DATA']);
% Tda = Tsy;
% Pda = Psy + 1e-6*(rand(size(Psy))-0.5)*2;
% Tda = [Tda, Tda+0.1e-3+max(Tda(:))];
% Pda = [Pda, Pda]*1e4;
% save([thisFolder,filesep,'TMP_simulated_data.mat'], 'Tda', 'Pda'); pause
% load([thisFolder,filesep,'TMP_simulated_data.mat']);
% NSda = size(Tda, 1);
% Dda = Dsy+0.01;
[Tda, Pda, Dda] = Load_Data_ATN(timestamp, nrun, pre_Pa/100, freq);
Pda = Pda';
NSda = size(Pda, 1);
Tda = repmat(Tda, NSda, 1);
Cda = repmat([0,0,0], NSda, 1);

% disp(['[',mfilename,'] Ratio between data amplitude and synthetic amplitude (i.e. recommended scaling of synthetics): ',sprintf('%.6e', range(Pda(:))/range(Psy(:))),'.']);

% stacc
tstacctime = {};
tstacctime{1} = [2.504205, 4.004375, 5.50361, 7.001905, 8.499265, 9.995695, 11.49119, 12.985755, 14.479395, 15.972095];
tstacctime{2} = [2.50493, 4.0051, 5.504335, 7.00263, 8.499995, 9.996425, 11.49192, 12.986485, 14.48012, 15.972825];
tstacctime{3} = [2.506775, 4.006945, 5.50618, 7.004475, 8.50184, 9.998265, 11.49376, 12.98833, 14.481965, 15.97467];
tstacctime{4} = [2.51049, 4.01066, 5.509895, 7.008185, 8.505555, 10.00198, 11.49748, 12.992045, 14.485675, 15.97838];
tstacctime{5} = [2.514205, 4.01437, 5.5136, 7.0119, 8.50926, 10.00569, 11.501185, 12.995745, 14.48939, 15.982085];
kind = 'bp';
%   fcut = [1150, 2923];
% fcut = [1050, 3495];
fcut = [958.7, 4085];
order = 2;
[tstack, pselsave, pstacc] = get_staccs(Tda, Pda, tstacctime, kind, fcut, order);
disp(['[',mfilename,'] Ratio between data amplitude and synthetic amplitude (i.e. recommended scaling of synthetics): ',sprintf('%.6e', range(pstacc{1})/range(Psy(1,:))),'.']);
% figure();
% i=4;
% for ti=1:numtstac
%   plot(tstack{i}, pselsave{i}(ti, :), 'linewidth', 2, 'color', [1,0.75,0.75]); hold on;
% end
% plot(tstack{i}, pstacc{i}, 'color', [1,0,0]*0.25);
% % select data
% % TODO eventually stack
% twindow = [-Inf, Inf];
% % twindow = [18e-3, 32e-3];
% NTda = []; NPda = [];
% for i = 1:NSda
%   selt_da = (Tda(i, :)>=min(twindow)) & (Tda(i, :)<=max(twindow));
%   NTda(i, :) = Tda(i, selt_da);
%   NPda(i, :) = Pda(i, selt_da);
% end
% Tda = NTda; Pda = NPda;
% Tda = Tda - min(Tda(:));

% Adjust timing w.r.t. synthetics.
% tshift = -(2.415-4.325)*1e-3;
% tshift = (12.2606-2514.205)*1e-3;
tshift = (1.74760004-2504.2)*1e-3;
% tshift = 0;

% prepare a scale to print
idsy = 4;
tamp = 1e-3;
amp = range(Psy(idsy, :));
amp = round(amp,abs(ceil(log10(amp)))+3);

% TIT = 'kek';
% [pre_t, fac_t] = prefix_factor_values({Tsy(:), Tda(:)});
[pre_t, fac_t] = prefix_factor_values({Tsy(:)});
XLAB = ['time [',pre_t,'s]'];
XLIM = [0,16e-3]; XTIC = (0:1:16)*1e-3;
YLAB = ['$r \mathrm{ [m]} + ',sprintf('%.2f',distOverPTP),' \mathrm{ [m/Pa]} \times p'' \mathrm{ [Pa]}$'];
fig_tvd = figure('units','normalized','outerposition',[0,0,1,1]);
tightAxes = tight_subplot(1, 1, [0,0], [0.12,0.05], [0.07, 0.02]);
LS = '-';
handlesToLines = [];
for i = 1:NSda
  dname = 'data';
  colour = Cda(i, :);
%   h = plot(fac_t*(Tda(i, :)+tshift), Dda(i) + distOverPTP * Pda(i, :), 'displayname', dname, 'color', colour,'linestyle',LS); hold on;
  
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
hscale = plot(fac_t*[1,1]*tamp, Dsy(idsy)+distOverPTP*[-1,1]*amp, 'linewidth', 5, 'color', [0,1,0]*0.5, 'displayname', [sprintf('%.2f',fa*amp),'~',pa,'Pa scale']);
xlabel(XLAB);
ylabel(YLAB);
xlim(fac_t*XLIM);
xticks(fac_t*XTIC);
legend([handlesToLines([1, NSda+1]), hscale], 'location', 'best');

if(do_save)
  name_fig = ['synth_v_data__f=',sprintf('%04d', freq), 'Hz_p=', sprintf('%04.0f', pre_Pa),'Pa'];
  customSaveFig(fig_tvd, [figfolder, name_fig], {'fig', 'eps', 'png', 'tex'}, 9999);
end
