clear all;
close all;
clc;

thisFolder = [regexprep(mfilename('fullpath'),mfilename,'')];
figfolder = [thisFolder, filesep, 'FIGURES', filesep];
addpath(genpath(thisFolder));

do_save = 0;

nrun = 322; pre_Pa = 20*100; freq = 2090;

switch(freq)
  case 2090
    OFD = [thisFolder,filesep,'OUTPUT_FILES_401255_2090']; rescaleFactor = 1.122997e+04;
  otherwise
    error('kek');
end

istattab = 1:5;
typeDisplay = 2;
distChoice = 1; % x, |x|, z, d
doGeometricAttenuation = 0;
subsample = 0;
subsample_wantedDt = 1;
distOverPTP = 3;

% synthetics
[Tsy, Psy, Dsy, Nsy, ~] = gji2020_loadSomeSynthetics(OFD, istattab, typeDisplay, 'BXZ', distChoice, doGeometricAttenuation, rescaleFactor, subsample, subsample_wantedDt);
nstat_synth = size(Tsy,1);
Dsy = Dsy.vals;
Csy = repmat([1,0,0], nstat_synth, 1);

% data
disp(['PRODUCING SIMULATED DATA']);
% Tda = Tsy;
% Pda = Psy + 1e-6*(rand(size(Psy))-0.5)*2;
% Tda = [Tda, Tda+0.1e-3+max(Tda(:))];
% Pda = [Pda, Pda]*1e4;
% save([thisFolder,filesep,'TMP_simulated_data.mat'], 'Tda', 'Pda'); pause
load([thisFolder,filesep,'TMP_simulated_data.mat']);
NSda = size(Tda, 1);
Dda = Dsy+0.01;
% [Tda, Pda, Dda] = Load_Data_ATN(nrun, pre_Pa/100, freq);
% NSda = size(dataMIC, 1);
% Tda = repmat(Tda, NSda, 1);
Cda = repmat([0,0,0], NSda, 1);

disp(['[',mfilename,'] Ratio between data amplitude and synthetic amplitude (i.e. recommended scaling of synthetics): ',sprintf('%.6e', range(Pda(:))/range(Psy(:))),'.']);

% select data
% TODO eventually stack
twindow = [18e-3, 32e-3];
NTda = []; NPda = [];
for i = 1:NSda
  selt_da = (Tda(i, :)>=min(twindow)) & (Tda(i, :)<=max(twindow));
  NTda(i, :) = Tda(i, selt_da);
  NPda(i, :) = Pda(i, selt_da);
end
Tda = NTda; Pda = NPda;
Tda = Tda - min(Tda(:));

% Adjust timing w.r.t. synthetics.
tshift = -(2.415-4.325)*1e-3;
% tshift = 0;

% prepare a scale to print
idsy = 4;
tamp = 6e-3;
amp = range(Psy(idsy, :));
amp = round(amp,abs(ceil(log10(amp)))+3);

% TIT = 'kek';
[pre_t, fac_t] = prefix_factor_values({Tsy(:), Tda(:)});
XLAB = ['time [',pre_t,'s]'];
XLIM = [0,16e-3]; XTIC = (0:1:16)*1e-3;
YLAB = ['$r \mathrm{ [m]} + ',scientific_latex_notation(distOverPTP,0),' \mathrm{ [m/Pa]} \times p'' \mathrm{ [Pa]}$'];
fig_tvd = figure('units','normalized','outerposition',[0,0,1,1]);
tightAxes = tight_subplot(1, 1, [0,0], [0.12,0.05], [0.07, 0.02]);
LS = '-';
handlesToLines = [];
for i = 1:NSda
  dname = 'data';
  colour = Cda(i, :);
  h = plot(fac_t*(Tda(i, :)+tshift), Dda(i) + distOverPTP * Pda(i, :), 'displayname', dname, 'color', colour,'linestyle',LS); hold on;
  handlesToLines = [handlesToLines, h];
end
for i = 1:nstat_synth
  dname = 'synthetics';
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