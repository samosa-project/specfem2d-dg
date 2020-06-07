clear all;
close all;
clc;

cFldr = [regexprep(mfilename('fullpath'),mfilename,'')];
figfolder = [cFldr, filesep, 'FIGURES', filesep];
addpath(genpath(cFldr));

% plot parameters
distOverPTP = 2;
do_save_plot = 1;
p0 = 20*100; f0 = 2090; % TODO: read from parfile

% synthetics parameters
OFD_w = [cFldr,filesep,'OUTPUT_FILES_402141_2090_lns_allatt']; OFD_wo = [cFldr,filesep,'OUTPUT_FILES_402830_2090_lns_noatt']; addend = '__lns';
% OFD_w = [cFldr, filesep, 'OUTPUT_FILES_402036_2090_fns_allatt']; OFD_wo = [cFldr,filesep,'OUTPUT_FILES_402879_2090_fns_noatt']; addend = '__fns';

istattab = 1:5; % should not change
typeDisplay = 2; % should not change
distChoice = 1; % should not change
rsclFctr = 1;
doGeometricAttenuation = 0;
subsample = 1;
subsample_wantedDt = 1e-6;

% load synthetics
[Tsy_w, Psy_w, Dsy_w, Nsy_w, ~] = gji2020_loadSomeSynthetics(OFD_w, istattab, typeDisplay, 'BXZ', distChoice, doGeometricAttenuation, rsclFctr, subsample, subsample_wantedDt);
[Tsy_wo, Psy_wo, Dsy_wo, Nsy_wo, ~] = gji2020_loadSomeSynthetics(OFD_wo, istattab, typeDisplay, 'BXZ', distChoice, doGeometricAttenuation, rsclFctr, subsample, subsample_wantedDt);
NSsy = size(Tsy_w,1);
Dsy = Dsy_w.vals';
Csy = repmat([1,0,0], NSsy, 1);

% make sure same scale
pstf_w = [OFD_w, filesep,'plot_source_time_function.txt'];
pstf_wo = [OFD_wo, filesep,'plot_source_time_function.txt'];
stf_w = importdata(pstf_w);
stf_wo = importdata(pstf_wo);
ratio = stf_wo(:,2)./stf_w(:,2);
wo_over_w = mean(ratio(not(isnan(ratio))));
Psy_w = Psy_w*wo_over_w;

% define integration bounds
i = 1;
twin_w(i, :) = [-Inf, 2.387]; twin_wo(i, :) = [-Inf, 2.394]; i = i+1;
twin_w(i, :) = [-Inf, 3.341]; twin_wo(i, :) = [-Inf, 3.351]; i = i+1;
twin_w(i, :) = [-Inf, 5.182]; twin_wo(i, :) = [-Inf, 5.196]; i = i+1;
twin_w(i, :) = [-Inf, 8.868]; twin_wo(i, :) = [-Inf, 8.887]; i = i+1;
twin_w(i, :) = [-Inf, 12.458]; twin_wo(i, :) = [-Inf, 12.528]; i = i+1;
% Compute Acoustic NRG.
ANRG_w = zeros(NSsy, 1);
ANRG_wo = zeros(NSsy, 1);
for i = 1:NSsy
%   window = [0,16]*1e-3; % [ms]
  window_w = twin_w(i,:)*1e-3; % [ms]
  sel_w = (Tsy_w(i, :)>=min(window_w)) & (Tsy_w(i, :)<=max(window_w));
  ANRG_w(i) = trapz(Tsy_w(i, sel_w), Psy_w(i, sel_w).^2);
  window_wo = twin_wo(i,:)*1e-3; % [ms]
  sel_wo = (Tsy_wo(i, :)>=min(window_wo)) & (Tsy_wo(i, :)<=max(window_wo));
  ANRG_wo(i) = trapz(Tsy_wo(i, sel_wo), Psy_wo(i, sel_wo).^2);
end
% ANRG_w = ANRG_w.^0.5;
% ANRG_wo = ANRG_wo.^0.5;

% compute amplitude of first peak
Afp_w = zeros(NSsy, 1);
Afp_wo = zeros(NSsy, 1);
for i = 1:NSsy
  [~, Afp_w(i), ~] = findFirstPeak(Tsy_w(i,:), max(Psy_w(i,:),0), 0.3);
  [~, Afp_wo(i), ~] = findFirstPeak(Tsy_wo(i,:), max(Psy_wo(i,:),0), 0.3);
end

% plot
fac_t = 1e3; pre_t = 'm';
YLAB = ['$r \mathrm{ [m]} + ',sprintf('%.2f',distOverPTP),' \mathrm{ [m/Pa]} \times p'' \mathrm{ [Pa]}$'];
Cwo = [0,0,1]*0.66;
Cw = min(Cwo+[1,1,1]*0.5,1);
fig_tvd = figure('units','normalized','outerposition',[0,0,1,1]);
tightAxes = tight_subplot(1, 1, [0, 0.05], [0.12,0.08], [0.07, 0.02]);
LS = '-';
handlesToLines = [];
% axes(tightAxes(1));
for i = 1:NSsy
  dnam = '';
  handlesToLines = [handlesToLines, plot(fac_t*Tsy_wo(i,:), Dsy(i) + distOverPTP * Psy_wo(i,:), 'color', Cwo, 'linestyle', LS, 'displayname', 'without attenuation')]; hold on;
  handlesToLines = [handlesToLines, plot(fac_t*Tsy_w(i,:), Dsy(i) + distOverPTP * Psy_w(i,:), 'color', Cw, 'linestyle', LS, 'displayname', 'with attenuation')]; hold on;
end
xlim(fac_t*[min(Tsy_w(:)), max(Tsy_w(:))]);
legend(handlesToLines([1,2]), 'location', 'northwest');
xlabel(['time [',pre_t,'s]']); ylabel(YLAB);
xticks(fac_t*(0:1:16)*1e-3);
if(do_save_plot)
  name_fig_tvd = ['synth_w_v_wo__f=',sprintf('%04d', f0), 'Hz_p=', sprintf('%04.0f', p0),'Pa'];
  customSaveFig(fig_tvd, [figfolder, name_fig_tvd, addend], {'fig', 'eps', 'png', 'tex'}, 9999);
end

% fit nrg
alpha_th = 0.319137;
fitt_nrg = 'c - 2*a*x - 2*b*log(x)'; fitfunc_nrg = fittype(fitt_nrg); op = fitoptions(fitfunc_nrg); fits_nrg={}; colfit={}; i = 1;
op.Lower = [0, 0, -Inf]; op.Upper = [0, Inf, Inf]; fits_nrg{i} = fit(Dsy', log(ANRG_wo), fitfunc_nrg, op); colfit{i}=Cwo; i=i+1; % fit without to get b
op.Lower = [0, fits_nrg{1}.b, -Inf]; op.Upper = [Inf, fits_nrg{1}.b, Inf]; fits_nrg{i} = fit(Dsy', log(ANRG_w), fitfunc_nrg, op); colfit{i}=Cw; i=i+1; % fit without to get b
op.Lower = [alpha_th, fits_nrg{1}.b, -Inf]; op.Upper = [alpha_th, fits_nrg{1}.b, Inf]; fits_nrg_th = fit(Dsy', log(ANRG_w), fitfunc_nrg, op); colfit_nrg_th=[1,0,0]; i=i+1; % theoretical
% fit amp
fitt_amp = 'c * exp(-a*x) / x.^b'; fitfunc_amp = fittype(fitt_amp); op = fitoptions(fitfunc_amp); fits_amp={}; colfit_amp={}; i = 1;
op.Lower = [0, 0, -Inf]; op.Upper = [0, Inf, Inf]; fits_amp{i} = fit(Dsy', Afp_wo, fitfunc_amp, op); colfit_amp{i}=Cwo; i=i+1; % fit without to get b
% op.Lower = [0, fits_nrg{1}.b, -Inf]; op.Upper = [0, fits_nrg{1}.b, Inf]; fits_amp{i} = fit(Dsy', Afp_wo, fitfunc_amp, op); colfit_amp{i}=Cwo; i=i+1; % fit with b from nrg
op.Lower = [0, fits_amp{1}.b, -Inf]; op.Upper = [Inf, fits_amp{1}.b, Inf]; fits_amp{i} = fit(Dsy', Afp_w, fitfunc_amp, op); colfit{i}=Cw; i=i+1; % fit without to get b
op.Lower = [alpha_th, fits_amp{1}.b, -Inf]; op.Upper = [alpha_th, fits_amp{1}.b, Inf]; fits_amp_th = fit(Dsy', Afp_w, fitfunc_amp, op); colfit_amp_th=[1,0,0]; i=i+1; % theoretical
% plot nrg
ms = 250;
fig_anrg = figure('units','normalized','outerposition',[0,0,1,1]);
tightAxes = tight_subplot(2, 1, [0.05, 0.05], [0.12,0.08], [0.09, 0.01]);
axes(tightAxes(1));
hfit=[];
x = linspace(min(Dsy), max(Dsy), 40);
hfit=[hfit,semilogy(x, exp(fits_nrg_th(x)), 'displayname', ['$\alpha=\alpha_\mathrm{th}=', sprintf('%.6f',alpha_th), '$, $\beta=',sprintf('%.2f', fits_nrg_th.b),'$'], 'color', colfit_nrg_th)]; hold on;
hpts = [];
hpts = [hpts, scatter(Dsy, ANRG_wo, ms, 'o', 'markerfacecolor', Cwo, 'markeredgecolor', 'none', 'displayname', 'without attenuation')]; hold on;
hpts = [hpts, scatter(Dsy, ANRG_w, ms, 'o', 'markerfacecolor', Cw, 'markeredgecolor', 'none', 'displayname', 'with attenuation')]; hold on;
for i=1:numel(fits_nrg)
  if(i>1)
    cfintcurfit = confint(fits_nrg{i});
    dnam = ['$\alpha=', sprintf('%.3f',fits_nrg{i}.a), '\pm',sprintf('%.3f',mean(abs(cfintcurfit(:,1)-fits_nrg{i}.a))),'$, $\beta=',sprintf('%.2f', fits_nrg{i}.b),'$'];
    semilogy(x, exp(fits_nrg{i}.c-2*max(cfintcurfit(1,1),0)*x-2*fits_nrg{i}.b*log(x)), 'displayname', dnam, 'linestyle', ':', 'color', colfit{i}); hold on;
    semilogy(x, exp(fits_nrg{i}.c-2*max(cfintcurfit(2,1),0)*x-2*fits_nrg{i}.b*log(x)), 'displayname', dnam, 'linestyle', ':', 'color', colfit{i}); hold on;
  else
    dnam = ['$\alpha=', sprintf('%.3f',fits_nrg{i}.a), '$, $\beta=',sprintf('%.2f', fits_nrg{i}.b),'$'];
  end
  hfit = [hfit, semilogy(x, exp(fits_nrg{i}(x)), 'displayname', dnam, 'color', colfit{i})]; hold on;
end
ylabel(['acoustic energy [Pa$^2\cdot$s]']);
legend([hpts, hfit],'location', 'eastoutside');
axes(tightAxes(2));
hfit=[];
x = linspace(min(Dsy), max(Dsy), 40);
hfit=[hfit,semilogy(x, fits_amp_th(x), 'displayname', ['$\alpha=\alpha_\mathrm{th}=', sprintf('%.6f',alpha_th), '$, $\beta=',sprintf('%.2f', fits_amp_th.b),'$'], 'color', colfit_amp_th)]; hold on;
hpts = [];
hpts = [hpts, scatter(Dsy, Afp_wo, ms, '^', 'markerfacecolor', Cwo, 'markeredgecolor', 'none', 'displayname', 'without attenuation')]; hold on;
hpts = [hpts, scatter(Dsy, Afp_w, ms, '^', 'markerfacecolor', Cw, 'markeredgecolor', 'none', 'displayname', 'with attenuation')]; hold on;
for i=1:numel(fits_amp)
%   dnam = ['$\alpha=', scientific_latex_notation(fits{i}.a,2,0), '$, $\beta=',sprintf('%.2f', fits{i}.b),'$, $c_1=',sprintf('%.2f', fits{i}.c),'$'];
%   dnam = ['$\alpha=', scientific_latex_notation(fits_amp{i}.a,2,0), '$, $\beta=',sprintf('%.2f', fits_amp{i}.b),'$'];
%   hfit=[hfit,semilogy(x, fits_amp{i}(x), 'displayname', dnam, 'color', colfit{i})]; hold on;
  if(i>1)
    cfintcurfit = confint(fits_amp{i});
    dnam = ['$\alpha=', sprintf('%.3f',fits_amp{i}.a), '\pm',sprintf('%.3f',mean(abs(cfintcurfit(:,1)-fits_amp{i}.a))),'$, $\beta=',sprintf('%.2f', fits_amp{i}.b),'$'];
    semilogy(x, fits_amp{i}.c*exp(-max(cfintcurfit(1,1),0)*x)./x.^(fits_amp{i}.b), 'displayname', dnam, 'linestyle', ':', 'color', colfit{i}); hold on;
    semilogy(x, fits_amp{i}.c*exp(-max(cfintcurfit(2,1),0)*x)./x.^(fits_amp{i}.b), 'displayname', dnam, 'linestyle', ':', 'color', colfit{i}); hold on;
  else
    dnam = ['$\alpha=', sprintf('%.3f',fits_amp{i}.a), '$, $\beta=',sprintf('%.2f', fits_amp{i}.b),'$'];
  end
  hfit = [hfit, semilogy(x, fits_amp{i}(x), 'displayname', dnam, 'color', colfit{i})]; hold on;
end
xlabel(['distance to the source $r$ [m]']);
ylabel(['amplitude [Pa]']);
legend([hpts, hfit],'location', 'eastoutside');
set(tightAxes, 'yscale', 'log');
set(tightAxes(1), 'ylim', [1e-7, 1e-5]);
set(tightAxes(1), 'xticklabel', {});
set(tightAxes(2), 'ylim', [1e-2, 2e-1]);
linkaxes(tightAxes, 'x');
xlim([0.25, 3.5]);
linkprop(tightAxes,'xtick');
tightAxes(1).Position([3]) = tightAxes(2).Position([3]);

if(do_save_plot)
  name_fig_nrg = ['synth_w_v_wo__f=',sprintf('%04d', f0), 'Hz_p=', sprintf('%04.0f', p0),'Pa__nrg'];
%   customSaveFig(fig_tvd, [figfolder, name_fig, addend], {'png'}, 9999);
  customSaveFig(fig_anrg, [figfolder, name_fig_nrg, addend], {'fig', 'eps', 'png', 'tex'}, 9999);
end


