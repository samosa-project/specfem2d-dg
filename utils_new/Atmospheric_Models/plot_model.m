% Author:        Léo Martire.
% Description:   Plots an atmospheric model main quantities.
% Notes:         Needs:
%                a) .m scripts and functions (if not alongside this
%                   script, recover via Léo):
%                  1) extract_atmos_model.m
%
% Usage:
%   plot_model(atmospheric_model_file, marker, colour, atmalts, maxalt, LRorDWUW, plot_nosave0_save1)
% with:
%   TODO.
% yields:
%   TODO.

function [] = plot_model(atmospheric_model_file, marker, colour, atmalts, maxalt, LRorDWUW, plot_nosave0_save1)
  if(nargin<4)
    error(['[',mfilename,', ERROR] Not enough input arguments. Needs ''atmospheric_model_file, marker, colour, atmalts'', with atmalts possibly [].']);
  end
  addpath('/home/l.martire/Documents/work/mars/mars_is'); % customSaveFig
  
  FS = 18;
  FSlab = FS*0.9;
  
  set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
  set(0, 'defaultTextFontSize', FS); set(0, 'defaultAxesFontSize', FS);
  set(0, 'DefaultTextInterpreter', 'latex');
  set(0, 'DefaultLegendInterpreter', 'latex');
  
  leftprop = [0.631, 0.0745, 0.180];
  rightprop = [0,176,146]/255;
  
  if(not(exist('maxalt','var')))
    maxalt_provided = 0;
  else
    maxalt_provided = 1;
  end
  
  if(not(exist('LRorDWUW','var')))
    % not found, make default: LR
    LRorDWUW='LR';
  else
    LRorDWUW = upper(LRorDWUW); % safeguard
    if(strcmp(LRorDWUW,'DWUW'))
      % found alright
    elseif(strcmp(LRorDWUW,'LR'))
      % found alright
    else
      % found but not alright, make default: LR
      LRorDWUW='LR';
    end
  end
  if(not(exist('plot_nosave0_save1','var')))
    plot_nosave0_save1 = 1;
  end
  
%   format compact;
%   set(0, 'DefaultLineLineWidth', 3); set(0, 'DefaultLineMarkerSize', 8);
%   set(0, 'defaultTextFontSize', 14); set(0, 'defaultAxesFontSize', 14);
%   set(0, 'DefaultTextInterpreter', 'latex');
%   set(0, 'DefaultLegendInterpreter', 'latex');
  
  [Z, RHO, T, C, ~, ~, ...
   ~, ~, KAPPA, MU, ~, ~, ~, W, CP, CV, GAMMA, FR, SVIB] = ...
   extract_atmos_model(atmospheric_model_file, 3, 0, 0);
  
  if(maxalt_provided)
    sel = (Z<=maxalt);
    Z = Z(sel);
    RHO = RHO(sel);
    T = T(sel);
    C = C(sel);
    KAPPA = KAPPA(sel);
    MU = MU(sel);
%     NSQ = NSQ(sel);
    W = W(sel);
    CP = CP(sel);
    CV = CV(sel);
    GAMMA = GAMMA(sel);
    if(not(isempty(FR)))
      FR = FR(sel);
      SVIB = SVIB(sel);
    end
  end
  
  [datestr, posstr, secondaryinfo] = extract_atmos_model_setup(atmospheric_model_file);
  
  T = T-273.15; % [°C]
  thresheqzero=1e-6;
  
%   D = differentiation_matrix(Z, 0);
%   DZW = D*W;
%   DZW(1)=DZW(2); % Correction hack.
  DZW = gradient(W,Z);
  DZW(1:2) = DZW(3); % Correction hack.
  
  fh=figure('units','normalized','outerposition',[0 0 1 1]);
  
  axxx = tight_subplot(2, 3, [0.12,0.03], [0.08,0.11], [0.05, 0.03]); % [l c gaph gapw margin_bot margin_top marg_left marg_right]
  
%   axxx(1)=subplot(231);
  i=1;
  axes(axxx(i)); i=i+1;
  myplot(RHO, Z, marker, colour, datestr, 'sx');
  xlabel('$\rho$ [kg/m$^3$]', 'fontsize', FSlab);
  ylabel('$z$ [m]', 'fontsize', FSlab); % xlim([0.1*min(RHO),10*max(RHO)]);
  title({'Density'});
  
%   axxx(2)=subplot(232);
  axes(axxx(i)); i=i+1;
  myplot(C, Z, marker, colour, '$c$', 'p');
  xlabel('$c$, $c_\mathrm{eff}$ [m/s]', 'fontsize', FSlab); hold on
  if(max(abs(W))>thresheqzero)
    % If wind==0, no need to plot effective sound speed.
    if(strcmp(LRorDWUW,'DWUW'))
      % plot downwind/upwind effective sound speed
      myplot(min(C+W,C-W), Z, marker, 'b', 'upwind $c_\mathrm{eff}$', 'p');
      myplot(max(C+W,C-W), Z, marker, 'r', 'downwind $c_\mathrm{eff}$', 'p');
    elseif(strcmp(LRorDWUW,'LR'))
      % plot left/right effective sound speed
      myplot(C+W, Z, marker, rightprop, 'rightward $c_\mathrm{eff}$', 'p');
      myplot(C-W, Z, marker, leftprop, 'leftward $c_\mathrm{eff}$', 'p');
    else
      error('kek');
    end
    legend('Location', 'best');
  end
  yticklabels([]);
  addatmosphericseparationlines([min(min(C+W),min(C-W)),max(max(C+W),max(C-W))], atmalts);
  if(max(abs(W))>thresheqzero)
    % If wind==0, no need to adjust.
    forcexlimminmax([min(min(C+W),min(C-W)),max(max(C+W),max(C-W))]);
  end
  
%   title({[posstr,', ',datestr],secondaryinfo,'', 'Sound Speed'});
  title({[posstr,', ',datestr,secondaryinfo], '', 'Sound Speed'});
  
%   axxx(3)=subplot(233);
  axes(axxx(i)); i=i+1;
%   myplot(NSQ, Z, marker, colour, datestr, 'p');
%   xlabel('$N^2$ [rad$^2$/s$^2$]'); % LIMX=NSQ; xlim([1.1*min(LIMX)-0.1*max(LIMX),1.1*max(LIMX)-0.1*min(LIMX)]);
%   yticklabels([]);
%   if(max(abs(NSQ-mean(NSQ)))>thresheqzero)
%     % If NSQ==cst, no need to adjust.
%     forcexlimminmax(NSQ);
%   else
%     xlim(max(NSQ)*[.9,1.1]);
%   end
  Npts = 100;
  if(numel(RHO)>Npts*5)
    lightened_sel = floor(linspace(1,numel(RHO),Npts*5));
  else
    lightened_sel = 1:numel(RHO);
  end
  freq = logspace(-2,1,Npts);
  [a_cl, a_rot, a_vib] = atmospheric_attenuation(freq, RHO, C, MU, KAPPA,CP, CV,GAMMA, FR, SVIB);
%   [a_cl, a_rot, a_vib] = atmospheric_attenuation(freq, RHO, C, MU, KAPPA,CP, CV,GAMMA, [], []);
  a_tot = a_cl+a_rot+a_vib;
  a_tot = a_vib;
  if(abs(max(max(a_tot)))==0)
    to_plot = a_tot;
    lab = '$\alpha_\mathrm{cl}+\alpha_\mathrm{rot}+\alpha_\mathrm{vib}$';
    cust_cmap = 1;
  else
    to_plot = log10(a_tot);
    lab = '$\log_{10}\left(\alpha_\mathrm{cl}+\alpha_\mathrm{rot}+\alpha_\mathrm{vib}\right)$';
    cust_cmap = 0;
  end
  contourf(freq, Z, to_plot, 15-1, 'ShowText', 'off', 'edgecolor', 'none');
  yticklabels([]);
  set(gca, 'xscale', 'log');
  xlabel(['frequency [Hz]'], 'fontsize', FSlab);
  shading interp;
  h = colorbar;
  ylabel(h, [lab], 'rotation', 90, 'interpreter', 'latex');
  if(cust_cmap)
    colormaps_custom([0,1], [1*[1,1,1];0*[1,1,1]], 1);
    caxis([0,1]);
  else
    colormaps_fromPython('YlOrRd', 1);
  end
  title({'Attenuation'});
  
%   axxx(4)=subplot(234);
  axes(axxx(i)); i=i+1;
  myplot(T, Z, marker, colour, datestr, 'p');
  xlabel('$T$ [$^\circ$C]', 'fontsize', FSlab);
  ylabel('$z$ [m]', 'fontsize', FSlab); % LIMX=WN; xlim([1.1*min(LIMX)-0.1*max(LIMX),1.1*max(LIMX)-0.1*min(LIMX)]);
  addatmosphericseparationlines(T, atmalts);
  if(max(abs(T-mean(T)))>thresheqzero)
    % If T==cst, no need to adjust.
    forcexlimminmax(T);
  end
  title({'Temperature'});
  
%   axxx(5)=subplot(235);
  axes(axxx(i)); i=i+1;
  myplot(W, Z, marker, colour, datestr, 'p');
  xlabel('horizontal wind, $w_\mathrm{h}$ [m/s]', 'fontsize', FSlab);
  line([0,0],[Z(1),Z(end)],'linestyle',':','color','k','linewidth',2);
  yticklabels([]);
  addatmosphericseparationlines(W, atmalts);
  if(max(abs(W-mean(W)))>thresheqzero)
    % If W==cst, no need to adjust.
    forcexlimminmax(W);
  end
  title({'Wind'});
  
%   axxx(6)=subplot(236);
  axes(axxx(i)); i=i+1;
  myplot(DZW, Z, marker, colour, datestr, 'p');
  xlabel('$\partial_zw_\mathrm{h}$ [1/s]', 'fontsize', FSlab);
  line([0,0],[Z(1),Z(end)],'linestyle',':','color','k','linewidth',2);
  yticklabels([]);
  if(max(abs(DZW))>thresheqzero)
    % If DZW==0 (W==cst), no need to adjust.
    forcexlimminmax(DZW);
  end
  title({'Wind Shear'});
  
  linkaxes(axxx,'y');
  linkprop(axxx,{'ytick','yscale'});
%   set(gca,'yscale','log')
  
%   % find closest multiple of 3 under log10(max(Z))
%   Nticks = 6;
%   nearestPow = log10(max(Z)) - mod(log10(max(Z)),3);
% %   yticks(linspace(min(Z),max(Z),7));
%   tentativeYticks = 10^nearestPow * floor(linspace(min(Z),max(Z),Nticks)/(10^nearestPow));
%   definitiveYticks = min(tentativeYticks) + min(diff(tentativeYticks))*linspace(0,Nticks-1,Nticks); % correct to have same dz between ticks
%   if(definitiveYticks(1)<min(Z))
%     definitiveYticks(1)=min(Z);
%   end
%   yticks(definitiveYticks);
%   set(axxx, 'yTickMode', 'auto', 'yTickLabelMode', 'auto'); % reset to auto in case user zooms
  set(axxx, 'yTickMode', 'auto'); % reset to auto in case user zooms
%   if(maxalt_provided)
%     ylim([0,maxalt]);
%   end
  
%   spl=split(atmospheric_model_file,'.');
%   spl{length(spl)+1}='jpg';
%   spl=join(spl,'.');
%   saveas(gcf,spl{1}, 'jpg');
%   disp(['[',mfilename,'] Plot of model saved to ''',spl{1},'''.']);
%   figure(fh.Number);
  prettyAxes(fh);
  if(plot_nosave0_save1)
    customSaveFig(fh,atmospheric_model_file);
  end
end

function myplot(XDATA, YDATA, MARKER, COLOUR, DISPLAYNAME, TYPE)
  if(strcmp(TYPE,'p'))
    plot(XDATA, YDATA, MARKER, 'Color', COLOUR, 'DisplayName', DISPLAYNAME); hold on;
  elseif(strcmp(TYPE,'sx'))
    semilogx(XDATA, YDATA, MARKER, 'Color', COLOUR, 'DisplayName', DISPLAYNAME); hold on;
  end
  ylim([min(YDATA), max(YDATA)]);
%   set(gca, 'TickLabelInterpreter','latex');
%   set(gca,'TickDir','both');
%   grid on;
end

function addatmosphericseparationlines(v, atmalts)
  for i=1:length(atmalts)
    line([min(v), max(v)],[atmalts(i),atmalts(i)],'linestyle',':','color','k','linewidth',1,'HandleVisibility','off');
  end
end

function forcexlimminmax(v)
  xlim([min(v), max(v)]);
end