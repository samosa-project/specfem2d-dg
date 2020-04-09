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
  if(nargin<1)
    %error(['[',mfilename,', ERROR] Not enough input arguments. Needs at least an ''atmospheric_model_file''.']);
    atmospheric_model_file = input(['[',mfilename,'] Path to an ''atmospheric_model.dat'' file? > '], 's');
  end
%   addpath('/home/l.martire/Documents/work/mars/mars_is'); % customSaveFig
  if(not(exist('marker', 'var')))
    marker = '-';
  end
  if(not(exist('colour', 'var')))
    colour = 'k';
  end
  if(not(exist('atmalts', 'var')))
    atmalts = [];
  end
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
  
  FS = 18;
  FSlab = FS*0.9;
  colour_leftprop = [0.631, 0.0745, 0.180];
  colour_rightprop = [0,176,146]/255;
  T_in_degC = 1;
  
  % load
  [Z, RHO, T, C, ~, ~, ...
   ~, ~, KAPPA, MU, ~, ~, ~, W, CP, CV, GAMMA, FR, SVIB] = ...
   extract_atmos_model(atmospheric_model_file, 3, 0, 0);
  [datestr, posstr, secondaryinfo] = extract_atmos_model_setup(atmospheric_model_file);
  mainTitle = [posstr,', ',datestr];
  
  % truncate to maxalt
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
  
  % trick for nicer yticks
  if( (maxalt_provided && (maxalt>=1e3)) || (not(maxalt_provided) && max(Z)>=1e3))
    Z=Z/1e3;
    unitz = ['[km]']; unitgrad = 'km';
  else
    unitz = ['[m]']; unitgrad = 'm';
  end
  % adjust margin left of subplots accordingly
  marg_l = 0.04 + log10(max(Z))*0.01;
  
  if(T_in_degC)
    T = T-273.15; % [°C]
    unit_T = '$^\circ$C';
  else
    unit_T = 'K';
  end
  thresheqzero = 1e-6;
  
%   D = differentiation_matrix(Z, 0);
%   DZW = D*W;
%   DZW(1)=DZW(2); % Correction hack.
  DZW = gradient(W, Z);
  if(numel(Z)>2)
    DZW(1:2) = DZW(3); % Correction hack.
  end
  
  % Prepare labels.
  i = 1;
  Title_Subplot{i}='Density'; DName_Subplot{i} = ['density [kg/m$^3]']; XLab_Subplot{i} = ['$\rho$ [kg/m$^3$]']; i=i+1;
  Title_Subplot{i}='Sound Speed'; DName_Subplot{i} = ['$c$']; XLab_Subplot{i} = ['$c$, $c_\mathrm{eff}$ [m/s]']; i=i+1;
  Title_Subplot{i}='Attenuation [1/m]'; DName_Subplot{i} = ['$\log_{10}\left(\mathrm{attenuation}\right)$']; XLab_Subplot{i} = []; i=i+1;
  Title_Subplot{i}='Temperature'; DName_Subplot{i} = ['temperature [',unit_T,']']; XLab_Subplot{i} = ['$T$ [',unit_T,']']; i=i+1;
  Title_Subplot{i}='Wind'; DName_Subplot{i} = ['horizontal wind [m/s]']; XLab_Subplot{i} = ['horizontal wind $w_\mathrm{h}$ [m/s]']; i=i+1;
  Title_Subplot{i}='Wind Shear'; DName_Subplot{i} = ['wind shear [(m/s)/',unitgrad,']']; XLab_Subplot{i} = ['$\partial_zw_\mathrm{h}$ [(m/s)/',unitgrad,']']; i=i+1;
  
  % Start figure.
  fh=figure('units','normalized','outerposition',[0 0 1 1]);
  axxx = tight_subplot(2, 3, [0.16,0.03], [0.10,0.12], [marg_l, 0.036]); % [l c gaph gapw margin_bot margin_top marg_left marg_right]
  
%   axxx(1)=subplot(231);
  xlab = [];
  i=1;
  axes(axxx(i));
  myplot(RHO, Z, marker, colour, DName_Subplot{i}, 'sx');
  xlab = [xlab, xlabel(XLab_Subplot{i}, 'fontsize', FSlab)];
  ylabel(['$z$ ',unitz], 'fontsize', FSlab); % xlim([0.1*min(RHO),10*max(RHO)]);
  title(Title_Subplot{i});
  i=i+1;
  
%   axxx(2)=subplot(232);
  axes(axxx(i));
  myplot(C, Z, marker, colour, DName_Subplot{i}, 'p');
  xlab = [xlab, xlabel(XLab_Subplot{i}, 'fontsize', FSlab)]; hold on
  if(max(abs(W))>thresheqzero)
    % If wind==0, no need to plot effective sound speed.
    if(strcmp(LRorDWUW,'DWUW'))
      % plot downwind/upwind effective sound speed
      myplot(min(C+W,C-W), Z, marker, 'b', 'upwind $c_\mathrm{eff}$', 'p');
      myplot(max(C+W,C-W), Z, marker, 'r', 'downwind $c_\mathrm{eff}$', 'p');
    elseif(strcmp(LRorDWUW,'LR'))
      % plot left/right effective sound speed
      myplot(C+W, Z, marker, colour_rightprop, 'rightward $c_\mathrm{eff}$', 'p');
      myplot(C-W, Z, marker, colour_leftprop, 'leftward $c_\mathrm{eff}$', 'p');
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
%  htit = title({[posstr,', ',datestr], 'Sound Speed'});
  htit = title(Title_Subplot{i});
  i=i+1;
  h = annotation('textbox', [0.5, 1, 0, 0], 'string',mainTitle,'interpreter', 'latex','fitboxtotext', true,'fontsize', htit.FontSize, 'linestyle', 'none','horizontalalignment', 'center');
  
%   axxx(3)=subplot(233);
  axes(axxx(i));
  Npts = 100;
  if(numel(RHO)>Npts*5)
    lightened_sel = floor(linspace(1,numel(RHO),Npts*5));
  else
    lightened_sel = 1:numel(RHO);
  end
  freq = logspace(-2,1,Npts);
  [a_cl, a_rot, a_vib] = atmospheric_attenuation(freq, RHO, T, C, MU, KAPPA,CP, CV,GAMMA, FR, SVIB);
%   [a_cl, a_rot, a_vib] = atmospheric_attenuation(freq, RHO, T, C, MU, KAPPA,CP, CV,GAMMA, [], []);
  a_tot = a_cl+a_rot+a_vib;
  %figure();subplot(1,2,1);contourf(freq,Z,real(a_tot));colorbar;subplot(1,2,2);contourf(freq,Z,real(a_tot));colorbar;pause
  a_tot = abs(a_tot);
%   a_tot = a_vib;
  if(abs(max(max(a_tot)))==0)
    to_plot = a_tot;
    lab = '$\alpha_\mathrm{cl}+\alpha_\mathrm{rot}+\alpha_\mathrm{vib}$';
    cust_cmap = 1;
  else
    to_plot = log10(a_tot);
    lab = '$\log_{10}\left(\alpha_\mathrm{cl}+\alpha_\mathrm{rot}+\alpha_\mathrm{vib}\right)$';
    cust_cmap = 0;
  end
  contourf(freq, Z, to_plot, 15-1, 'ShowText', 'off', 'edgecolor', 'none', 'displayname', DName_Subplot{i});
  yticklabels([]);
  set(gca, 'xscale', 'log');
  xlab = [xlab, xlabel(['frequency [Hz]'], 'fontsize', FSlab)];
  shading interp;
  h = colorbar;
  ylabel(h, [lab], 'rotation', 90, 'interpreter', 'latex');
  if(cust_cmap)
    colormaps_custom([0,1], [1*[1,1,1];0*[1,1,1]], 1);
    caxis([0,1]);
  else
    colormaps_fromPython('YlOrRd', 1);
  end
  title(Title_Subplot{i});
  i=i+1;
  
%   axxx(4)=subplot(234);
  axes(axxx(i));
  myplot(T, Z, marker, colour, DName_Subplot{i}, 'p');
  xlab = [xlab, xlabel(XLab_Subplot{i}, 'fontsize', FSlab)];
  ylabel(['$z$ ',unitz], 'fontsize', FSlab); % LIMX=WN; xlim([1.1*min(LIMX)-0.1*max(LIMX),1.1*max(LIMX)-0.1*min(LIMX)]);
  addatmosphericseparationlines(T, atmalts);
  if(max(abs(T-mean(T)))>thresheqzero)
    % If T==cst, no need to adjust.
    forcexlimminmax(T);
  end
  title(Title_Subplot{i});
  i=i+1;
  
%   axxx(5)=subplot(235);
  axes(axxx(i));
  myplot(W, Z, marker, colour, DName_Subplot{i}, 'p');
  xlab = [xlab, xlabel(XLab_Subplot{i}, 'fontsize', FSlab)];
  line([0,0],[Z(1),Z(end)],'linestyle',':','color','k','linewidth',2);
  yticklabels([]);
  addatmosphericseparationlines(W, atmalts);
  if(max(abs(W-mean(W)))>thresheqzero)
    % If W==cst, no need to adjust.
    forcexlimminmax(W);
  end
  title(Title_Subplot{i});
  i=i+1;
  
%   axxx(6)=subplot(236);
  axes(axxx(i));
  myplot(DZW, Z, marker, colour, DName_Subplot{i}, 'p');
  xlab = [xlab, xlabel(XLab_Subplot{i}, 'fontsize', FSlab)];
  line([0,0],[Z(1),Z(end)],'linestyle',':','color','k','linewidth',2);
  yticklabels([]);
  if(max(abs(DZW))>thresheqzero)
    % If DZW==0 (W==cst), no need to adjust.
    forcexlimminmax(DZW);
  end
  title(Title_Subplot{i});
  i=i+1;
  
  % Post-treatment.
  linkaxes(axxx,'y');
  linkprop(axxx,{'ytick','yscale'});
  set(axxx, 'yTickMode', 'auto'); % reset to auto in case user zooms
  
  prettyAxes(fh);
  
  for i=1:numel(xlab)
    xlab(i).Units = 'normalized';
    xlab(i).Position = xlab(i).Position + [0, 0.04, 0];
  end
  
  if(plot_nosave0_save1)
    customSaveFig(fh,atmospheric_model_file);
    disp(['[] Use the following to save in other formats.']);
    disp(['[] TEX:   customSaveFig(gcf,''',atmospheric_model_file,''',{''tex''});']);
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