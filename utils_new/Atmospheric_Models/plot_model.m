% Author:        Léo Martire.
% Description:   Plots an atmospheric model main quantities.
% Notes:         Needs:
%                a) .m scripts and functions (if not alongside this
%                   script, recover via Léo):
%                  1) extract_atmos_model.m
%
% Usage:
%   plot_model(atmospheric_model_file, marker, colour, atmalts, LRorDWUW)
% with:
%   TODO.
% yields:
%   TODO.

function [] = plot_model(atmospheric_model_file, marker, colour, atmalts, LRorDWUW, plot_nosave0_save1)
  if(nargin<4)
    error(['[',mfilename,', ERROR] Not enough input arguments. Needs ''atmospheric_model_file, marker, colour, atmalts'', with atmalts possibly [].']);
  end
  addpath('/home/l.martire/Documents/work/mars/mars_is'); % customSaveFig
  
  if(not(exist('LRorDWUW')))
    % not found, make default: LR
    LRorDWUW='LR';
  else
    if(strcmp(LRorDWUW,'DWUW'))
      % found alright
    elseif(strcmp(LRorDWUW,'LR'))
      % found alright
    else
      % found but not alright, make default: LR
      LRorDWUW='LR';
    end
  end
  if(not(exist('plot_nosave0_save1')))
    plot_nosave0_save1 = 1;
  end
  
%   format compact;
%   set(0, 'DefaultLineLineWidth', 3); set(0, 'DefaultLineMarkerSize', 8);
%   set(0, 'defaultTextFontSize', 14); set(0, 'defaultAxesFontSize', 14);
%   set(0, 'DefaultTextInterpreter', 'latex');
%   set(0, 'DefaultLegendInterpreter', 'latex');
  
  [Z, RHO, T, C, ~, ~, ...
   ~, NSQ, ~, ~, ~, ~, ~, W, ~, ~, ~] = ...
   extract_atmos_model(atmospheric_model_file, 3, 0, 0);

  [datestr, posstr, secondaryinfo] = extract_atmos_model_setup(atmospheric_model_file);
  
  T=T-273.15;
  thresheqzero=1e-6;
  
%   D = differentiation_matrix(Z, 0);
%   DZW = D*W;
%   DZW(1)=DZW(2); % Correction hack.
  DZW=gradient(W,Z);
  DZW(1:2)=DZW(3); % Correction hack.
  
  fh=figure('units','normalized','outerposition',[0 0 1 1]);
  
  axxx = tight_subplot(2, 3, [0.08,0.01], [0.08,0.1], [0.05, 0.01]); % [l c gaph gapw margin_bot margin_top marg_left marg_right]
  
%   axxx(1)=subplot(231);
  i=1;
  axes(axxx(i)); i=i+1;
  myplot(RHO, Z, marker, colour, datestr, 'sx');
  xlabel('$\rho$ [kg/m$^3$]');
  ylabel('$z$ [m]'); % xlim([0.1*min(RHO),10*max(RHO)]);
  
%   axxx(2)=subplot(232);
  axes(axxx(i)); i=i+1;
  myplot(C, Z, marker, colour, '$c$', 'p'); xlabel('$c$ [m/s]'); hold on
  if(max(abs(W))>thresheqzero)
    % If wind==0, no need to plot effective sound speed.
    if(strcmp(LRorDWUW,'DWUW'))
      % plot downwind/upwind effective sound speed
      myplot(min(C+W,C-W), Z, marker, 'b', 'upwind $c_e$', 'p');
      myplot(max(C+W,C-W), Z, marker, 'r', 'downwind $c_e$', 'p');
    elseif(strcmp(LRorDWUW,'LR'))
      % plot left/right effective sound speed
      myplot(C+W, Z, marker, [0.9290    0.6940    0.1250], 'rightward $c_e$', 'p');
      myplot(C-W, Z, marker, [0.6350    0.0780    0.1840], 'leftward $c_e$', 'p');
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
  
  title({[posstr,', ',datestr],secondaryinfo,''});
  
%   axxx(3)=subplot(233);
  axes(axxx(i)); i=i+1;
  myplot(NSQ, Z, marker, colour, datestr, 'p');
  xlabel('$N^2$ [rad$^2$/s$^2$]'); % LIMX=NSQ; xlim([1.1*min(LIMX)-0.1*max(LIMX),1.1*max(LIMX)-0.1*min(LIMX)]);
  yticklabels([]);
  if(max(abs(NSQ-mean(NSQ)))>thresheqzero)
    % If NSQ==cst, no need to adjust.
    forcexlimminmax(NSQ);
  else
    xlim(max(NSQ)*[.9,1.1]);
  end
  
%   axxx(4)=subplot(234);
  axes(axxx(i)); i=i+1;
  myplot(T, Z, marker, colour, datestr, 'p');
  xlabel('$T$ [$^\circ$C]');
  ylabel('$z$ [m]'); % LIMX=WN; xlim([1.1*min(LIMX)-0.1*max(LIMX),1.1*max(LIMX)-0.1*min(LIMX)]);
  addatmosphericseparationlines(T, atmalts);
  if(max(abs(T-mean(T)))>thresheqzero)
    % If T==cst, no need to adjust.
    forcexlimminmax(T);
  end
  
%   axxx(5)=subplot(235);
  axes(axxx(i)); i=i+1;
  myplot(W, Z, marker, colour, datestr, 'p');
  xlabel('$w$ [m/s]');
  line([0,0],[Z(1),Z(end)],'linestyle',':','color','k','linewidth',2);
  yticklabels([]);
  addatmosphericseparationlines(W, atmalts);
  if(max(abs(W-mean(W)))>thresheqzero)
    % If W==cst, no need to adjust.
    forcexlimminmax(W);
  end
  
%   axxx(6)=subplot(236);
  axes(axxx(i)); i=i+1;
  myplot(DZW, Z, marker, colour, datestr, 'p');
  xlabel('$\partial_zw$ [1/s]');
  line([0,0],[Z(1),Z(end)],'linestyle',':','color','k','linewidth',2);
  yticklabels([]);
  if(max(abs(DZW))>thresheqzero)
    % If DZW==0 (W==cst), no need to adjust.
    forcexlimminmax(DZW);
  end
  
  linkaxes(axxx,'y');
  
%   spl=split(atmospheric_model_file,'.');
%   spl{length(spl)+1}='jpg';
%   spl=join(spl,'.');
%   saveas(gcf,spl{1}, 'jpg');
%   disp(['[',mfilename,'] Plot of model saved to ''',spl{1},'''.']);
  figure(fh.Number);
  if(plot_nosave0_save1)
    customSaveFig(atmospheric_model_file);
  end
end

function myplot(XDATA, YDATA, MARKER, COLOUR, DISPLAYNAME, TYPE)
  if(strcmp(TYPE,'p'))
    plot(XDATA, YDATA, MARKER, 'Color', COLOUR, 'DisplayName', DISPLAYNAME); hold on;
  elseif(strcmp(TYPE,'sx'))
    semilogx(XDATA, YDATA, MARKER, 'Color', COLOUR, 'DisplayName', DISPLAYNAME); hold on;
  end
  ylim([min(YDATA), max(YDATA)]);
  set(gca, 'TickLabelInterpreter','latex');
  set(gca,'TickDir','both');
  grid on;
end

function addatmosphericseparationlines(v, atmalts)
  for i=1:length(atmalts)
    line([min(v), max(v)],[atmalts(i),atmalts(i)],'linestyle',':','color','k','linewidth',1,'HandleVisibility','off');
  end
end

function forcexlimminmax(v)
  xlim([min(v), max(v)]);
end