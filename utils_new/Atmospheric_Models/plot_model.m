% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [] = plot_model(output_file, marker, colour, atmalts)
  format compact;
  set(0, 'DefaultLineLineWidth', 3); set(0, 'DefaultLineMarkerSize', 8);
  set(0, 'defaultTextFontSize', 26); set(0, 'defaultAxesFontSize', 28);
  set(0, 'DefaultTextInterpreter', 'latex');
  set(0, 'DefaultLegendInterpreter', 'latex');
  
  [Z, RHO, T, C, ~, ~, ...
   ~, NSQ, ~, ~, ~, ~, ~, W, ~, ~, GAMMA] = ...
   extract_data(output_file, 3, 0, 0);

  [datestr, posstr, ~, d, s, ~, ~, F107A, F107, AP] = extract_data_setup(output_file);
  apf107str=strcat("F10.7 avg. = ", sprintf("%.1f",F107A), ", F10.7 = ", sprintf("%.1f",F107), ", AP = ", sprintf("%.1f",AP));
  
  T=T-273.15;
  
%   D = differentiation_matrix(Z, 0);
%   DZW = D*W;
%   DZW(1)=DZW(2); % Correction hack.
  DZW=gradient(W,Z);
  DZW(1:2)=DZW(3); % Correction hack.
  
  figure();
  
  ax(1)=subplot(231);
  myplot(RHO, Z, marker, colour, datestr, 'sx');
  xlabel('$\rho$ (kg/m$^3$)');
  ylabel('$z$ (m)'); % xlim([0.1*min(RHO),10*max(RHO)]);
  
  ax(2)=subplot(232);
  myplot(C, Z, marker, colour, '$c$', 'p'); xlabel('$c$ (m/s)'); hold on
  myplot(min(C+W,C-W), Z, marker, 'b', 'upwind $c_e$', 'p');
  myplot(max(C+W,C-W), Z, marker, 'r', 'downwind $c_e$', 'p');
  legend('Location', 'best');
  yticklabels([]);
  
  title({[posstr,', ',datestr],apf107str,''});
  
  ax(3)=subplot(233);
  myplot(NSQ, Z, marker, colour, datestr, 'p');
  xlabel('$N^2$ (rad$^2$/s$^2$)'); % LIMX=NSQ; xlim([1.1*min(LIMX)-0.1*max(LIMX),1.1*max(LIMX)-0.1*min(LIMX)]);
  yticklabels([]);
  
  ax(4)=subplot(234);
  myplot(T, Z, marker, colour, datestr, 'p');
  xlabel('$T$ ($^\circ$C)');
  ylabel('$z$ (m)'); % LIMX=WN; xlim([1.1*min(LIMX)-0.1*max(LIMX),1.1*max(LIMX)-0.1*min(LIMX)]);
  for i=1:length(atmalts)
    line([min(T), max(T)],[atmalts(i),atmalts(i)],'linestyle',':','color','k','linewidth',2);
  end
  
  ax(5)=subplot(235);
  myplot(W, Z, marker, colour, datestr, 'p');
  xlabel('$w$ (m/s)');
  line([0,0],[Z(1),Z(end)],'linestyle',':','color','k','linewidth',2);
  yticklabels([]);
  
  ax(6)=subplot(236);
  myplot(DZW, Z, marker, colour, datestr, 'p');
  xlabel('$\partial_zw$ (1/s)');
  line([0,0],[Z(1),Z(end)],'linestyle',':','color','k','linewidth',2);
  yticklabels([]);
  
  linkaxes(ax,'y');
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