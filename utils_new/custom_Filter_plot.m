% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Plots the Bode diagram of a filter.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [] = custom_Filter_plot(freq,filt)
  set(0, 'DefaultTextInterpreter', 'latex');
  set(0, 'DefaultLegendInterpreter', 'latex');
  set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
  set(0, 'defaultTextFontSize', 14); set(0, 'defaultAxesFontSize', 14);
  
  figure();
  
  axx(1)=subplot(211);
  magnitude=10*log10(abs(filt));
  semilogx(freq,magnitude);
  ylim([floor(1.1*min(magnitude)), ceil(1.1*max(magnitude))]);
  title('Bode Diagram');
  ylabel('magnitude [dB]');
  grid on;
  set(gca,'ticklabelinterpreter','latex');
  set(gca,'tickdir','both');
  
  axx(2)=subplot(212);
  semilogx(freq,angle(filt)/pi*180);
  ylabel('phase [deg]');
  xlabel('frequency [Hz]');
  grid on;
  set(gca,'ticklabelinterpreter','latex');
  set(gca,'tickdir','both');
  
  linkaxes(axx,'x');
  xlim([min(abs(freq)),max(freq)]);
end

