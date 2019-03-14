% Author:        Léo Martire.
% Description:   Plots an atmospheric model absorption parameters.
% Notes:         Needs:
%                a) .m scripts and functions (if not alongside this
%                   script, recover via Léo):
%                  1) frsvib2tausigtaueps.m
%                  2) extract_atmos_model.m
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function [fh] = plot_model_frsvib(modelfile, zmax_interest)
  addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/Atmospheric_Models/tools');
  
  [Z, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, FR, SVIB] = extract_atmos_model(modelfile, 3, 0, 0);
  
  [TSIG, TEPS] = frsvib2tausigtaueps(FR, SVIB);
  sel = find(Z <= zmax_interest);

  axx = [];
  fh = figure('units','normalized','outerposition',[0 0 1 1]);
  subplot(131); axx = [axx, gca()];
  semilogx(FR(sel),Z(sel),'.');
%     xlim(rad2deg*[0,2*pi]);
  ylabel(['$z$ [m]']);
  xlabel(['$f_{r}$ [Hz]']);
%     xticks(rad2deg*2*pi*linspace(0,1,9));
  grid on; box on; set(gca,'ticklabelinterpreter','latex'); set(gca,'tickdir','both');

  subplot(132); axx = [axx, gca()];
  plot(SVIB(sel),Z(sel),'.');
  xlabel(['$S_{vib}$ [m/s]']);
  yticklabels([]);
  grid on; box on; set(gca,'ticklabelinterpreter','latex'); set(gca,'tickdir','both');

  subplot(133); axx = [axx, gca()];
  semilogx(TSIG(sel),Z(sel),'.','displayname','$\tau_\sigma$'); hold on;
  semilogx(TEPS(sel),Z(sel),'.','displayname','$\tau_\epsilon$'); hold on;
  xlabel(['$\tau_{\sigma,\epsilon}$ [s]']);
  yticklabels([]);
  legend('location','best');
  grid on; box on; set(gca,'ticklabelinterpreter','latex'); set(gca,'tickdir','both');
  linkaxes(axx,'y');
end

