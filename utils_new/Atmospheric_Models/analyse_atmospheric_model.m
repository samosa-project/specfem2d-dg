% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

clear all;
close all;
clc;
set(0, 'DefaultLineLineWidth', 1.5); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 6); % Default at 6.
%set(0, 'defaultTextFontSize', 20);
set(0, 'defaultAxesFontSize', 12); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCD Model Loading.                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interpolate = 0;
interp_delta = 4000; % Interpolation step (m).
save_plots = 0;
% If one has to cut data, choose maximum altitude here. Put Inf if all altitudes are to be considered.
% maxalt=Inf;
maxalt=150e3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Regularisation method. %%%%%%
% method='integrate';
method='bruteforce';
% method='metaheuristic';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set datafile.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp(['[INPUT] Path (relative or absolute) to datafile to be loaded? You''re in ''', pwd ,'''.']); DATAFILE = input(' > ');
% DATAFILE = 'Mars/Mars_atmosphere_models/MARS_MCD_AGW_model_INSIGHT_1_Oct_2016_12h.dat'; headerlines=1;
% DATAFILE = 'Mars/Mars_atmos_models_LEO/msise_Sumatra_model_DG.dat'; headerlines=1; % WARNING: DATA FORMAT IS NOT THE SAME AS THE ONE IN MARS MODELS!
% DATAFILE = "Earth/wrapper/msisehwm_model_output"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/0_100000_111_0.00000_0.00000_0_0_0.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/0_100000_111_45.00000_0.00000_0_0_0.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/0_100000_111_-45.00000_0.00000_0_0_0.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/0_100000_11_0.00000_0.00000_0_0_0.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/0_100000_21_0.00000_0.00000_0_0_0.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/0_100000_21_-45.00000_0.00000_0_0_0.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/0_60000_111_-45.00000_0.00000_0_0_0.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/tests/0_300000_51_66.56306_0.00000_0_356_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_ATMOSPHERIC_SELECTED/0_300000_301_66.56306_0.00000_0_173_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_ATMOSPHERIC_SELECTED/0_300000_301_66.56306_0.00000_0_356_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_ATMOSPHERIC_SELECTED/0_300000_301_45.00000_0.00000_0_173_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_ATMOSPHERIC_SELECTED/0_300000_301_45.00000_0.00000_0_356_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_ATMOSPHERIC_SELECTED/0_300000_301_23.43694_0.00000_0_173_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_ATMOSPHERIC_SELECTED/0_300000_301_23.43694_0.00000_0_356_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_ATMOSPHERIC_SELECTED/0_300000_301_0.00000_0.00000_0_173_0.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_ATMOSPHERIC_SELECTED/0_300000_301_0.00000_0.00000_0_173_43200.00000_0.00000"; headerlines=3;

DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/stratospheric/0_200000_301_66.56306_0.00000_0_173_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/stratospheric/0_200000_301_66.56306_0.00000_0_173_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/stratospheric/0_200000_301_66.56306_0.00000_0_356_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/stratospheric/0_200000_301_45.00000_0.00000_0_173_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/stratospheric/0_200000_301_45.00000_0.00000_0_356_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/stratospheric/0_200000_301_23.43694_0.00000_0_173_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/stratospheric/0_200000_301_23.43694_0.00000_0_356_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/stratospheric/0_200000_301_0.00000_0.00000_0_173_0.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/stratospheric/stratospheric/0_200000_301_0.00000_0.00000_0_173_43200.00000_0.00000"; headerlines=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Loading.");
[datestr, posstr, ~, ~, ~, ~, ~] = extract_data_setup(DATAFILE);
[Z, RHO, TEMP, SOUNDSPEED, P, LOCALPRESSURESCALE, ...
 G, NBVSQ, KAPPA, VISCMU, MUVOL, WNORTH, WEAST, W, CP, CV, GAMMA] = ...
 extract_data(DATAFILE, headerlines, interpolate, interp_delta);
N = NBVSQ .^ 0.5 / (2 * pi);

plot_model(DATAFILE, 230, '-', 'k', 'no');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deduce other parameters.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NBV = NBVSQ.^0.5 / (2 * pi); % Unused.
% NA = NCUTA / (2 * pi); % Unused.
% coef_A = 1;
% coef_B = interp_SVIB ./ (2 * pi * interp_VIBRATTENUATION);
% coef_C = - (2 * pi * interp_VIBRATTENUATION) .^ (- 2);
% interp_tau_sig = (- coef_B + (coef_B .^ 2.0 - 4.0 * coef_A * coef_C) .^ 0.5) / (2.0 * coef_A);
% interp_tau_eps = interp_tau_sig + coef_B;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots.                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure();
% semilogx(interp_tau_eps, interp_ALTITUDE, 'o', interp_tau_sig, interp_ALTITUDE, '.', abs(interp_tau_eps-interp_tau_sig), interp_ALTITUDE);
% xlabel('relaxation time (s)'); ylabel('altitude (m)');
% legend('\tau_\epsilon','\tau_\sigma','|\tau_\epsilon - \tau_\sigma|', 'Location', 'best');
% title('Relaxation Times for Martian Atmosphere');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hydrostatic Equilibrium Treatment.                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eventually cut data.        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_maxalt=find(abs(Z-maxalt)==min((abs(Z-maxalt))), 1, 'last');
Z=Z(1:ind_maxalt);
RHO=RHO(1:ind_maxalt);
TEMP=TEMP(1:ind_maxalt);
P=P(1:ind_maxalt);
G=G(1:ind_maxalt);
N=N(1:ind_maxalt);
KAPPA=KAPPA(1:ind_maxalt);
VISCMU=VISCMU(1:ind_maxalt);
WEAST=WEAST(1:ind_maxalt);
WNORTH=WNORTH(1:ind_maxalt);
W=W(1:ind_maxalt);
nz=length(Z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(" ");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check Evaluation            %
% Coefficients.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tit_plus={posstr, datestr};
D=differentiation_matrix(Z, 0);

hydrostatic_ratio = @(D, P, RHO, G) (D * P) ./ (-RHO .* G);
east_Richardson = Richardson_number(WEAST, D, N);
north_Richardson = Richardson_number(WNORTH, D, N);
proj_Richardson = Richardson_number(W, D, N);
HR=hydrostatic_ratio(D, P, RHO, G);

if(any(proj_Richardson<1))
  disp(["Projected wind's Richardson number is < 1 somewhere. Unstability can occur."]);
end
if(any(proj_Richardson<0.25))
  disp(["[WARNING] Projected wind's Richardson number is < 0.25 somewhere. Unstability will probably occur."]);
end

disp(strcat(['Maximum relative gap to hydrostatic state (via ratio): ', num2str(max(abs(HR-1))*100), '%.']));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots.                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
semilogx(abs(D * (VISCMU .* (D * W))), Z, 'r'); hold on;
semilogx(abs(D * (KAPPA .* (D * TEMP) + W .* VISCMU .* (D * W))), Z, 'g');
semilogx(abs(D * P + RHO .* G), Z, 'b');
tmp_valmat=cell2mat(get(get(gca, 'children'), 'XData')); tmp_plot_maxval=max(max(tmp_valmat)); tmp_valmat(tmp_valmat==0)=Inf; tmp_plot_minval=min(min(tmp_valmat));
xlim([0.5*tmp_plot_minval, 2*tmp_plot_maxval]);
xlabel({'amplitude of hydrostatic unbalance terms', '(projected wind)'}); ylabel('altitude (m)');
legend('$\left|\partial_z\left(\mu\partial_zw\right)\right|$', '$\left|\partial_z\left(\kappa\partial_zT+w\mu\partial_zw\right)\right|$', '$\left|\partial_zp + \rho g_z\right|$', 'Location', 'best');
title(tit_plus);
if save_plots == 1
  saveas(gcf, strcat(DATAFILE,'__unbalance_terms.png'));
end

% figure();
% plot(WEAST, Z, WNORTH, Z);
% % xlim([0.5 * min([east_Richardson; north_Richardson]), 2 * max([east_Richardson; north_Richardson])]);
% xlabel('wind (m/s)'); ylabel('altitude (m)');
% legend('eastward wind', 'northward wind', 'Location', 'NorthWest');
% title(tit_plus);
% if save_plots == 1
%   saveas(gcf, strcat(DATAFILE,'__winds.png'));
% end

figure();
% plot(HR(2:end), Z(2:end), ones(size(Z)), Z, 'k:');
plot(HR, Z, ones(size(Z)), Z, 'k:');
xlabel('$\partial_zp/(-\rho g)$'); ylabel('altitude (m)');
% legend('hydrostatic ratio - 1', 'Location', 'best');
title(tit_plus);

figure();
semilogx(east_Richardson, Z, north_Richardson, Z, proj_Richardson, Z, ones(size(Z)), Z, 'k:', 0.25*ones(size(Z)), Z, 'k');
xlim([0.5 * min([east_Richardson; north_Richardson]), 2 * max([east_Richardson; north_Richardson])]);
xlabel('Richardson number'); ylabel('altitude (m)');
legend('eastward Richardson number', 'northward Richardson number', 'projected wind Richardson number', 'Location', 'NorthWest');
title(tit_plus);
if save_plots == 1
  saveas(gcf, strcat(DATAFILE,'__richardson.png'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(" ");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularise model.           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(["Trying regularisation."]);
modify_atmospheric_model; % See script.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%