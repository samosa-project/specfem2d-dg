% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Loads an atmospheric model. Analyses hydrostatic balance
%                and other stability coefficients. Tries regularisations.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

clear all;
close all;
clc;
set(0, 'DefaultLineLineWidth', 2); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 8); % Default at 6.
set(0, 'defaultTextFontSize', 12);
set(0, 'defaultAxesFontSize', 12); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
richardson_advice='Should be >1, or at least >0.25 to prevent instability.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Loading.                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Base parameters.
interpolate = 0; % Set to 1 to activate interpolation.
interp_delta = 4000; % Interpolation step (m).
save_plots = 0; % Set to 1 to save plots.
maxalt=Inf; % If one has to cut data, choose maximum altitude here. Put Inf if all altitudes are to be considered.
% maxalt=140e3; % If one has to cut data, choose maximum altitude here. Put Inf if all altitudes are to be considered.

% Regularisation method.
method='bruteforce_rho'; % Bruteforce $\rho = -\partial_z{P} / g_z$.
% method='bruteforce_rho_log'; % Bruteforce $\rho = -(\partial_z{\log_10(P)} . P) / g_z$ (rewriting, given sensibly the same results as those of 'bruteforce_rho').
% method='integrate'; % Use an iterative scheme (using the inverse of the differentiation matrix) to obtain $P = \int (\rho g_z) dz$ (in fact, use the alternative variable $log_10(P)$ to do so more efficiently).
% method='metaheuristic'; % Use metaheuristics optimisation methods to minimise $\int_{z_{min}}^{z_{max}} |1-\partial_z{P}/(\rho g_z)| dz$, ie the integral of the absolute gap to hydrostatic ratio.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set DATAFILE.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp(['[INPUT] Path (relative or absolute) to datafile to be loaded? You''re in ''', pwd ,'''.']); DATAFILE = input(' > ', 's'); disp(['[INPUT] Number of header lines of datafile to be loaded? Default is 3.']); headerlines = input(' > ');
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/atmospheric/stratospheric/0_150000_1501_66.56306_0.00000_0_173_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/Atmospheric_Models/Earth/msiseecmwf_tests/msiseecmwf_2016_149__-37.00000_280.00000_0_70000_4000_0.00000"; headerlines=3;
DATAFILE=input(['[',mfilename,'] Input atmospheric model file to use > '],'s'); headerlines=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['[',mfilename,'] Loading.']);
disp(['[',mfilename,'] > File:']);
disp(strcat("  '",DATAFILE, "'."));
[datestr, posstr, ~, ~, ~, ~, ~] = extract_atmos_model_setup(DATAFILE);
[Z, RHO, TEMP, SOUNDSPEED, P, LOCALPRESSURESCALE, ...
 G, NBVSQ, KAPPA, VISCMU, MUVOL, WNORTH, WEAST, W, CP, CV, GAMMA] = ...
 extract_atmos_model(DATAFILE, headerlines, interpolate, interp_delta);
% NBVSQ: rad^2/s^2
% Nf = NBVSQ .^ 0.5 / (2 * pi); % 1/s
N = NBVSQ .^ 0.5; % rad/s
plot_model(DATAFILE, '-', 'k', []);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eventually cut data.        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_maxalt=find(abs(Z-maxalt)==min((abs(Z-maxalt))), 1, 'last'); original_zmax=Z(end);
Z=Z(1:ind_maxalt); RHO=RHO(1:ind_maxalt); TEMP=TEMP(1:ind_maxalt); SOUNDSPEED=SOUNDSPEED(1:ind_maxalt);
P=P(1:ind_maxalt); LOCALPRESSURESCALE=LOCALPRESSURESCALE(1:ind_maxalt); G=G(1:ind_maxalt); N=N(1:ind_maxalt);
KAPPA=KAPPA(1:ind_maxalt); VISCMU=VISCMU(1:ind_maxalt); MUVOL=MUVOL(1:ind_maxalt); WEAST=WEAST(1:ind_maxalt);
WNORTH=WNORTH(1:ind_maxalt); W=W(1:ind_maxalt); CP=CP(1:ind_maxalt); CV=CV(1:ind_maxalt);
GAMMA=GAMMA(1:ind_maxalt); nz=length(Z);
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

disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hydrostatic Equilibrium Treatment.                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['[',mfilename,'] Checking stability.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check Evaluation            %
% Coefficients.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tit_plus={posstr, datestr};
D=differentiation_matrix(Z, 0);

hydrostatic_ratio = @(D, P, RHO, G) (D * P) ./ (-RHO .* G);
Mach_number = @(W, C) abs(W)./C;
east_Richardson = Richardson_number(WEAST, D, N);
north_Richardson = Richardson_number(WNORTH, D, N);
proj_Richardson = Richardson_number(W, D, N);
east_Mach = Mach_number(WEAST, SOUNDSPEED);
north_Mach = Mach_number(WNORTH, SOUNDSPEED);
proj_Mach = Mach_number(W, SOUNDSPEED);
HR=hydrostatic_ratio(D, P, RHO, G);

disp(['[',mfilename,'] > Minimum Richardson number: ',num2str(min(proj_Richardson)), ' (@z=', num2str(Z(proj_Richardson==min(proj_Richardson))), ' m). ',richardson_advice]);
if(any(proj_Richardson<0.25))
  disp(['[',mfilename,', WARNING] > Projected wind''s Richardson number is < 0.25 somewhere. Unstability will probably occur.']);
else
  if(any(proj_Richardson<1))
    disp(['[',mfilename,'] > Projected wind''s Richardson number is < 1 somewhere. Unstability can occur.']);
  end
end

disp(['[',mfilename,'] > Maximum relative gap to hydrostatic state (via ratio): ', num2str(max(abs(HR-1))*100), '%.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots.                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_hydrostat_unbalance(Z, RHO, TEMP, P, G, KAPPA, VISCMU, W, tit_plus, save_plots);

% figure();
% plot(WEAST, Z, WNORTH, Z);
% % xlim([0.5 * min([east_Richardson; north_Richardson]), 2 * max([east_Richardson; north_Richardson])]);
% xlabel('wind (m/s)'); ylabel('altitude (m)');
% legend('eastward wind', 'northward wind', 'Location', 'northeast');
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

% figure();
% semilogx(east_Richardson, Z, north_Richardson, Z, proj_Richardson, Z, ones(size(Z)), Z, 'k:', 0.25*ones(size(Z)), Z, 'k');
% xlim([0.5 * min([east_Richardson; north_Richardson; proj_Richardson]), 2*max([east_Richardson; north_Richardson; proj_Richardson])]);
% xlabel('Richardson number'); ylabel('altitude (m)');
% legend('eastward Richardson number', 'northward Richardson number', 'projected wind Richardson number', 'Location', 'northeast');
% title(tit_plus);
% if save_plots == 1
%   saveas(gcf, strcat(DATAFILE,'__richardson.png'));
% end

% plot_model_effective_soundspeed(DATAFILE,0);
% if save_plots == 1
%   saveas(gcf, strcat(DATAFILE,'__effective_sound_speed.png'));
% end

% figure();
% semilogx(east_Mach, Z, north_Mach, Z, proj_Mach, Z, ones(size(Z)), Z, 'k:');
% xlim([0.5 * min([east_Mach; north_Mach; proj_Mach]), 2 * max([east_Mach; north_Mach; proj_Mach])]);
% xlabel('Mach number ($|w|/c$)'); ylabel('altitude (m)');
% legend('eastward Mach number', 'northward Mach number', 'projected wind Mach number', 'Location', 'northeast');
% title(tit_plus);
% if save_plots == 1
%   saveas(gcf, strcat(DATAFILE,'__mach.png'));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularise model.           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['[',mfilename,'] Trying regularisation.']);
modify_atmos_model; % See script.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display modifications.      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['[',mfilename,'] Displaying modifications.']);

tit_plus={strcat(['Regularised model ("', regexprep(method, '_', '\\_'), '" method)']), tit_plus{1}, tit_plus{2}};
if(strcmp(method, 'bruteforce_rho'))
  disp(['[',mfilename,'] > Used "', method, '" method.']);
  nRHO=bruteforced_RHO;
%   nRHO=smoooth(nRHO,10); disp(['[',mfilename,']   [WARNING] RHO WAS SMOOTHED.']); % Smooth using sliding Gaussian on 10 points.
%   i1=find(abs(Z-7.2e4)==min(abs(Z-7.2e4))); i2=find(abs(Z-7.32e4)==min(abs(Z-7.32e4))); work_RHO=bruteforced_RHO; work_RHO(i1:i2)=polyval(polyfit([Z(i1), Z(i2)], [work_RHO(i1), work_RHO(i2)], 1),Z(i1:i2)); nRHO=work_RHO;
%   spline=fit(Z,nRHO,'smoothingspline','smoothingparam',1e-11); nRHO=spline(Z); disp(['[',mfilename,']   [WARNING] RHO WAS SMOOTHED.']); % Smooth using splines.

  disp(['[',mfilename,'] > Maximum relative difference on RHO: ', num2str(max(abs(nRHO-RHO)./RHO)*100), '%.']);
  
  nHR=hydrostatic_ratio(D, P, nRHO, G);
  if(max(abs(nHR-1))>1e-12)
    figure(); plot(nHR, Z, ones(size(Z)), Z, 'k:');
  end
  disp(['[',mfilename,'] > New maximum relative gap to hydrostatic state (via ratio): ', num2str(max(abs(nHR-1))*100), '%.']);
  
  disp(['[',mfilename,'] > Recomputing sound speed.']);
  nSOUNDSPEED=sqrt(GAMMA.*P./bruteforced_RHO);
%   nSOUNDSPEED=smoooth(nSOUNDSPEED,20); disp(['[',mfilename,']   [WARNING] SOUNDSPEED WAS SMOOTHED (Gaussian).']); % Smooth using sliding Gaussian on 20 points.
%   spline=fit(Z,nSOUNDSPEED,'smoothingspline','smoothingparam',1e-11); nSOUNDSPEED=spline(Z); disp(['[',mfilename,']   [WARNING] SOUNDSPEED WAS SMOOTHED (spline).']); % Smooth using splines.
  spline=fit(Z,smoooth(nSOUNDSPEED,20),'smoothingspline','smoothingparam',1e-10); nSOUNDSPEED=spline(Z); disp(['[',mfilename,', WARNING] > SOUNDSPEED WAS SMOOTHED (Gaussian + spline).']); % Smooth using sliding Gaussian on 20 points, and then splines.
%   figure();
%   semilogx(SOUNDSPEED, Z, nSOUNDSPEED, Z);
%   xlim([0.5 * min([SOUNDSPEED;nSOUNDSPEED]), 2 * max([SOUNDSPEED;nSOUNDSPEED])]);
%   xlabel('$c$'); ylabel('altitude (m)');
%   legend('old $c$', 'new $c$', 'Location', 'northeast');
%   title(tit_plus);
  
  disp(['[',mfilename,'] > Recomputing Brunt-Väisälä frequency.']);
  nN=sqrt((GAMMA-1).*(G./nSOUNDSPEED).^2); % rad/s
%   figure();
%   semilogx(N.^2, Z, nN.^2, Z);
%   xlim([0.5 * min([N.^2;nN.^2]), 2 * max([N.^2;nN.^2])]);
%   xlabel('$N^2$ (rad/s)'); ylabel('altitude (m)');
%   legend('old $N^2$', 'new $N^2$', 'Location', 'best');
%   title(tit_plus);
  
  nTEMP=TEMP;nP=P;nLOCALPRESSURESCALE=LOCALPRESSURESCALE;nG=G;nKAPPA=KAPPA;nVISCMU=VISCMU;nMUVOL=MUVOL;nWNORTH=WNORTH;nWEAST=WEAST;nW=W;nCP=CP;nCV=CV;nGAMMA=GAMMA; % Not modified.

  plot_hydrostat_unbalance(Z, nRHO, nTEMP, nP, nG, nKAPPA, nVISCMU, nW, tit_plus, save_plots);
  
  Richard=Richardson_number(W, D, N);
  nRichard=Richardson_number(W, D, nN);
  disp(['[',mfilename,'] > New minimum Richardson number: ',num2str(min(nRichard)), ' (@z=', num2str(Z(nRichard==min(nRichard))), ' m). ', richardson_advice]);
  figure();
  semilogx(Richard, Z, nRichard, Z, ones(size(Z)), Z, 'k:', 0.25*ones(size(Z)), Z, 'k');
  xlim([0.5 * min([Richard;nRichard]), 2 * max([Richard;nRichard])]);
  xlabel('Richardson number'); ylabel('altitude (m)');
  legend('old Richardson number', 'new Richardson number', 'Location', 'northeast');
  title(tit_plus);
  if save_plots == 1
    saveas(gcf, strcat(DATAFILE,'__new_richardson.png'));
  end
  
%   Mach=Mach_number(W, SOUNDSPEED);
%   nMach=Mach_number(W, nSOUNDSPEED);
%   figure();
%   semilogx(Mach, Z, nMach, Z, ones(size(Z)), Z, 'k:');
%   xlim([0.5 * min([Mach;nMach]), 2 * max([Mach;nMach])]);
%   xlabel('Mach number'); ylabel('altitude (m)');
%   legend('old Mach number', 'new Mach number', 'Location', 'northeast');
%   title(tit_plus);
%   if save_plots == 1
%     saveas(gcf, strcat(DATAFILE,'__new_mach.png'));
%   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output regularised model to %
% file.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['[',mfilename,'] Outputting regularised model to file.']);

if(maxalt<original_zmax)
  disp(['[',mfilename,', WARNING] maxalt = ', num2str(maxalt),' < original max_alt = ', num2str(original_zmax),', the regularised model might not go as far up as the original model. Be careful, or re-run script with maxalt=Inf.']);
end
if(interpolate~=0)
  disp(['[',mfilename,', WARNING] Interpolation occured, the regularised model might not have the same resolution as the original model. Be careful, or re-run script with interpolate=0.']);
end
decision=-1;
while(not(ismember(decision,[0,1])))
  decision=input(['[',mfilename,'] > Output regularised model to another file? (0 for no, 1 for yes) > ']);
end
if(decision==0)
  disp(['[',mfilename,'] > Outputting cancelled, stopping script.']); return;
end

decision_w0atz0=-1;
while(not(ismember(decision_w0atz0,[0,1])))
  decision_w0atz0=input(['[',mfilename,'] > Force w(z=0)=0? (0 for no, 1 for yes) > ']);
end
if(decision_w0atz0)
%   width=20*max(diff(Z)); % Spread over 20 steps.
  width=input(['[',mfilename,'] > Height of apodisation (in meters)? > ']);
%   b=-1.82138636771844967304021031862099524348122888360095; % 1e-2 level.
  b=-2.75106390571206079614551316854267817109677902559646; % 1e-3 level.
%   b=-3.45891073727950002215092763595756951991566980804288674707621013; % 1e-6 level.
  a=-2*b/width;
  apoWind=erf(a*Z+b)/2+0.5;
  apoWind(Z==0)=0;
  nW=apoWind.*nW;
end

decision_mu=-1;
while(not(ismember(decision_mu,[0,1])))
  decision_mu=input(['[',mfilename,'] > Force mu=0? (0 for no, 1 for yes) > ']);
end
decision_kap=-1;
while(not(ismember(decision_kap,[0,1])))
  decision_kap=input(['[',mfilename,'] > Force kappa=0? (0 for no, 1 for yes) > ']);
end

if(decision_mu && not(decision_kap))
  nVISCMU=0*nVISCMU;prefix='reg+mu0_'; disp(['[',mfilename,']   [WARNING] MU WAS SET TO ZERO.']);
elseif(decision_kap && not(decision_mu))
  nKAPPA=0*nKAPPA;prefix='reg+kappa0_'; disp(['[',mfilename,']   [WARNING] KAPPA WAS SET TO ZERO.']);
elseif(decision_mu && decision_kap)
  nVISCMU=0*nVISCMU;nKAPPA=0*nKAPPA;prefix='reg+mukappa0_'; disp(['[',mfilename,']   [WARNING] MU AND KAPPA WERE SET TO ZERO.']);
else
  prefix='reg_';
end
SPL=split(DATAFILE, '/'); SPL(end)=strcat(prefix, SPL(end)); nDATAFILE=join(SPL, '/');nDATAFILE=nDATAFILE{1};
rewrite_atmos_model(nDATAFILE, DATAFILE, Z, nRHO, nTEMP, nSOUNDSPEED, nP, nLOCALPRESSURESCALE, nG, nN.^2, nKAPPA, nVISCMU, nMUVOL, nWNORTH, nWEAST, nW, nCP, nCV, nGAMMA);
plot_model(nDATAFILE, '-', 'k', []);
if(decision_w0atz0)
  nnRichard=Richardson_number(nW, D, nN);
  disp(['[',mfilename,'] > New minimum Richardson number after setting w(z=0)=0: ',num2str(min(nnRichard)), ' (@z=', num2str(Z(nnRichard==min(nnRichard))), ' m). ', richardson_advice]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sQ=smoooth(Q,n)
  sQ=smoothdata(Q,'gaussian',n);
%   disp(strcat("[WARNING] Smoothed quantity ", inputname(1), " by a Gaussian filter over ", num2str(n), " elements."));
end

function [term1, term2, term3]=hydrostat_unbalance(Z, RHO, TEMP, P, G, KAPPA, VISCMU, W)
  D=differentiation_matrix(Z, 0);
  term1 = abs(D * (VISCMU .* (D * W)));
  term2 = abs(D * (KAPPA .* (D * TEMP) + W .* VISCMU .* (D * W)));
  term3 = abs(D * P + RHO .* G);
end

function plot_hydrostat_unbalance(Z, RHO, TEMP, P, G, KAPPA, VISCMU, W, tit_plus, save_plots)
  figure();
  [term1, term2, term3]=hydrostat_unbalance(Z, RHO, TEMP, P, G, KAPPA, VISCMU, W);
  semilogx(term1, Z, 'r'); hold on;
  semilogx(term2, Z, 'g');
  semilogx(term3, Z, 'b');
  tmp_valmat=cell2mat(get(get(gca, 'children'), 'XData')); tmp_plot_maxval=max(max(tmp_valmat)); tmp_valmat(tmp_valmat==0)=Inf; tmp_plot_minval=min(min(tmp_valmat));
  xlim([0.5*tmp_plot_minval, 2*tmp_plot_maxval]);
  xlabel({'amplitude of hydrostatic unbalance terms', '(projected wind)'}); ylabel('altitude (m)');
  legend('$\left|\partial_z\left(\mu\partial_zw\right)\right|$', '$\left|\partial_z\left(\kappa\partial_zT+w\mu\partial_zw\right)\right|$', '$\left|\partial_zp + \rho g_z\right|$', 'Location', 'best');
  title(tit_plus);
  if save_plots == 1
    saveas(gcf, strcat(DATAFILE,'__unbalance_terms.png'));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Older DATAFILE paths.       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/atmospheric/stratospheric/0_200000_301_66.56306_0.00000_0_173_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/atmospheric/stratospheric/0_200000_301_66.56306_0.00000_0_356_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/atmospheric/stratospheric/0_200000_301_45.00000_0.00000_0_173_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/atmospheric/stratospheric/reg_0_200000_301_45.00000_0.00000_0_173_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/atmospheric/stratospheric/0_200000_301_45.00000_0.00000_0_356_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/atmospheric/stratospheric/0_200000_301_23.43694_0.00000_0_173_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/atmospheric/stratospheric/0_200000_301_23.43694_0.00000_0_356_43200.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/atmospheric/stratospheric/0_200000_301_0.00000_0.00000_0_173_0.00000_0.00000"; headerlines=3;
% DATAFILE = "/home/l.martire/Documents/SPECFEM/Ongoing_Work/atmospheric/stratospheric/0_200000_301_0.00000_0.00000_0_173_43200.00000_0.00000"; headerlines=3;
