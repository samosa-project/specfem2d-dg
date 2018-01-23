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
interp_delta = 1000; % Interpolation step (m).
save_plots = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set datafile.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp(['[INPUT] Path (relative or absolute) to datafile to be loaded? You''re in ''', pwd ,'''.']); DATAFILE = input(' > ');
% DATAFILE = 'Mars/Mars_atmosphere_models/MARS_MCD_AGW_model_INSIGHT_1_Oct_2016_12h.dat'; headerlines=1;
% DATAFILE = 'Mars/Mars_atmos_models_LEO/msise_Sumatra_model_DG.dat'; headerlines=1; % WARNING: DATA FORMAT IS NOT THE SAME AS THE ONE IN MARS MODELS!
DATAFILE = "Earth/wrapper/msisehwm_model_output"; headerlines=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and store data.        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(DATAFILE);
if(fid==-1)
  error(strcat("Cannot open file ", DATAFILE,').'))
end
stop=0;
while(stop==0)
  line = fgetl(fid);
  year=str2num(regexprep(regexprep(line, 'day.*', ''), 'year', ''));
  daysincenewyear=str2num(regexprep(regexprep(line, 'seconds.*', ''), 'year *[0-9]+ day', ''));
  secondssincenewday=str2num(regexprep(regexprep(line, 'lat.*', ''), 'year *[0-9]+ day *[0-9]+ seconds', ''));
  datestr=[num2str(daysincenewyear), 'th day of ', num2str(year), ' at ' num2str(floor(secondssincenewday/3600)),':',num2str(floor((secondssincenewday - floor(secondssincenewday/3600)*3600)/60)), ' UT'];
  
  lat=str2num(regexprep(regexprep(line, 'year *[0-9]+ day *[0-9]+ seconds *[0-9]+\.[0-9]+ *lat', ''), 'lon.*', ''));
  lon=str2num(regexprep(line, 'year *[0-9]+ day *[0-9]+ seconds *[0-9]+\.[0-9]+ *lat *[0-9]+\.[0-9]+ *lon', ''));
  posstr=['lat. ',num2str(lat), ', lon. ',num2str(lon)];
  stop=1;
end
fclose('all');

DATA = importdata(DATAFILE, ' ', headerlines);
ALTITUDE = DATA.data(:, 1); %                                     (z)
DENSITY = DATA.data(:, 2); %                                      (rhoat)
TEMPERATURE = DATA.data(:, 3); %                                  (T)
SOUNDSPEED = DATA.data(:, 4); %                                   (c)
PRESSURE = DATA.data(:, 5); %                                     (p)
LOCALPRESSURESCALE = DATA.data(:, 6); % Unused.                   (H)
G = DATA.data(:, 7); %                                            (gravity)
NBVSQ = DATA.data(:, 8); %                                        (Nsqtab)
KAPPA = DATA.data(:, 9); %                                        (kappa)
MU = DATA.data(:, 10); % dynamic viscosity                        (mu)
MUVOL = DATA.data(:, 11); % volumic viscosity?                    (mu_vol)
NORTHWIND = DATA.data(:, 12); %                                   (w_M)
EASTWIND = DATA.data(:, 13); %                                    (w_Z)
CP = DATA.data(:, 14); %                                          (c_p)
CV = DATA.data(:, 15); %                                          (c_v)
GAMMA = DATA.data(:, 16); %                                       (gamma)
% NCUTA = DATA.data(:, 9); %                                        (Ncuttab)
% MUVOLROT = DATA.data(:, 13); % rotational volumic viscosity?      (mu_volrottab)
% NORTHWIND = DATA.data(:, 16); %                                   (Wind(1,iz))
% EASTWIND=NORTHWIND;
% EASTWIND = DATA.data(:, 17); %                                    (Wind(2,iz))
% CP = DATA.data(:, 18); % Unused.                                (Cp)
% CV = DATA.data(:, 19); % Unused.                                (Cv)
% GAMMA = DATA.data(:, 20); % Unused.                             (gammatab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(interpolate==1)
  % Setup.
  % ALTITUDE = ALTITUDE - ALTITUDE(1); % Set first altitude to 0.
  if (ALTITUDE(2) - ALTITUDE(1)) < interp_delta
    disp('[WARNING] Interpolation step is greater than data step. Undersampling will occur.');
  end
  % interp_ALTITUDE=[ALTITUDE(1):interp_delta:ALTITUDE(end)]; % Using this, it is possible that interp_ALTITUDE(end)<ALTITUDE(end).
  interp_ALTITUDE=[ALTITUDE(1):interp_delta:ALTITUDE(1)+ceil((ALTITUDE(end)-ALTITUDE(1)) / interp_delta)*interp_delta]; % Using this, interp_ALTITUDE(end)>=ALTITUDE(end), always.

  % Logarithmically interpolated quantities.
  interp_DENSITY = exp(interp1(ALTITUDE, log(DENSITY), interp_ALTITUDE))';                  % (rhoat)
  interp_PRESSURE = exp(interp1(ALTITUDE, log(PRESSURE), interp_ALTITUDE))';                % (P)
  % interp_VIBRATTENUATION = exp(interp1(ALTITUDE, log(VIBRATTENUATION), interp_ALTITUDE))';  % (fr)

  % Linearly interpolated quantities.
  interp_TEMPERATURE = interp1(ALTITUDE, TEMPERATURE, interp_ALTITUDE)';                    % (T)
  interp_SOUNDSPEED = interp1(ALTITUDE, SOUNDSPEED, interp_ALTITUDE)';                      % (v)
  % interp_LOCALPRESSURESCALE = interp1(ALTITUDE, LOCALPRESSURESCALE, interp_ALTITUDE)';    % (Htab)
  interp_G = interp1(ALTITUDE, G, interp_ALTITUDE)';                                        % (gravity)
  interp_NBVSQ = interp1(ALTITUDE, NBVSQ, interp_ALTITUDE)';                                % (Nsqtab)
  % interp_NCUTA = interp1(ALTITUDE, NCUTA, interp_ALTITUDE)';                                % (Ncuttab)
  interp_KAPPA = interp1(ALTITUDE, KAPPA, interp_ALTITUDE)';                                % (Kappatab)
  interp_MU = interp1(ALTITUDE, MU, interp_ALTITUDE)';                                      % (MUtab)
  interp_MUVOL = interp1(ALTITUDE, MUVOL, interp_ALTITUDE)';                                % (MUvoltab)
  % interp_MUVOLROT = interp1(ALTITUDE, MUVOLROT, interp_ALTITUDE)';                          % (MUvolrottab)
  % interp_SVIB = interp1(ALTITUDE, SVIB, interp_ALTITUDE)';                                  % (Svib)
  interp_NORTHWIND = interp1(ALTITUDE, NORTHWIND, interp_ALTITUDE)';                        % (Wind(1,iz))
  interp_EASTWIND = interp1(ALTITUDE, EASTWIND, interp_ALTITUDE)';                          % (Wind(2,iz))
  interp_CP = interp1(ALTITUDE, CP, interp_ALTITUDE)';                                    % (Cp)
  interp_CV = interp1(ALTITUDE, CV, interp_ALTITUDE)';                                    % (Cv)
  interp_GAMMA = interp1(ALTITUDE, GAMMA, interp_ALTITUDE)';                              % (gammatab)
  
  ALTITUDE = interp_ALTITUDE;
  DENSITY = interp_DENSITY;
  TEMPERATURE = interp_TEMPERATURE;
  SOUNDSPEED = interp_SOUNDSPEED;
  PRESSURE = interp_PRESSURE;
  LOCALPRESSURESCALE = interp_LOCALPRESSURESCALE;
  G = interp_G;
  NBVSQ = interp_NBVSQ;
  KAPPA = interp_KAPPA;
  MU = interp_MU;
  MUVOL = interp_MUVOL;
  NORTHWIND = interp_NORTHWIND;
  EASTWIND = interp_EASTWIND;
  CP = interp_CP;
  CV = interp_CV;
  GAMMA = interp_GAMMA;
end
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
% Data.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxalt=Inf; % If one has to cut data, choose maximum altitude here. Put Inf if all altitudes are to be considered.
Z = ALTITUDE;
RHO = DENSITY;
TEMP = TEMPERATURE;
P = PRESSURE;
G = G;
N = NBVSQ .^ 0.5 / (2 * pi);
KAPPA = KAPPA;
VISCMU = MU;
WEAST = EASTWIND;
WNORTH = NORTHWIND;
% Z = ALTITUDE;
% RHO = DENSITY;
% TEMP = TEMPERATURE;
% P = PRESSURE;
% G = G;
% N = NBVSQ .^ 0.5 / (2 * pi);
% KAPPA = KAPPA;
% VISCMU = MU;
% WEAST = EASTWIND;
% WNORTH = NORTHWIND;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differentiation Matrix.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nz=length(Z);
D=differentiation_matrix(Z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check Evaluation            %
% Coefficients.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tit_plus={posstr, datestr};

figure();
semilogx(abs(D * (VISCMU .* (D * WEAST))), Z, 'r'); hold on;
semilogx(abs(D * (KAPPA .* (D * TEMP) + WEAST .* VISCMU .* (D * WEAST))), Z, 'g');
semilogx(abs(D * P + RHO .* G), Z, 'b');
tmp_valmat=cell2mat(get(get(gca, 'children'), 'XData')); tmp_plot_maxval=max(max(tmp_valmat)); tmp_valmat(tmp_valmat==0)=Inf; tmp_plot_minval=min(min(tmp_valmat));
xlim([0.5*tmp_plot_minval, 2*tmp_plot_maxval]);
xlabel({'amplitude of hydrostatic unbalance terms', '(projection on zonal winds)'}); ylabel('altitude (m)');
legend('$\left|\partial_z\left(\mu\partial_zw\right)\right|$', '$\left|\partial_z\left(\kappa\partial_zT+w\mu\partial_zw\right)\right|$', '$\left|\partial_zp + \rho g_z\right|$', 'Location', 'best');
title(tit_plus);
if save_plots == 1
  saveas(gcf, strcat(DATAFILE,'__unbalance_terms.png'));
end

figure();
plot(WEAST, Z, WNORTH, Z);
% xlim([0.5 * min([east_Richardson; north_Richardson]), 2 * max([east_Richardson; north_Richardson])]);
xlabel('wind (m/s)'); ylabel('altitude (m)');
legend('eastward wind', 'northward wind', 'Location', 'NorthWest');
title(tit_plus);
if save_plots == 1
  saveas(gcf, strcat(DATAFILE,'__winds.png'));
end

hydrostatic_ratio = @(D, P, RHO, G) (D * P) ./ (-RHO .* G);
east_Richardson = Richardson_number(WEAST, D, N);
north_Richardson = Richardson_number(WNORTH, D, N);

figure();
HR=hydrostatic_ratio(D, P, RHO, G);
plot(HR(2:end), Z(2:end), ones(size(Z)), Z, 'k:');
xlabel('$\partial_zp/(-\rho g)$'); ylabel('altitude (m)');
% legend('hydrostatic ratio - 1', 'Location', 'best');
title(tit_plus);

figure();
semilogx(east_Richardson, Z, north_Richardson, Z, ones(size(Z)), Z, 'k:');
xlim([0.5 * min([east_Richardson; north_Richardson]), 2 * max([east_Richardson; north_Richardson])]);
xlabel('Richardson number'); ylabel('altitude (m)');
legend('eastward Richardson number', 'northward Richardson number', 'Location', 'NorthWest');
title(tit_plus);
if save_plots == 1
  saveas(gcf, strcat(DATAFILE,'__richardson.png'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularise model.           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method='integrate';
% method='bruteforce';
modify_atmospheric_model % See script.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%