% Author:        Léo Martire.
% Description:   Loads an atmospheric model from a atmospheric_model.dat
%                formatted (SPECFEM-DG) file.
% Notes:         TODO.
%
% Usage:
%   [Z, RHO, TEMP, C, P, H, G, NBVSQ, KAPPA, MU, MUVOL, NORTHWIND, ...
%    EASTWIND, WIND, CP, CV, GAMMA, FR, SVIB] = ...
%   extract_atmos_model(DATAFILE, headerlines, interpolate, interp_delta)
% with:
%   TODO.
% yields:
%   TODO.

function [Z, RHO, TEMP, C, P, H, G, NBVSQ, KAPPA, MU, MUVOL, NORTHWIND, EASTWIND, WIND, CP, CV, GAMMA, FR, SVIB] = extract_atmos_model(DATAFILE, headerlines, interpolate, interp_delta)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Load and store data.        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  DATA = importdata(DATAFILE, ' ', headerlines);
  frsvibfound=0;
  if(size(DATA.data,2)<17)
    error('not enough columns (must be either 17 or 19)');
  elseif(size(DATA.data,2)==17)
    % ok
    frsvibfound=0;
  elseif(size(DATA.data,2)>17)
    if(size(DATA.data,2)==19)
      % ok
      frsvibfound=1;
    else
      error('wrong number of columns (must be either 17 or 19)');
    end
  end
  Z = DATA.data(:, 1); %                                     (z)
  RHO = DATA.data(:, 2); %                                      (rhoat)
  TEMP = DATA.data(:, 3); %                                  (T)
  C = DATA.data(:, 4); %                                   (c)
  P = DATA.data(:, 5); %                                     (p)
  H = DATA.data(:, 6); % Unused.                   (H)
  G = DATA.data(:, 7); %                                            (gravity)
  NBVSQ = DATA.data(:, 8); %                                        (Nsqtab)
  KAPPA = DATA.data(:, 9); %                                        (kappa)
  MU = DATA.data(:, 10); % dynamic viscosity                        (mu)
  MUVOL = DATA.data(:, 11); % volumic viscosity?                    (mu_vol)
  NORTHWIND = DATA.data(:, 12); %                                   (w_M) meridional = towards north
  EASTWIND = DATA.data(:, 13); %                                    (w_Z) zonal = towards east
  WIND = DATA.data(:, 14); %                                        (w_P)
  CP = DATA.data(:, 15); %                                          (c_p)
  CV = DATA.data(:, 16); %                                          (c_v)
  GAMMA = DATA.data(:, 17); %                                       (gamma)
  % NCUTA = DATA.data(:, 9); %                                        (Ncuttab)
  % MUVOLROT = DATA.data(:, 13); % rotational volumic viscosity?      (mu_volrottab)
  % NORTHWIND = DATA.data(:, 16); %                                   (Wind(1,iz))
  % EASTWIND=NORTHWIND;
  % EASTWIND = DATA.data(:, 17); %                                    (Wind(2,iz))
  % CP = DATA.data(:, 18); % Unused.                                (Cp)
  % CV = DATA.data(:, 19); % Unused.                                (Cv)
  % GAMMA = DATA.data(:, 20); % Unused.                             (gammatab)
  if(frsvibfound)
    FR = DATA.data(:, 18); %                                        (f_r)
    SVIB = DATA.data(:, 19); %                                      (s_vib)
  else
    FR = [];
    SVIB = [];
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Interpolate.                %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(interpolate==1)
    % Setup.
    % ALTITUDE = ALTITUDE - ALTITUDE(1); % Set first altitude to 0.
    if min(diff(Z)) < interp_delta
      disp('[WARNING] Interpolation step is greater than data step. Undersampling will occur.');
    end
    % interp_ALTITUDE=[ALTITUDE(1):interp_delta:ALTITUDE(end)]; % Using this, it is possible that interp_ALTITUDE(end)<ALTITUDE(end).
    interp_ALTITUDE=[Z(1):interp_delta:Z(1)+ceil((Z(end)-Z(1)) / interp_delta)*interp_delta]; % Using this, interp_ALTITUDE(end)>=ALTITUDE(end), always.

    % Logarithmically interpolated quantities.
    interp_DENSITY = exp(interp1(Z, log(RHO), interp_ALTITUDE))';                  % (rhoat)
    interp_PRESSURE = exp(interp1(Z, log(P), interp_ALTITUDE))';                % (P)
    % interp_VIBRATTENUATION = exp(interp1(ALTITUDE, log(VIBRATTENUATION), interp_ALTITUDE))';  % (fr)

    % Linearly interpolated quantities.
    interp_TEMPERATURE = interp1(Z, TEMP, interp_ALTITUDE)';                    % (T)
    interp_SOUNDSPEED = interp1(Z, C, interp_ALTITUDE)';                      % (v)
    interp_LOCALPRESSURESCALE = interp1(Z, H, interp_ALTITUDE)';    % (Htab)
    interp_G = interp1(Z, G, interp_ALTITUDE)';                                        % (gravity)
    interp_NBVSQ = interp1(Z, NBVSQ, interp_ALTITUDE)';                                % (Nsqtab)
    % interp_NCUTA = interp1(ALTITUDE, NCUTA, interp_ALTITUDE)';                                % (Ncuttab)
    interp_KAPPA = interp1(Z, KAPPA, interp_ALTITUDE)';                                % (Kappatab)
    interp_MU = interp1(Z, MU, interp_ALTITUDE)';                                      % (MUtab)
    interp_MUVOL = interp1(Z, MUVOL, interp_ALTITUDE)';                                % (MUvoltab)
    % interp_MUVOLROT = interp1(ALTITUDE, MUVOLROT, interp_ALTITUDE)';                          % (MUvolrottab)
    % interp_SVIB = interp1(ALTITUDE, SVIB, interp_ALTITUDE)';                                  % (Svib)
    interp_NORTHWIND = interp1(Z, NORTHWIND, interp_ALTITUDE)';                        % (Wind(1,iz))
    interp_EASTWIND = interp1(Z, EASTWIND, interp_ALTITUDE)';                          % (Wind(2,iz))
    interp_WIND = interp1(Z, WIND, interp_ALTITUDE)';                          % projected wind
    interp_CP = interp1(Z, CP, interp_ALTITUDE)';                                    % (Cp)
    interp_CV = interp1(Z, CV, interp_ALTITUDE)';                                    % (Cv)
    interp_GAMMA = interp1(Z, GAMMA, interp_ALTITUDE)';                              % (gammatab)

    Z = interp_ALTITUDE;
    RHO = interp_DENSITY;
    TEMP = interp_TEMPERATURE;
    C = interp_SOUNDSPEED;
    P = interp_PRESSURE;
    H = interp_LOCALPRESSURESCALE;
    G = interp_G;
    NBVSQ = interp_NBVSQ;
    KAPPA = interp_KAPPA;
    MU = interp_MU;
    MUVOL = interp_MUVOL;
    NORTHWIND = interp_NORTHWIND;
    EASTWIND = interp_EASTWIND;
    WIND = interp_WIND;
    CP = interp_CP;
    CV = interp_CV;
    GAMMA = interp_GAMMA;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end