% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [ALTITUDE, DENSITY, TEMPERATURE, SOUNDSPEED, PRESSURE, LOCALPRESSURESCALE, G, NBVSQ, KAPPA, MU, MUVOL, NORTHWIND, EASTWIND, WIND, CP, CV, GAMMA, FR, SVIB] = extract_atmos_model(DATAFILE, headerlines, interpolate, interp_delta)
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
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Interpolate.                %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(interpolate==1)
    % Setup.
    % ALTITUDE = ALTITUDE - ALTITUDE(1); % Set first altitude to 0.
    if min(diff(ALTITUDE)) < interp_delta
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
    interp_LOCALPRESSURESCALE = interp1(ALTITUDE, LOCALPRESSURESCALE, interp_ALTITUDE)';    % (Htab)
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
    interp_WIND = interp1(ALTITUDE, WIND, interp_ALTITUDE)';                          % projected wind
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
    WIND = interp_WIND;
    CP = interp_CP;
    CV = interp_CV;
    GAMMA = interp_GAMMA;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end