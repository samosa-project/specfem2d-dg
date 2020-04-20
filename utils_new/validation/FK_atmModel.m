if(externalDGAtmosModel)
  % External DG atmospheric model, load it.
  disp(['[',mfilename,'] Loading atmospheric model from external file in ''',OFd,'''.']);
  [ALT, RHO, ~, SOUNDSPEED, ~, ~, GRA, NSQ, ~, ~, ~, ~, ~, WIND, ~, ~, ~] = extract_atmos_model(atmmodelfile, 3, 0, 0); % plot_model(atmmodelfile,'-','k',[]);
%   H        = (ALT(2)-ALT(1))/log(RHO(1)/RHO(2))
  % Set H.
  H        = (-ALT(2))/log(RHO(2)/RHO(1)); % NEEDED FOR NSQ AND ANALYTIC SOLUTION
  disp(['[',mfilename,'] Scale height computed,         = ',num2str(H),'.']);
%   rho_0    = RHO(1); % used nowhere
%   Nsqt      = -(GRA(1)*( -1/H - GRA(1)/(velocity(1)^2) )); % used nowhere
  % Set NSQ.
  error(['[',mfilename,', ERROR] N^2 computation for external models not implemented yet.']);
  % Set wind_x.
  if(max(abs(diff(WIND)))==0)
    wind_x = WIND(1);
    disp(['[',mfilename,'] Wind_x loaded,                 = ',num2str(wind_x),'. Constant with altitude.']);
  else
    wind_x = WIND(1);
    disp(['[',mfilename,'] Wind_x loaded,                 = ',num2str(wind_x),'. Not constant with altitude, setting wind_x as wind at altitude 0.']);
  end
else
  % Internal model, load it.
  disp(['[',mfilename,'] Building atmospheric model from parfile in ''',OFd,'''.']);
  % Set H.
  H = readExampleFiles_extractParam(parfile, 'SCALE_HEIGHT', 'float'); % NEEDED FOR NSQ AND ANALYTIC SOLUTION
  disp(['[',mfilename,'] Scale height loaded,           = ',num2str(H),'.']);
  USE_ISOTHERMAL_MODEL = readExampleFiles_extractParam(parfile, 'USE_ISOTHERMAL_MODEL', 'bool');
  if(USE_ISOTHERMAL_MODEL)
    % isothermal case
    disp(['[',mfilename,'] Isothermal case.']);
%     error(['[',mfilename,', ERROR] Isothermal case not implemented yet.']);
    GAM = readExampleFiles_extractParam(parfile, 'constant_p', 'float') / readExampleFiles_extractParam(parfile, 'constant_v', 'float');
    GRA = readExampleFiles_extractParam(parfile, 'gravity', 'float'); % NEEDED FOR NSQ
    SOUNDSPEED = sqrt(GAM .* GRA .* H); % c^2=\gamma*p/\rho, but p=rho*g*H (hydrostaticity) thus c^2=\gamma*rho*g*h/rho=\gamma*g*h
  else
    % isobaric case
    disp(['[',mfilename,'] Isobaric case.']);
    SOUNDSPEED = readExampleFiles_extractParam(parfile, 'sound_velocity', 'float');
    disp(['[',mfilename,'] Speed of sound loaded,         = ',num2str(SOUNDSPEED),'.']);
    GAM = []; % not needed, but keep for completeness
    GRA = 0; % NEEDED FOR NSQ
    disp(['[',mfilename,'] Gravity set for isobaric case, = ',num2str(GRA),'.']);
  end
  % Set NSQ.
  NSQ = -(GRA*( -1/H - GRA/(SOUNDSPEED^2) )); % NOT SURE ABOUT THAT FORMULA, CHECK
  disp(['[',mfilename,'] N^2 computed,                  = ',num2str(NSQ),'.']);
  % Set wind_x.
  wind_x = readExampleFiles_extractParam(parfile, 'wind', 'float');
  disp(['[',mfilename,'] Wind_x loaded,                 = ',num2str(wind_x),'.']);
end
disp(['[',mfilename,'] Setting wind_y = wind_z = 0.']);
wind_y = 0;
wind_z = 0;