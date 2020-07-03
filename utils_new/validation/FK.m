% Author:        Quentin Brissaud, Léo Martire.
% Description:   TODO.
% Notes:         [Fritts and Alexander, 2003, Brissaud et al., 2016].
%
% Usage:
%   Use this file to plot on seismograms traces exiting from specfem2D in 
%   ASCII filesand a synthetic signal waited at the station in case of 
%   simpler homogeneous medium with a variable density and a constant 
%   velocity.
%   You will choose the directory in which the ASCII files are.
%   Then choose the number of receivers with the velocity of the signal and 
%   finally run this m-file.
% with:
%   TODO.
% yields:
%   TODO.

% [Brissaud et al., 2016] Brissaud, Q., Martin, R., Garcia, R. F., & Komatitsch, D. (2016). Finite-difference numerical modelling of gravitoacoustic wave propagation in a windy and attenuating atmosphere. Geophysical Journal International, 206(1), 308–327. https://doi.org/10.1093/gji/ggw121

clear all;
% close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SPCFMloc = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/';
addpath(genpath([SPCFMloc,'utils_new']));
% factor_err = 100; % factor by which multiply difference in plots.
prefix = 'validation_lns_fk';
% fignames
savefigname_base = ['compare2FKanalytic'];
savefigpath = ['/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/',prefix,'_info/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rootd=[SPCFMloc,'EXAMPLES/',prefix,'_1.00dx_1.00dt/']; % EXAMPLE path
% rootd=[SPCFMloc,'EXAMPLES/',prefix,'_1.00dx_0.50dt/']; % EXAMPLE path
% rootd=[SPCFMloc,'EXAMPLES/',prefix,'_0.50dx_0.50dt/']; % EXAMPLE path
% OFd = [rootd,'OUTPUT_FILES/'];
% OFd = [rootd,'OUTPUT_FILES_isobaric_LNS_st1']; % ISOBARIC CASE
% OFd = [rootd,'OUTPUT_FILES_isothermal_LNS_st1'];

% rootd=[SPCFMloc,'EXAMPLES/',prefix,'_isobaric/']; % EXAMPLE path
rootd=[SPCFMloc,'EXAMPLES/',prefix,'_isothermal/']; % EXAMPLE path
% rootd=[SPCFMloc,'EXAMPLES/',prefix,'_isothermal_shorter/']; % EXAMPLE path
OFd = [rootd, 'OUTPUT_FILES'];
% OFd = [rootd, 'OUTPUT_FILES_LNS'];
% OFd = [rootd, 'OUTPUT_FILES_FNS'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beging program.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subsample = 0; % Subsample synthetics? See sumsamble_dt below, near T_0.
% Pretreatment.
simuType = 1; % 1 for forcing and 2 for attenuation
if(not(OFd(end)=='/')); OFd = [OFd, '/']; end;
% atmmodelfile=[directory,'1D_iso_enhanced.txt'];
atmmodelfile = [OFd, 'input_atmospheric_model.dat'];
parfile = [OFd, 'input_parfile'];
int_file = [OFd, 'input_interfaces'];

% Automatically load parameters from simulation's OUTPUT_FILES directory.
disp([' ']); disp(['[',mfilename,', INFO] Loading simulation parameters.']);
% [zmin, zmax] = readExampleFiles_zMinMaxInterfacesFile(int_file);
[~, ~, zmin, zmax] = readExampleFiles_meshfem_mesh(int_file);
dt = readExampleFiles_extractParam(parfile, 'DT', 'float');
TYPE_FORCING = readExampleFiles_extractParam(parfile, 'TYPE_FORCING', 'int');
P_0 = readExampleFiles_extractParam(parfile, 'main_spatial_period', 'float'); % main spatial period
T_0 = readExampleFiles_extractParam(parfile, 'main_time_period', 'float'); % main time period
X_0 = readExampleFiles_extractParam(parfile, 'forcing_initial_loc', 'float'); % forcing initial location
t_0 = readExampleFiles_extractParam(parfile, 'forcing_initial_time', 'float'); % forcing initial time
FORCING_FACTOR = readExampleFiles_extractParam(parfile, 'FORCING_DG_FACTOR', 'float');
seismotype = readExampleFiles_extractParam(parfile, 'seismotype', 'int');
xmin = readExampleFiles_extractParam(parfile, 'xmin', 'float');
xmax = readExampleFiles_extractParam(parfile, 'xmax', 'float');
ymin = 0; ymax = 0; % Setting this for 2D simulations. For 3D simulations, another call to readExampleFiles_extractParam would be needed.
nx = readExampleFiles_extractParam(parfile, 'nx', 'int');
subsample_dt = T_0/100; % if subsample == 1, wanted dt for subsampling: ~30 pts/period or higher is safe
% Atmospheric model.
switch(readExampleFiles_extractParam(parfile, 'MODEL', 'string'))
  case 'default'
    externalDGAtmosModel = 0;
  case 'external_DG'
    externalDGAtmosModel = 1;
  otherwise
    error('MODEL not implemented.');
end
% Deduce other useful parameters.
dx = (xmax-xmin)/nx;
dy = dx; disp(['[',mfilename,'] Setting dy=dx. Make sure this is the case for your simulation.']);
t0 = 0.0;
% Safeguards.
if(any([readExampleFiles_extractParam(parfile, 'dynamic_viscosity', 'float'), readExampleFiles_extractParam(parfile, 'thermal_conductivity', 'float')]~=0))
  error(['[',mfilename,', ERROR] Simulation was viscous (dynamic_viscosity and/or thermal_conductivity != 0), this case is not implemented yet.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Atmospheric Parameters.
disp([' ']); disp(['[',mfilename,', INFO] Preparing atmospheric model.']);
FK_atmModel;

% Loading Results.
disp([' ']); disp(['[',mfilename,', INFO] Loading synthetics.']);
% data = readAndSubsampleSynth(OFd, 2, 'BXZ', 'semv', 1, 1, 0); sig_t = data(:,1)'; sig_v = data(:,2)'; figure(); plot(sig_t, sig_v); % test
[~, signal_type, vname_long, vunit] = seismotypeNames(seismotype, readExampleFiles_extractParam(parfile, 'USE_DISCONTINUOUS_METHOD', 'bool'));
[x_stat, z_stat, y_stat, ~] = loadStations(OFd);
pos_stat = [y_stat, z_stat, x_stat]; % try to do something to stick with rest of code
istattab = 1:numel(x_stat); nstat = length(istattab); % Build istattab.
ystattab = pos_stat(istattab, 1); zstattab = pos_stat(istattab, 2); xstattab = pos_stat(istattab, 3);
if(0)
  % TODO: A switch for this test case.
  % Selection with respect to propagation angle. See [Brissaud et al.,
  % 2016, Section 6.1].
%   NBV             = sqrt(NSQ(1))*10;
  NBV           = sqrt(NSQ(1));
  beta          = acos( 2*pi*(1/T_0)/NBV );
  beta_stat_tot = atan(pos_stat(:, 2) ./ (pos_stat(:,1)-X_0));
  beta_stat     = beta_stat_tot(istattab);
  istattab_ok   = istattab .* (beta_stat<beta)';
end
% Read synthetics' files.
for locStatNumber = 1:nstat
    gloStatNumber = istattab(locStatNumber);
    [tmpDat, nt] = readAndSubsampleSynth(OFd, gloStatNumber, 'BXX', ['sem', signal_type], subsample, subsample_dt, locStatNumber);
    Xtime(locStatNumber,1:nt) = tmpDat(:, 1)';
    Xamp(locStatNumber,1:nt) = tmpDat(:, 2)';
    [tmpDat, nt] = readAndSubsampleSynth(OFd, gloStatNumber, 'BXZ', ['sem', signal_type], subsample, subsample_dt, locStatNumber);
    Ztime(locStatNumber, 1:nt) = tmpDat(:, 1)';
    Zamp(locStatNumber, 1:nt) = tmpDat(:, 2)';
end
clear('tmpDat');
subsampledDt = Ztime(locStatNumber, 2)-Ztime(locStatNumber, 1);
nSamp = nt;
%   dt = subsampledDt; % update dt
disp(['[',mfilename,'] ',num2str(size(Ztime,1)),' synthetics loaded: ', num2str(nSamp), ' samples each.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of the analytical solution.
disp([' ']); disp(['[',mfilename,', INFO] Building analytical solution.']);
mult_tSpan = 1; % multiplier for time span (used for analytical solution), must be whole numbers
mult_xSpan = 1; % multiplier for x span (used for analytical solution), must be whole numbers
mult_ySpan = 1; % multiplier for y span (used for analytical solution), must be whole numbers, if 2D simulation leave this to 1
dt_anal = min(diff(Ztime(1,:))); % Choose a dt corresponding to the dt of the simulation.
[t, x, y, TFM0, KX, KZ] = FK_buildAnalytical(t_0, dt_anal, ...
                                             xmin, xmax, ymin, ymax, dx, dy, ...
                                             zstattab, ...
                                             externalDGAtmosModel, USE_ISOTHERMAL_MODEL, SOUNDSPEED, H, GRA, NSQ, wind_x, wind_y, GAM, ...
                                             TYPE_FORCING, T_0, FORCING_FACTOR, mult_tSpan, mult_xSpan, mult_ySpan);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Brissaud et al., 2016, Section 5.1]: "Multiplication, in the Fourier
% domain, of the forcing function with a complex filter based on the
% representation of the solution in the case of an harmonic source or
% forcing term (see Appendix B for more details)."

%         z_station = zstattab(istat);
%         shifttime=0.0-real(KZ*z_station)./Omega;
%         KZ(1,1,:) = 0;
% remove the waves moving with a shifttime < 0
%         filt=exp(i*(KZ*z_station));
%       filt=exp(i*(KZ*z_station)).*(shifttime>=0);
%         filt(1,1,:)=0.0;
% remove waves assumed to propagate downward
%       indw=find((real(KZ)>0.0));
%      indw=find((real(KZ)<0.0));
%       filt(indw)=0.0;
      % try to remove aliasing
%         ind1=find(shifttime>(maxwind/2));
%        ind=find(abs(shifttime)>(maxwind/2));
%        filt(ind1)=exp(-abs(shifttime(ind1))/(maxwind/2)).*filt(ind1);
%         Mz=1*exp(z_station/(2*H))*ifftn((filt.*TFMo));
% compute X (horizontal component) assuming UX(f)=-(KZ/KX)*UZ
% equations 2.8 and 2.16 of Nappo, 2002 (not affected by wind)
% but frequencies close to (wi*KX-Omega)=0 create very large KZ...
%         TFMx=TFMo.*(KZ./KX).*filt;
%         ind1=find(isnan(TFMx));
%         ind2=find(isinf(TFMx));
%         TFMx(ind1)=0.0+i*0.0;
%         TFMX(ind2)=0.0+i*0.0;
%         Mx=0.0-exp(z_station/(2*H))*ifft2(TFMx);
%        xstattab_analytic = size(xstattab);
%        ystattab_analytic = size(ystattab);
synf = [];
for locStatNumber = 1:nstat
  gloStatNumber = istattab(locStatNumber);
  x_station = xstattab(gloStatNumber);
  z_station = zstattab(gloStatNumber);

  % MAYBE THIS FILTER TO TO BRING BACK TO ACTUAL TIME SERIES, would make sense
  %filt = exp(1i*(KZ*z_station));
  if(externalDGAtmosModel)
    error('no filter implemented for external atmos model');
  else
    if(USE_ISOTHERMAL_MODEL)
      % isothermal
%       filt = exp(1i*(KX*x_station + KZ*z_station));
%       filt = exp(1i*(KX*x_station + KZ*z_station) + z_station/H);
      filt = exp(1i*(KX*x_station + KZ*z_station) + z_station/(2*H));
%           filt = exp(1i*(KX*x_station + KZ*z_station) - 2*z_station/H);
%           filt = exp(1i*(KX*x_station + KZ*z_station) - z_station/(2*H));
%           filt = exp(1i*(KX*x_station + KZ*z_station) + z_station/H + 2*z_station/(SOUNDSPEED^2));
%           filt = exp(1i*(KX*x_station + KZ*z_station) - z_station/H);
    else
      % isobaric
      filt = exp(1i*(KX*x_station + KZ*z_station));
    end
  end
%       filt = exp(1i*(KX*x_station + KZ*z_station));
%       filt = exp(1i*(KX*x_station + KZ*z_station) + z_station/H);
%       shifttime=0.0-real(KZ*z_station)./Omega;
%       filt = exp(1i*KZ*z_station).*(shifttime>=0);
  % what
%       filt(1, 1, :) = 0.0; % why
  % remove waves assumed to propagate downward
%       indw=find((real(KZ)>0.0));
%       indw=find((real(KZ)<0.0));
%       filt(indw)=0.0;
%       Mz = exp(-z_station/(2*H)) * ifftn(filt.*TFMo);
%       Mz = exp(z_station/(2*H)) * ifftn(filt.*TFMo);
%       Mz = exp(z_station/(H)) * ifftn(filt.*TFMo);
%       Mz = 0.5*exp(z_station/(2*H)) * ifftn(filt.*TFMo);
  %Mz = exp(z_station/69109) * ifftn(filt.*TFMo);
%       Mz = exp(z_station/68719) * ifftn(filt.*TFMo);
  Mz = ifftn(filt.*TFM0);
%       Mz = ifftn(filt.*TFMo) * exp(2*z_station/SOUNDSPEED);
%           xcoord(istat) = xstattab(istat);
  ix = round((xstattab(gloStatNumber)-xmin)/dx) + 1; % first guess for x
  iy = round((ystattab(gloStatNumber)-zmin)/dy) + 1; % first guess for y
%           iy=round((ystattab(istat)-ymin)/dy) + 1;
  % GET X LOCATION OF STATION
%       ACHECK = xstattab(gloStatNumber); positionInMeshgrid = X(1, ix, 1); positionInMeshgrid_p1 = X(1, ix+1, 1); positionInMeshgrid_m1 = X(1, ix-1, 1);
  toCheck = xstattab(gloStatNumber); positionInMeshgrid = x(ix); positionInMeshgrid_p1 = x(ix+1); positionInMeshgrid_m1 = x(ix-1);
  if(abs(positionInMeshgrid-toCheck) > abs(positionInMeshgrid_p1-toCheck))
    ix = ix + 1; % if real x is closer to ix+1 than to ix, choose this ix+1
  elseif(abs(positionInMeshgrid-toCheck) > abs(positionInMeshgrid_m1-toCheck));
    ix = ix - 1; % or if real x is closer to ix-1 than to ix, choose this ix-1
  end  
  xstattab_analytic(gloStatNumber) = x(ix); %X(1, ix, 1);

  % GET Y LOCATION OF STATION
  if(numel(y)>1)
    % This if is to prevent trying to find another possible y when
    % simulation is 2D (because in the 2D case, the array 'y' only
    % contains one element, 0).
%       ACHECK = ystattab(gloStatNumber); positionInMeshgrid = Y(iy, 1, 1); positionInMeshgrid_p1 = Y(iy+1, 1, 1); positionInMeshgrid_m1 = Y(iy-1, 1, 1);
    toCheck = ystattab(gloStatNumber); positionInMeshgrid = y(iy); positionInMeshgrid_p1 = y(iy+1); positionInMeshgrid_m1 = y(iy-1);
    if(abs(positionInMeshgrid-toCheck) > abs(positionInMeshgrid_p1-toCheck))
      iy = iy + 1; % if real y is closer to iy+1 than to iy, choose this iy+1
    elseif(abs(positionInMeshgrid-toCheck) > abs(positionInMeshgrid_m1-toCheck))
      iy = iy - 1; % or if real y is closer to iy-1 than to iy, choose this iy-1
    end
  end
  ystattab_analytic(gloStatNumber) = y(iy); %Y(iy,1,1)

  % If out of mesh, bring back to closest mesh point.
  if(ix<=0); ix = 1; end
  if(iy<=0); iy = 1; end

  % Take real part as analytic solution.
  synf(locStatNumber, :) = real(Mz(iy, ix, :));
  
  if(z_station>0)
    % Remove constant offset (the analytical solution should start valued 0).
    synf(locStatNumber, :) = synf(locStatNumber, :)-synf(locStatNumber, 1);
  end
end

% test to see what forms has the analytic solution.
disp([' ']); disp(['[',mfilename,', INFO] Checking analytical solution.']);
FK_checkAnalytical;

% TEST ONLY
% synf = Zamp; 
% test only !!!!
% synf_X=synf;
% clear TFMx;
% clear Mx ;
% clear Mz ;
timef = t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finalise.
factor_err = 1;
allStations = 1:nstat;
stationsToPlot = allStations;
nstatToPlot = numel(stationsToPlot);

disp([' ']);
disp(['[',mfilename,'] Computing errors.']);
FK_computeError;

disp([' ']);
disp(['[',mfilename,'] Plotting everything.']);
FK_plot;