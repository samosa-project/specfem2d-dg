clear all;
close all;
clc;

debugfig = 0;

station_ids = [1,2];
setupLocal;

time_stf_slant_IP = 0.059000000;
time_stf_slant_RP = 0.082800001;
time_stf_slant_RS = 0.097199999;
time_stf_slant_TP = 0.097599998;

time_stf_ortho_IP = 0.024800001;
time_stf_ortho_RP = 0.049800001;
time_stf_ortho_TP = 0.064400002;

time_fts_slant_IP = 0.128399998;
time_fts_slant_TP = 0.172600001;
time_fts_slant_RP = 0.180000007;
time_fts_slant_TS = 0.203999996;

time_fts_ortho_IP = 0.054000001;
time_fts_ortho_TP = 0.093400001;
time_fts_ortho_RP = 0.108000003;

stf_slant_the = zeros(3, numel(nvals)); % TPP, RPP, RPS
stf_slant_res = zeros(3, numel(nvals)); % TPP, RPP, RPS
stf_ortho_the = zeros(3, numel(nvals)); % TPP, RPP, RPS=nan
stf_ortho_res = zeros(3, numel(nvals)); % TPP, RPP, RPS=nan
fts_slant_the = zeros(3, numel(nvals)); % TPP, TPS, RP
fts_slant_res = zeros(3, numel(nvals)); % TPP, TPS, RP
fts_ortho_the = zeros(3, numel(nvals)); % TPP, TPS=nan, RP
fts_ortho_res = zeros(3, numel(nvals)); % TPP, TPS=nan, RP

for c = 1:numel(folderz)
% for c = 4
  folder = folderz{c};
  nel = nelz{c};
  ival = find(nel==nvals);
  icas = idcasez{c};
  disp(['[',mfilename,'] Considering folder ''',regexprep(folder,thisFolder,''),'''.']);
  
  OFD=dir([folder, filesep, 'OUTPUT_FILES_*']);
  switch(numel(OFD))
    case 0
      error(['found no output files folder']);
    case 1
      % ok
    otherwise
      error(['more than one output files folder']);
  end
  OFD = [OFD.folder,filesep,OFD.name,filesep];
  parfile = [OFD, 'input_parfile'];
  [t, amp] = load_synthetics(OFD, parfile, station_ids);
  t = t(1,:);
  dp = squeeze(amp(1, 1, :));
  vx = squeeze(amp(1, 2, :));
  vz = squeeze(amp(2, 2, :));
  [~, vr] = cart2pol(vx, vz);
  
  [rho__1, alpha__1, rho__2, alpha__2, beta__2] = get_models(parfile);
  Z1 = alpha__1*rho__1;
  [i1_i2_j2t_j2r] = get_predicted_angles_deg(ic_rad, alpha__1, alpha__2, beta__2)* pi/180; % cast to [rad]
  
  if(cases{icas}.fts0_stf1)
    if(cases{icas}.ortho0_slant1)
      i=mt(t,time_stf_slant_TP,dp); r_TP = dp(i);
      i=mt(t,time_stf_slant_IP,vr); [t_IP, r_IP] = cart2pol(vx(i), vz(i));
      i=mt(t,time_stf_slant_RP,vr); [t_RP, r_RP] = cart2pol(vx(i), vz(i));
      i=mt(t,time_stf_slant_RS,vr); [t_RS, r_RS] = cart2pol(vx(i), vz(i));
      if(debugfig)
        figure(); plot(t, dp*range(vz)/range(dp)); hold on; plot(t, vx); plot(t, vz);
        scatter(t(dp==r_TP), r_TP*range(vz)/range(dp), 'x'); scatter(t(vr==r_IP), vz(vr==r_IP), 'x'); scatter(t(vr==r_RP), vz(vr==r_RP), 'x'); scatter(t(vr==r_RS), vz(vr==r_RS), 'x');
      end
      
      t_IP = pi/2-t_IP; % measure it clockwise so remove it from 90 deg
      t_RP = -(pi/2-t_RP); % measure it from the other side so flip sign
      t_RS = -(pi/2 - pi/2-t_RS); % because S wave, displacement direction is 90 deg from propagation direction, plus measure it from other side
      
      [R, T] = RTCoefs(cases{icas}.fts0_stf1, alpha__1, rho__1, alpha__2, beta__2, rho__2, ic_rad);
      stf_slant_the(:, ival) = [T, abs(R)];
      stf_slant_res(:, ival) = [r_TP/Z1, r_RP, r_RS]/r_IP; % divide by impedance to get velocity in air
      
    else
      i=mt(t,time_stf_ortho_TP,dp); r_TP = dp(i);
      i=mt(t,time_stf_ortho_IP,vr); [t_IP, r_IP] = cart2pol(vx(i), vz(i));
      i=mt(t,time_stf_ortho_RP,vr); [t_RP, r_RP] = cart2pol(vx(i), vz(i));
      if(debugfig)
        figure(); plot(t, dp*range(vz)/range(dp)); hold on; plot(t, vx); plot(t, vz);
        scatter(t(dp==r_TP), r_TP*range(vz)/range(dp), 'x'); scatter(t(vr==r_IP), vz(vr==r_IP), 'x'); scatter(t(vr==r_RP), vz(vr==r_RP), 'x');
      end
      
      t_IP = pi/2-t_IP; % measure it clockwise so remove it from 90 deg
      t_RP = -(pi/2-t_RP); % measure it from the other side so flip sign
      
      [R, T] = RTCoefs(cases{icas}.fts0_stf1, alpha__1, rho__1, alpha__2, beta__2, rho__2, 0);
      stf_ortho_the(:, ival) = [T, abs(R)]; 
      stf_ortho_res(:, ival) = [r_TP/Z1, r_RP, nan]/r_IP; % divide by impedance to get velocity in air
      
    end
  else
    if(cases{icas}.ortho0_slant1)
      i=mt(t,time_fts_slant_IP,dp); r_IP = dp(i);
      i=mt(t,time_fts_slant_RP,dp); r_RP = dp(i);
      i=mt(t,time_fts_slant_TP,vr); [t_TP, r_TP] = cart2pol(vx(i), vz(i));
      i=mt(t,time_fts_slant_TS,vr); [t_TS, r_TS] = cart2pol(vx(i), vz(i));
      if(debugfig)
        figure(); plot(t, dp*range(vz)/range(dp)); hold on; plot(t, vx); plot(t, vz);
        scatter(t(dp==r_IP), r_IP*range(vz)/range(dp), 'x'); scatter(t(dp==r_RP), r_RP*range(vz)/range(dp), 'x'); scatter(t(vr==r_TP), vz(vr==r_TP), 'x'); scatter(t(vr==r_TS), vz(vr==r_TS), 'x');
      end
      
      t_TP = pi/2-t_TP; % measure it clockwise so remove it from 90 deg
      t_TS = pi/2+t_TS+pi/2; % because S wave, displacement direction is 90 deg from propagation direction, plus measure it from other side
      
      [R, T] = RTCoefs(cases{icas}.fts0_stf1, alpha__1, rho__1, alpha__2, beta__2, rho__2, ic_rad);
      fts_slant_the(:, ival) = [abs(T), R]; 
      fts_slant_res(:, ival) = [r_TP, r_TS, r_RP/Z1]/(r_IP/Z1); % divide by impedance to get velocity in air
      
    else
      i=mt(t,time_fts_ortho_IP,dp); r_IP = dp(i);
      i=mt(t,time_fts_ortho_RP,dp); r_RP = dp(i);
      i=mt(t,time_fts_ortho_TP,vr); [t_TP, r_TP] = cart2pol(vx(i), vz(i));
      if(debugfig)
        figure(); plot(t, dp*range(vz)/range(dp)); hold on; plot(t, vx); plot(t, vz);
        scatter(t(dp==r_IP), r_IP*range(vz)/range(dp), 'x'); scatter(t(dp==r_RP), r_RP*range(vz)/range(dp), 'x'); scatter(t(vr==r_TP), vz(vr==r_TP), 'x');
      end
      
      t_TP = pi/2+t_TP; % measure it clockwise so remove it from 90 deg
      
      [R, T] = RTCoefs(cases{icas}.fts0_stf1, alpha__1, rho__1, alpha__2, beta__2, rho__2, 0);
      fts_ortho_the(:, ival) = [abs(T), R]; 
      fts_ortho_res(:, ival) = [r_TP, nan, r_RP/Z1]/(r_IP/Z1); % divide by impedance to get velocity in air
      
    end
  end
      
  if(nel==max(nvals))
    ff=diag_time_series(icas, t, dp, vx, vz, R, T, Z1, r_IP, i1_i2_j2t_j2r);
    customSaveFig(ff, [plotFolder,filesep,cases{icas}.code], extToSave, 9999);
  end
  
end

sf_slant_err = 100* abs(stf_slant_res-stf_slant_the)./stf_slant_the;
sf_ortho_err = 100* abs(stf_ortho_res-stf_ortho_the)./stf_ortho_the;
fs_slant_err = 100* abs(fts_slant_res-fts_slant_the)./fts_slant_the;
fs_ortho_err = 100* abs(fts_ortho_res-fts_ortho_the)./fts_ortho_the;

[fig] = diag_errors(sf_ortho_err, sf_slant_err, fs_ortho_err, fs_slant_err);
customSaveFig(fig, [plotFolder,filesep,'summaryFluidSolidTransmission'], {'fig','tex'}, 9999);

function mintid = mt(t, t0, y)
%   mintid = find(abs(t-t0)==min(abs(t-t0)));
  mintid = find(y==max(y(abs(t-t0)<0.005)));
end