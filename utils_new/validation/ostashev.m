% Author:        Léo Martire.
% Description:   Validates 2 simulations against [Ostashev et al., 2005].
% Notes:         [Ostashev et al., 2005] Ostashev, V. E., Wilson, D. K.,
%                  Liu, L., Aldridge, D. F., Symons, N. P., and Marlin, D.
%                  (2005).  Equations for finite-difference, time-domain
%                  simulation of sound propagation in moving inhomogeneous
%                  media and numerical implementation. The Journal of the
%                  Acoustical Society of America, 117(2):503–517.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

clear all;
% close all;
clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'defaultTextFontSize', 16); set(0, 'defaultAxesFontSize', 16);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools'); % truncToShortest, readAndSubsampleSynth
addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/standalone');
SPCFMEXloc = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/';
rootd = strcat(SPCFMEXloc,'validation_lns_ostashev/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate stations or plot results?
plotResults0orgenstats1orgengmsh2=0;

% params for plotting restulsts
NSTATIONS=15;

% ofdweq0='OUTPUT_FILES_rho_M0_dx1_FNS/'; ofdwneq0=ofdweq0; wind=0; simullab='FNS $\rho$'; % pas ouf
% ofdweq0='OUTPUT_FILES_rho_M0_dx1_FNS/'; ofdwneq0='OUTPUT_FILES_rho_M.3_dx1_FNS/'; wind=102; simullab='FNS $\rho$'; % pas ouf
% ofdweq0='OUTPUT_FILES_E_M0_dx1_FNS/'; ofdwneq0=ofdweq0; wind=0; simullab='FNS $E$';
% ofdweq0='OUTPUT_FILES_E_M0_dx1_FNS/'; ofdwneq0='OUTPUT_FILES_E_M.3_dx1_FNS/'; wind=102; simullab='FNS $E$';

% ofdweq0='OUTPUT_FILES_rho_M0_dx1/'; ofdwneq0=ofdweq0; wind=0; simullab='LNS $\rho$'; % pas ouf
% ofdweq0='OUTPUT_FILES_rho_M0_dx1/'; ofdwneq0='OUTPUT_FILES_rho_M.15_dx1/'; wind=51; simullab='LNS $\rho$'; % pas ouf
% ofdweq0='OUTPUT_FILES_rho_M0_dx1/'; ofdwneq0='OUTPUT_FILES_rho_M.3_dx1/'; wind=102; simullab='LNS $\rho$'; % pas ouf
% ofdweq0='OUTPUT_FILES_rho_M0_dx1/'; ofdwneq0='OUTPUT_FILES_rho_M.3_dx1_testformulaE/'; wind=102; simullab='LNS $\rho$ test $E$'; % pas ouf
% ofdweq0='OUTPUT_FILES_rho_M0_dx1_testformulaM/'; ofdwneq0='OUTPUT_FILES_rho_M.3_dx1_testformulaM/'; wind=102; simullab='LNS $\rho$ test $\rho v$'; % pas ouf
% ofdweq0='OUTPUT_FILES_rho_M0_dx.5/'; ofdwneq0='OUTPUT_FILES_rho_M.3_dx.5/'; wind=102; simullab='LNS $\rho$ $\Delta x=.5$'; % pas ouf
% ofdweq0='OUTPUT_FILES_rho_M0_dx1gmsh/'; ofdwneq0=ofdweq0; wind=0; simullab='LNS GMSH $\rho$'; % alright
% ofdweq0='OUTPUT_FILES_rho_M0_dx1gmsh/'; ofdwneq0='OUTPUT_FILES_rho_M.3_dx1gmsh/'; wind=102; simullab='LNS GMSH $\rho$'; % pas ouf

% ofdweq0='OUTPUT_FILES_rho_M0_dx1_RK6/'; ofdwneq0=ofdweq0; wind=0; simullab='LNS $\rho$ 6RK4'; % pas ouf
% ofdweq0='OUTPUT_FILES_rho_M0_dx1_RK6/'; ofdwneq0='OUTPUT_FILES_rho_M.3_dx1_RK6/'; wind=102; simullab='LNS $\rho$ 6RK4'; % pas ouf

% ofdweq0='OUTPUT_FILES_rho_M0_dx1_fullE/'; ofdwneq0=ofdweq0; wind=0; simullab='LNS $\rho$ fullE'; % pas ouf
% ofdweq0='OUTPUT_FILES_rho_M0_dx1_fullE/'; ofdwneq0='OUTPUT_FILES_rho_M.3_dx1_fullE/'; wind=102; simullab='LNS $\rho$ fullE'; % pas ouf

% ofdweq0='OUTPUT_FILES_E_M0_dx1/'; ofdwneq0=ofdweq0; wind=0; simullab='LNS $E$';
ofdweq0='OUTPUT_FILES_E_M0_dx1/'; ofdwneq0='OUTPUT_FILES_E_M.3_dx1/'; wind=102; simullab='LNS $E$';
% ofdweq0='OUTPUT_FILES_E_M0_dx1_reredone/'; ofdwneq0='OUTPUT_FILES_E_M.3_dx1_reredone/'; wind=102; simullab='LNS $E$ reredone'; % did not change anything

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% treatment %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f0=100;
c=340;
M=wind/c;
% According to Ostashev D:
% dx=dy=dr<lambda*(1-M)/Npts=c*(1-M)/(Npts*f0)
% dt<(1-M)/((1+M)*Npts*f0)
Npts=2;
dxmax=c*(1-M)/(Npts*f0);
dtmax=(1-M)/((1+M)*Npts*f0);

kr=20;
alpha=linspace(0,pi,NSTATIONS);
k=2*pi*f0/c;
r=kr/k;
% A=1;
% p_th= A*(sqrt(1-M^2*sin(alpha).^2) - M*cos(alpha)) ./ (sqrt(2*pi*k*r) * (1-M^2) * (1-M^2*sin(alpha).^2).^0.75);
% p0_th= A ./ (sqrt(2*pi*k*r));
p_th = p(k, r, alpha, M)
p0_th = p(k, r, 0, 0)

switch(plotResults0orgenstats1orgengmsh2)
  case (0)
    figPath = [rootd,'comp__',regexprep(ofdweq0,'/',''),'__vs__',regexprep(ofdwneq0,'/','')];
    % load results
    OFd_w0 = strcat(rootd, ofdweq0);
    OFd_w = strcat(rootd, ofdwneq0);
    % load reference
    [data, ~] = readAndSubsampleSynth(OFd_w0, 1, 'BXZ', 'semv', 0, -1, 1); Zamp0(1,:)=data(:,2)';
    % load test stations
    ZampW=zeros(NSTATIONS,size(Zamp0,2));
    for istat = 1:NSTATIONS
      [data, ~] = readAndSubsampleSynth(OFd_w, istat, 'BXZ', 'semv', 0, -1, istat); ZampW(istat,:)=data(:,2)';
    end
    [xstattab, ystattab, stations_data] = loadStations(OFd_w); % Load stations data (first try OUTPUT folder, then if not found, try parent DATA folder).
    
    % get reference
    p0_exp=Zamp0; % Use this after having ran synth_load on Mach=0 simulation and selected the right station.
    p0_exp_ampl=max(p0_exp)-min(p0_exp);
    disp(['[',mfilename,'] Reference amplitude: ',num2str(p0_exp_ampl),' [Pa].']);
    % get normalised
    p_exp=ZampW;
    p_exp_ampl=(max(p_exp')-min(p_exp'));
    angglle=atan(ystattab./xstattab)*180/pi;
    angglle(angglle<0)=angglle(angglle<0)+180;
    angglle(end)=180;
    disp(['[',mfilename,'] Angles and windy (w=',num2str(wind),' [m/s]) amplitudes [Pa]:']);
    disp([angglle';p_exp_ampl]);
    
    % actuallay plot
    fighandlll = figure('outerposition',[0 0 750 750]);
    
    plot(alpha*180/pi,p_th/p0_th,'k','displayname','Ostashev''s (72)'); hold on
    
    maxdiffpercent=max(abs(p_th/p0_th - p_exp_ampl/p0_exp_ampl)./(p_th/p0_th))*100;
    
    plot(angglle,p_exp_ampl/p0_exp_ampl,'displayname',[simullab, ' (max. diff. = ',sprintf('%.3g',maxdiffpercent),'\%)']);
    
%     set(gca,'ticklabelinterpreter','latex');
    set(gca,'xtick',180*[0,0.25,0.5,0.75,1]);
    xlim([0,180]);
    title(['$f_0=',num2str(f0),'$ Hz, $c=',num2str(c),'$, $w=',num2str(wind),'$, $M=',num2str(M),'$, $k=',num2str(k),'$, $r=',num2str(r),'$, $kr=',num2str(k*r),'$']);
    ylabel(['normalised sound pressure amplitude']);
    xlabel(['azimuth [deg], counter-clockwise from wind vector']);
    legend('location','northwest');
    prettyAxes(fighandlll);
    customSaveFig(figPath,{'jpg','eps'});
    
    % manual plot
    if(0)
      p0=Zamp; % Use this after having ran synth_load on Mach=0 simulation and selected the right station.
      p0_ampl=max(p0)-min(p0)

      p=Zamp; % Use this after having ran synth_load on relevant simulation.
      p_ampl=(max(p')-min(p'))
      angglle=atan(ystattab./xstattab)*180/pi;
      angglle(angglle<0)=angglle(angglle<0)+180;
      angglle(end)=180;

      plot(angglle,p_ampl/p0_ampl,'displayname','LNS');
      plot(angglle,p_ampl/p0_ampl,'displayname','FNS');
    end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case (1)
    clc;
    disp(['nreceiversets = ',num2str(numel(alpha))]);
    disp(['# Orientation']);
    disp(['anglerec                        = 0.0d0          # Angle to rotate components at receivers.']);
    disp(['rec_normal_to_surface           = .false.        # Base anglerec normal to surface (external mesh and curve file needed).']);
    for ang=alpha
      disp(['# az ', num2str(ang*180/pi),'°']);
      disp(['nrec = 1']);
      disp(['xdeb = ',char(sprintf('%.16g',r*cos(ang)))]);
      disp(['zdeb = ',char(sprintf('%.16g',r*sin(ang)))]);
      disp(['xfin = 0.']);
      disp(['zfin = 0.']);
      disp(['record_at_surface_same_vertical = .false.']);
    end
    disp(' ');
    disp(['[',mfilename,'] Check if dx<',num2str(dxmax),' and dt<',num2str(dtmax),'']);
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case (2)
    % GMSH point definition
    for angi=1:length(alpha)
      disp(['Point(',num2str(4+angi),') = {',char(sprintf('%.18g',r*cos(alpha(angi)))),',',char(sprintf('%.18g',r*sin(alpha(angi)))),',0};']);
    end
    for angi=1:length(alpha)
      disp(['Point{',num2str(4+angi),'} In Surface{1};']);
    end
    
  otherwise
    error('kek')
end

function P=p(k, r, alp, mach)
  A = 1;
  P = A * (sqrt(1-mach^2*sin(alp).^2) - mach*cos(alp)) ./ (sqrt(2*pi*k*r) * (1-mach^2) * (1-mach^2*sin(alp).^2).^0.75);
end