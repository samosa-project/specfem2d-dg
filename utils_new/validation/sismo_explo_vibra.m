%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                       seismoplot4.m                             %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use this file to plot on seismograms traces exiting from specfem2D in 
% ASCII filesand a synthetic signal waited at the station in case of 
% simpler homogeneous medium with a variable density and a constant 
% velocity.
% You will choose the directory in which the ASCII files are.
% Then choose the number of receivers with the velocity of the signal and 
% finally run this m-file.
% This program is made for a maximum number of 9 stations.

clear all ; clc ;
% close all;
format long ;
                        
%%%%%%%%%%%%%%%%%%%%r(istat,1:NF
% Kind of simulation
simu = 9;               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collect simulation parameters from par_file
seismotype = 1 ;        % write the number chosen in the Par_file for the 
                        % seismotype.
                        
%%%%%%%%%%%%%%%%%%%%%%%
% Compare simulations ?
compare = false;          
% compare = false; 
comparison_mod = false;



% Parameters which have to be chosen
dt = 1*10^(-2) ;
if(compare)
   dt_2 = 1*10^(-2) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequency domain initialization
fs = 1/dt ;
fmax = 1 ;
fmin = -fmax ;

% Amplitude (to be removed)
A = 1.0;

%%%%%%%%%%%%%%%%%%%%
% Density parameters
rho_cst      = 1.2 ;
% Impose constant density or not
rho_cst_bool = false;

MODEL_QUENTIN = false;

%%%%%%%%%%%%%%%
% Viscosity eta
eta_cst = 1;
% eta_cst = (4/3)*eta_cst;
% Impose constant viscosity or not
eta_cst_bool = true;
% Is there attenuation in the external file ?
any_eta = false;

%%%%%%%%%%%%%%%%
% Sound velocity
speed_cst_bool = false;
% Value
speed_cst = 342;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if only analytic max time iteration
nt_max = 2000;

%%%%%%%%%%%%%%%%%%%%%%%
% Forcing signal period
perio = 1;
to    = 1.2;%-0.06;

%%%%%%%%%%%%%%%%%%%%
% Plot eta/rho ratio
plot_ratio_etarho = false;

plot_stability = false;

%%%%%%%%%%%%%%%%
% Plot CFL conf.
plot_CFL = false;
plot_cfl_vit = false;
freq_CFL = 1/60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot aborption coefficient alpha
plot_aborb = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot inverse quality factor (% to absorption)
quality_factor = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity data in external model
any_gravity = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure data in external model
any_pressure = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relaxation times
tau_eps = 1%1.12;
tau_sig = 1%0.016;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External file containing model data
%external_file = 'msise_11032011_model10821.txt';
% external_file = '1D_isothermal_thermosphere_model.txt';
 external_file = 'atmos_model/1D_iso_enhanced.txt';
% external_file = '1D_isothermal_atmosphere_model_N2const.txt';
% external_file = '1D_isothermal_thermosphere_model_N2const.txt';
% external_file = 'MSISE_attenuation.txt';
% external_file = '1D_isothermal_atmosphere_model.txt';
% external_file = 'MSISE_attenuation_enhanced.txt';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation files (in Documents/ directory)
%   simu_files = 'Results.semd.FDTD';
% simu_files = 'Results.semd.gravi_para16';
% simu_files = '~/Documents/FD/FD_21';
simu_files = '~/Documents/GIT_3D_withinversion/FD_3';
% simu_files = 'OLD_ACOUS_ATTEN/Results.FD_2';
% simu_files = 'GJI_PAPER/1Diso_rhovar_grav_noatten_nowind_acousForc_RK4';
% simu_files = 'GJI_PAPER/1Diso_rhovar_grav_noatten_nowind_acousForc_RK4';
% simu_files = 'GJI_PAPER/1Diso_rhovar_nograv_noatten_nowind_acousForc_RK4';

if(compare)
% simu_files_2 = 'Results.FD_2(copy).save';
simu_files_2 = 'Results.FD_2';
% simu_files_2 = 'GJI_PAPER/1Diso_rhocte_nograv_noatten_nowind_acousForc';
end

%%%%%%%%%%%%%%%%%%
% Graphs requested
analytic   = true;
modeled    = true;
comparison = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if we want y-axis restricted to [-1, 1]
axis_bounded = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if we want to plot phase velocity
plot_phase = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legends for simulation graphs
if(seismotype == 1)
    variable = ' Displacement along' ;
    signal_type = 'd' ; % -> 'a' for acceleration
                        % -> 'd' for displacement or sqrt(density) *
                        %    displacement
                        % -> 'v' for velocity or sqrt(density) * velocity
                        % Watch the right letter at the end of the output
                        % ascii file from specfem2D
else if(seismotype == 2)
        variable = 'Velocity along' ;
        signal_type = 'v' ;
else if (seismotype == 3)
        variable = 'Acceleration along' ;
        signal_type = 'a' ;
else if (seismotype == 4)
        variable = 'Pressure along ' ;
        signal_type = 'p' ;
else if (seismotype == 5)
        variable = 'Curl of displacement along ' ;
        signal_type = 'c' ;
else if (seismotype == 6)
        variable = 'Fluid potential along ' ;
        signal_type = 'p' ;                  
else if (seismotype == 7)
        variable = '\rho^{1/2} cdot {\bf u} along' ;
        signal_type = 'd' ;
    end
    end
    end
    end
    end
    end
end
% signal_type = 'p' ;
%%%%%%%%%%%%%%%%%%%%%%
% Load simulation data
directory = strcat(simu_files) ;
if(compare)
    directory_2 = strcat('../',simu_files_2,'/') ;
end

%%%%%%%%%%%%%%%%%%%%%%
% Stations data
pos_stat = load(strcat(simu_files,'/stations.txt')) ;

istattab = [3];

xstattab = pos_stat(istattab,1);
zstattab = pos_stat(istattab,2);
nstat    = size(xstattab,1) ;
% istattab = 1:size(xstattab,1) ;

% zstattab(1) = zstattab(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gathering external model data
atmos = load(strcat(...
    simu_files,'/',external_file)) ;

% Id of parameter in file
k = 1;

% initialization
variable_velocity = false; 

alti_atm = atmos(:,k) ;
k = k+1;

alti_diff = max(abs(alti_atm(1:end-1) - alti_atm(2:end)));

if(rho_cst_bool)
    rho = ones(size(atmos(:,k)))*rho_cst ;
else
    rho = atmos(:,k);
end
rho_0 = rho(1) ;
k = k+1;

%%% MODIF A SUP
if(MODEL_QUENTIN)
rho = 1.2*exp(-(abs(9.81))*alti_atm/(speed_cst^2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If eta among the external file data
if(eta_cst_bool)
    eta = eta_cst*ones(size(atmos(:,1)));
    if(any_eta)
      k = k+1;
    end
else
    eta = atmos(:, k);
    k = k+1;
end

velocity = atmos(:,k) ;
k = k+1;
V = velocity;%3.42857143e+02 ;       % must be in m/s
% Constant velocity ?
if(speed_cst_bool)
    V = speed_cst*ones(size(atmos(:,1)));
    velocity = V;
end
% Check if there's different velocities
if(find(V ~= V(1)))
   variable_velocity = true; 
end

%%%%%%%%%%%%%%%%%%%%%%
% If gravity influence
if(any_gravity)
    gravity = atmos(:,k) ;
    k = k+1;
end

%%%%%%%%%%%%%%%%%%%%%%
% Brunt-Valsala number
Nsq = atmos(:,k) ;
k = k+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If pressure among the external file data
if(any_pressure)
    Pressure = atmos(:,k) ;
    k = k+1;
end

%%%%%%%%%%%%%
% Main titles
if(not(rho_cst_bool))
    rho_text = strcat('variable density profile based on " ',external_file,' " ');
else
    rho_text = strcat('constant density profile = ',num2str(rho_cst));
end
if(not(eta_cst_bool))
    eta_text = strcat('variable viscosity profile based on " ',external_file,' " ');
else
    eta_text = strcat('constant viscosity profile = ',num2str(eta_cst));
end
if(variable_velocity)
    veloc_text = 'variable velocity';
else
    veloc_text = 'velocity equal to ',num2str(velocity(1)),' m.s^{-1}';
end
subtitle = strcat('Forcing with the whole surface moving at the bottom of a fluid medium');
subtitle2 = strcat([' with a ',rho_text,', a ',eta_text,' and a ', veloc_text ]) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of ratio eta/rho (density on viscosity)
if(plot_ratio_etarho)
    
figure;
plot(alti_atm, eta./rho)
hold on
plot(alti_atm, 27296, '-r')

Max_etarho2 = max(eta./(rho.*rho))

x = find(eta./rho > 27296);
y = ylim;
% If there's a eta/rho ratio that overshoot stability limit
%if(any(x))
%    x = x(1)
%    line([alti_atm(x);alti_atm(x)],repmat(y(:),1,numel(x)),'color','r','linestyle','--')
%end
    
title(strcat('eta/rho, eta and rho based on " ',external_file,' ") against elevation.'))
xlabel('Elevation (m)','FontSize',14.3)
legend('\fontsize{13}{  eta/rho}','\fontsize{13}{  Numerical threshold of blow-up }','Location','West')

end %if(plot_ratio_etarho)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Profile of density with a constant scaling H
% Thermosphere model
%H = 6.568879214275763e+03;
% Atmosphere model
%H = 8.6611888899e+03;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read traces of simualtion
for istat = 1 : nstat
    
    istat_REAL = istattab(istat);
    
    display(istat_REAL)
    
%     coef = 0;
%     if(istat_REAL >= 4) 
%         coef = 250;
%     end
    
%     alti(istat) = ((istat-1)*Zpath+Z0)*1000 + coef ;
    alti(istat) = zstattab(istat);
    
    if(modeled)
    
        if(istat_REAL > 9)
            cmplt = '00';
        else
            cmplt = '000';
        end
        
    % reading of the horizontal component
%     file = strcat(directory,'S',cmplt,num2str(istat_REAL),'.AA.BXX.sem',signal_type) ;
%     file = strcat(directory,'S',cmplt,num2str(istat_REAL),'.AA.BXX.semp',signal_type) ;
%     data = load(file) ;
%     %if (istat == 1)
%         nt = max(size(data)) ;
%     %end
%     Xtime(istat,1:nt) = data(1:nt,1)' ;
%     Xamp(istat,1:nt) = data(1:nt,2)' ;
    
    % reading of the vertical component
%     file = strcat(directory,'/OUTPUT_SISMO_MARS_atten/AA.S',cmplt,num2str(istat_REAL),'.p.sem',signal_type) ;
%    file = strcat(directory,'/OUTPUT_SISMO_MARS_noatten_nowind/AA.S',cmplt,num2str(istat_REAL),'.p.sem',signal_type) ;
   file = strcat(directory,'/OUTPUT_SISMO/AA.S',cmplt,num2str(istat_REAL),'.p.sem',signal_type) ;
%     file = strcat(directory,'S',cmplt,num2str(istat_REAL),'.AA.PRE.sem',signal_type) 
    data = load(file) ;
    nt=max(size(data));
    Ztime_save(istat,1:nt) = data(1:nt,1)' ;
    Zamp(istat,1:nt) = data(1:nt,2)' ;
    Xtime(istat,1:nt) = data(1:nt,1)' ;
    Xamp(istat,1:nt) = data(1:nt,2)' ;
    
    mult = 10;
    dt = 0.05
%     temp = 0:dt:Ztime_save(istat,end)*5;%linspace(0,Ztime_save(istat,end)*mult,nt);
    
%     temp = [0:dt:1000];
    Ztime(istat,:) = Ztime_save(istat,:);
    
%     test(istat)=max(Zamp(istat,:))
%     Ztime(istat,:) = [0:dt:nt_max]' ;
    else %if modeled
   
    Ztime(istat,:) = [0:dt:nt_max]' ;
    nt = size(Ztime) ;
    
    end

end    


% nt = size(Ztime) ;


if(analytic)
    %%%%%%%%%%%%%%%%
    % Initialization
    syn  = zeros(nstat,length(Ztime(nstat,:))) ;
    
    NFFT  = 2*2^nextpow2(length(syn(nstat,:))); % Next power of 2 from length of syn(istat,:)
    
    
       k  = zeros(nstat,NFFT);
       f0 = zeros(nstat, 1);
       absorb = zeros(nstat,NFFT);
       filter = ones(nstat,NFFT);
    
    t0 = zeros(nstat,1) ;
    
    i = complex(0,1);
    
    
%    nstat=2
    
    for istat = 1:nstat
        
        display(istat)
        
        % Starting time
        t0(istat) = 0.0;
        
        % Vectorial indice of the station position
        ii_temp = find(abs(alti_atm - alti(istat)) < alti_diff);
        ii = min(ii_temp);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Computation of the signal
%         syn(istat,:) =  (exp(-((Ztime(istat,:)-(to-perio/4)-t0(istat))/(perio/4)).^2)...
%            - exp(-((Ztime(istat,:)-(to+perio/4)-t0(istat))/(perio/4)).^2));
       
        t = Ztime(istat,:);
        a = pi*pi*((1/perio)^2);
%         syn(istat,:) = -2.d0*a*( exp(-a*(t-to).^2) + -2.d0*a*((t-to).^2).*exp(-a*(t-to).^2) );
        AA = 12;
        syn(istat,:) = AA*2.d0*a*(t-to).*exp(-a*(t-to).^2) ;
        
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % FILTERING FOR THE ANALYTICAL SOLUTION 
    if(seismotype == 1)
        
        % Forcing signal FFT(1:NFFT/2+1)
        SYN(istat,:) = fft(syn(1,:),NFFT);
                
        % Initialization
        f = (fs/2.0)*linspace(0,fmax,NFFT/2+1);
%         f = fs*(0:(NFFT/2))/(NFFT);
%         filter(istat, :) = ones(size(SYN(istat,:)));
        sumat=zeros(1,NFFT/2+1);
        Amp1=ones(1,NFFT/2+1,1);
        Amp2=ones(1,NFFT/2+1,1);
        omega = -2*pi*f;
        
        % Setup physical param.
        c = 340;
        rho = 1;
        tau_epsilon = 1.%0.12;
        tau_sigma = 1%0.0008;
        
        % Complex modulus
%         Mc = rho * (c^2) * ( 1 + i*omega*tau_epsilon ) ./ ( 1 + i*omega*tau_sigma );
%         % Complex velocity
% %         Vc = sqrt(Mc/rho);
%         Vc = (c) * sqrt(( 1 + i*omega*tau_epsilon ) ./ ( 1 + i*omega*tau_sigma ));
%         
%         Mc2 = rho * (c^2) * ( 1 - i*omega(NFFT/2:-1:2)*tau_epsilon ) ./ ( 1 - i*omega(NFFT/2:-1:2)*tau_sigma );
% %         Vc2 = sqrt(Mc2/rho);
%         Vc2 = (c) * sqrt(( 1 - i*omega(NFFT/2:-1:2)*tau_epsilon ) ./ ( 1 - i*omega(NFFT/2:-1:2)*tau_sigma ));
        
        
%         Mc(1:NFFT/2+1) = rho * (c^2) * ( 1 + i*omega*tau_epsilon ) ./ ( 1 + i*omega*tau_sigma );
% %         Mc(1:NFFT/2+1) = conj(Mc(1:NFFT/2+1));
%         Mc(NFFT/2+2:NFFT) = rho * (c^2) * ( 1 - i*omega(NFFT/2:-1:2)*tau_epsilon ) ./ ( 1 - i*omega(NFFT/2:-1:2)*tau_sigma );
        
        Vc(1:NFFT/2+1) = (c) * sqrt(( 1 + i*omega*tau_epsilon ) ./ ( 1 + i*omega*tau_sigma ));
        Vc(1:NFFT/2+1) = conj(Vc(1:NFFT/2+1));
        Vc(NFFT/2+2:NFFT) = (c) * sqrt(( 1 - i*omega(NFFT/2:-1:2)*tau_epsilon ) ./ ( 1 - i*omega(NFFT/2:-1:2)*tau_sigma ));
%         Vc(NFFT/2+2:NFFT) = conj(Vc(NFFT/2+2:NFFT));
        
        % Green response function
%         for omega_i = 1:length(omega)
%         eta = 0.;
%         t0  = to+0;
%         omega0 = 2*pi*1/(perio);
%         eps = 1;
        
        % Carcione 1998
%         F = pi*(pi/eta)*(1/omega0)*...
%             exp(i*omega*t0).*(exp(-((pi*pi)/eta)*( eps/2 - omega/omega0 ).^2) + exp(-((pi*pi)/eta)*( eps/2 + omega/omega0 ).^2));
%         
%         filter(istat,1:NFFT/2+1) = -i*pi*besselh(0,2,R*omega./Vc).*F;

      % Ostashev
        wx = 0;%-40.5;
%         M = wx./Vc;
%         A = 1;
%         k = omega./Vc;
        alpha = 0;
%         filter(istat,1:NFFT/2+1) = A*exp(i*(sqrt(1 - M.*M*sin(alpha)^2)-M*cos(alpha)).*k.*R./(1 - M.*M) + i*pi/4).*...
%             (sqrt(1-M.*M*sin(alpha)^2) - M*cos(alpha))./(sqrt(2*pi*k*R).*(1-M.*M).*(1 - M.*M*sin(alpha)^2)).^(3/4);
%         
%         filter(istat,1:NFFT/2+1) = filter(istat,1:NFFT/2+1).*omega/i;
%         filter(istat,NFFT/2+2:NFFT) = conj(filter(istat,NFFT/2:-1:2));
%         filter(istat,1)=0.0;
        
        k(istat,1:NFFT/2+1)    = 0.0+(omega./(Vc(1:NFFT/2+1))) ;
        k(istat,1:NFFT/2+1)    = conj(k(istat,1:NFFT/2+1));
%         k(istat,NFFT/2+2:NFFT) = -omega(NFFT/2:-1:2)./(Vc(NFFT/2:-1:2));
        k(istat,NFFT/2+2:NFFT) = -omega(NFFT/2:-1:2)./(Vc(NFFT/2+2:NFFT));
        k(istat,NFFT/2+2:NFFT) = conj(k(istat,NFFT/2+2:NFFT));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Af(1:NFFT/2+1)    = 0.0+A*omega*c*c./((Vc(1:NFFT/2+1).^2)*i) ;
        Af(1:NFFT/2+1)    = conj(Af(1:NFFT/2+1));
        Af(NFFT/2+2:NFFT) = -A*omega(NFFT/2:-1:2)*c*c./((Vc(NFFT/2+2:NFFT).^2)*i);
        Af(NFFT/2+2:NFFT) = conj(Af(NFFT/2+2:NFFT));
        
%         Af(1:NFFT/2+1)    = omega./(i.*(Vc(1:NFFT/2+1).^2)) ;
%         Af(1:NFFT/2+1)    = conj(Af(1:NFFT/2+1));
%         Af(NFFT/2+2:NFFT) = -omega(NFFT/2:-1:2)./(i.*(Vc(NFFT/2+2:NFFT).^2));
        R  = xstattab(istat)-500%istat%alti(istat);%xstattab(istat) - xo
        
%         Af(1:NFFT/2+1)    = 0.0+(rho.*Vc(1:NFFT/2+1)./(2*rho*c^2)).*omega/(1*i) ;
%         Af(1:NFFT/2+1)    = conj(Af(1:NFFT/2+1));
%         Af(NFFT/2+2:NFFT) = -(rho*Vc(NFFT/2+2:NFFT)./(2*rho*c^2)).*omega(NFFT/2:-1:2)/(1*i);
%         Af(NFFT/2+2:NFFT)   = conj(Af(NFFT/2+2:NFFT));
        
        M(1:NFFT/2+1)    = wx./Vc(1:NFFT/2+1);
        M(1:NFFT/2+1)    = conj(M(1:NFFT/2+1));
        M(NFFT/2+2:NFFT) = wx./Vc(NFFT/2+2:NFFT);
        M(NFFT/2+2:NFFT) = conj(M(NFFT/2+2:NFFT));
        
        alpha = 0%acos((xstattab(istat)-xo_explo)/R);
%         filter(istat,1:NFFT/2+1)    = exp(i*(k(istat,1:NFFT/2+1)*alti(istat)))  ;
%         filter(istat,1)             = 0 ;
%         filter(istat,NFFT/2+2:NFFT) = exp(i*(k(istat,NFFT/2+2:NFFT)*alti(istat)))  ;
       
%         t0_explo = 0;

%         Computation of the frequency filter
% %         filter(istat,1:NFFT/2+1)    = -( Af(1:NFFT/2+1) ./sqrt(2*pi*k(istat,1:NFFT/2+1)*R ) ).*...
%         filter(istat,1:NFFT/2+1)    = -( Af(1:NFFT/2+1) ./(4*pi*R) ).*...    
%             ( (sqrt(1 - M(1:NFFT/2+1).*M(1:NFFT/2+1)*sin(alpha)^2) - M(1:NFFT/2+1)*cos(alpha))./( (1-M(1:NFFT/2+1).*M(1:NFFT/2+1)).*(1-M(1:NFFT/2+1).*M(1:NFFT/2+1)*sin(alpha)^2).^(3/4) ) )...
%             .*exp((i./( 1 - M(1:NFFT/2+1).*M(1:NFFT/2+1) )).*(sqrt( 1 - M(1:NFFT/2+1).*M(1:NFFT/2+1)*sin(alpha).^2 ) - M(1:NFFT/2+1)*cos(alpha)).*(k(istat,1:NFFT/2+1)*R   ))  ;
%         filter(istat,1)             = 0 ;
% %         filter(istat,NFFT/2+2:NFFT) = ( Af(NFFT/2+2:NFFT) ./sqrt(2*pi*k(istat,NFFT/2+2:NFFT)*R ) ).*...
%         filter(istat,NFFT/2+2:NFFT) = ( Af(NFFT/2+2:NFFT) ./(4*pi*R ) ).*...
%             ( (sqrt(1 - M(NFFT/2+2:NFFT).*M(NFFT/2+2:NFFT)*sin(alpha)^2) - M(NFFT/2+2:NFFT)*cos(alpha))./( (1-M(NFFT/2+2:NFFT).*M(NFFT/2+2:NFFT)).*(1-M(NFFT/2+2:NFFT).*M(NFFT/2+2:NFFT)*sin(alpha)^2).^(3/4) ) )...
%             .*exp((i./( 1 - M(NFFT/2+2:NFFT).*M(NFFT/2+2:NFFT) )).*(sqrt( 1 - M(NFFT/2+2:NFFT).*M(NFFT/2+2:NFFT)*sin(alpha).^2 ) - M(NFFT/2+2:NFFT)*cos(alpha)).*(k(istat,NFFT/2+2:NFFT)*R) )  ;
        
        filter(istat,1:NFFT/2+1)    = ( Af(1:NFFT/2+1) ./sqrt(2*pi*k(istat,1:NFFT/2+1)*R) ).*...    
            ( (sqrt(1 - M(1:NFFT/2+1).*M(1:NFFT/2+1)*sin(alpha)^2) - M(1:NFFT/2+1)*cos(alpha))./( (1-M(1:NFFT/2+1).*M(1:NFFT/2+1)).*(1-M(1:NFFT/2+1).*M(1:NFFT/2+1)*sin(alpha)^2).^(3/4) ) )...
            .*exp((i./( 1 - M(1:NFFT/2+1).*M(1:NFFT/2+1) )).*(sqrt( 1 - M(1:NFFT/2+1).*M(1:NFFT/2+1)*sin(alpha).^2 ) - M(1:NFFT/2+1)*cos(alpha)).*(k(istat,1:NFFT/2+1)*R + i*pi/4  ))  ;
        filter(istat,1)             = 0 ;
%         filter(istat,NFFT/2+2:NFFT) = ( Af(NFFT/2+2:NFFT) ./sqrt(2*pi*k(istat,NFFT/2+2:NFFT)*R ) ).*...
        filter(istat,NFFT/2+2:NFFT) = ( Af(NFFT/2+2:NFFT) ./sqrt(2*pi*k(istat,NFFT/2+2:NFFT)*R ) ).*...
            ( (sqrt(1 - M(NFFT/2+2:NFFT).*M(NFFT/2+2:NFFT)*sin(alpha)^2) - M(NFFT/2+2:NFFT)*cos(alpha))./( (1-M(NFFT/2+2:NFFT).*M(NFFT/2+2:NFFT)).*(1-M(NFFT/2+2:NFFT).*M(NFFT/2+2:NFFT)*sin(alpha)^2).^(3/4) ) )...
            .*exp((i./( 1 - M(NFFT/2+2:NFFT).*M(NFFT/2+2:NFFT) )).*(sqrt( 1 - M(NFFT/2+2:NFFT).*M(NFFT/2+2:NFFT)*sin(alpha).^2 ) - M(NFFT/2+2:NFFT)*cos(alpha)).*(k(istat,NFFT/2+2:NFFT)*R + i*pi/4) )  ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Computation of the final signal => signal * filter
        SYNf(istat,:) = SYN(istat,:) .* filter(istat,:) ;
        synf2(istat,:) = real(ifft(SYNf(istat,:),NFFT)) ;      % Solution analytique
        synf(istat,:) = synf2(istat,1:length(Ztime(nstat,:)));
      
      end % seismotype == 1
     
   end % end for istat
    
    ampf = ( max(syn(1,:)) - min(syn(1,:)) ) / ...
            ( max(real(synf(1,:))) - min(real(synf(1,:))) ) ;
        
end % if(analytic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure of the nstat vertical components and synthetic signals against time
figure
set(0,'defaulttextinterpreter','latex')
for istat = 1 : nstat

    subplot(nstat,1,nstat-istat+1) ; 
    hold on
    
    if(analytic)
        plot(Ztime(istat,:),-real(synf(istat,1:length(Ztime(istat,:)))),'Color',[0 0 0],'LineWidth',1)
    end
    
    if(modeled)
        plot(Ztime_save(istat,1:nt),-Zamp(istat,1:nt),'-.b','LineWidth',2)
    end
    
%     plot(Ztime_save(istat,1:nt),Zamp(istat,1:nt)-real(synf(istat,1:length(Ztime(istat,:)))),'-.b','LineWidth',2)
    
    if (istat == nstat)
        %title({[forcing_way];[variable,'z  for test in homogeneous medium with'];...
        title({[subtitle];[subtitle2];['']},...
            'FontSize',24)
        text(0.943,1.03,['Altitude : ',num2str(alti(istat)),'km']) %; amplitude ratio =',num2str(ampZ(istat))])
    else
        if(alti(istat)/1000 < 10)
            text(0.96*max(Ztime(istat,:)),1.03,['Altitude : ',num2str(alti(istat)/1000),'km'])   
        else if(alti(istat)/1000 < 100)
            text(0.951*max(Ztime(istat,:)),1.03,['Altitude : ',num2str(alti(istat)/1000),'km'])        
            else
            text(0.943*max(Ztime(istat,:)),1.03,['Altitude : ',num2str(alti(istat)/1000),'km'])
            end
        end
    end
    
    if (istat == round(nstat/2))
        ylabel([variable,' z-axis (m)'] ,'FontSize',24.3)
    end
    if (istat == 1)
        xlabel('time (s)','FontSize',24.3)
    end
            
%     if(analytic)
%     h_legend = legend('Analytical','Modeled','Modeled 2','Location','West') 
%     end
    
    % Limit of y axis
    if(axis_bounded)
        ylim([-1 1])
    end
        
end