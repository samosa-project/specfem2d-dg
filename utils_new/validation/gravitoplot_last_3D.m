%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                       gravitoplot4.m                             %%%%%
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

clear all; clc ; 
% close all;
format long ;
clear DIFF
clear synf

%%%%%%%%%%%%%%%%%%%%
% Kind of simulation and seismotype (Found in Par_file)
seismotype = 1 ;  
simu = 1;               % 1 for forcing and 2 for attenuation

%%%%%%%%%%%%%%%%%%%%%
% Stations parameters
% nstat      = 7 ;             % Number of station.
% dx_station = 60000; %58333.3333333;   % x-distance between stations
% z_station  = 100250.0;        % Height along z of stations
% % x_0        = 150000.0; % first station along x
% x_0        = 600250.0; % first station along x

% istattab=[6:7:48]
% xstattab=[450250.0:50000.0:750250.0]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time and space domain variables
% xmin = 350000;
% xmax = 1250000.0;
xmin = 400000;
xmax = 600000.0;
ymin = xmin;
ymax = xmax;

dx   = 5000.0;
dy   = dx;
t0 = 0.0;
% dt   = 25.0 ;
% syn_duration=400000 ;
dt   = 100.0 ;
extent_time = 3.5;
extent_x = 1;
extent_y = 1;

% Frequency domain
fs   = 1/dt ;
fmax = 1 ;
fmin = -fmax ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical model variables
%  xo     = 750000.0;
% xo     = 250000.0;
% lambdo = 10000.0;
lambdo = 80000.0;
V      =  639.52000;       % must be in m/s
A      = 1.0;
% rho_cst = 0.02 ;
% Time function parameters
perio  = 1600.0;
to     = 1400.0;
% perio  = 50.0;
% to     = 65.0;
% xo = 600000.0;
% wind in x direction
% uo=99.999
% %  uo=4.9
uo = 10
% uo = 0.0

%%%%%%%%%
% Display
subtitle = strcat('Forcing with the whole surface moving at the bottom of a fluid medium');
subtitle2 = strcat([' with a variable density profile and a velocity equal to ',num2str(V),' m.s^{-1}']) ;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory of semd files
% directory = strcat('/home/garcia/SATGRAVI/Brissaud/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4/') ;
% directory = strcat('/home/garcia/SATGRAVI/Brissaud/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_LDDRK/') ;
% directory2 = strcat('/home/garcia/SATGRAVI/Brissaud/test_solution_analytique_gravi/Results.semd/') ;
% directory = strcat('/home/garcia/SATGRAVI/Brissaud/files_wind_quentin/gravi_forcing_nowind/') ; 
%directory = strcat('/home/garcia/SATGRAVI/Brissaud/files_wind_quentin/forcing_gravi_1Diso_windcte_noatten/') ; 
%  directory = strcat('/home/garcia/SATGRAVI/Brissaud/files_wind_quentin/FG1_1Diso_wind_noatten/') ; 
% directory = strcat('/home/garcia/SATGRAVI/Brissaud/files_wind_quentin/FG1_1Diso_nowind_noatten/') ; 
% directory = strcat('/home/qbrissaud/Documents/Results/GJI_PAPER/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4_3D/');
% directory = strcat('/home/qbrissaud/Documents/Results/LAST_GRAVI_Roland/200PROCS/');
directory = strcat('./');
% directory = strcat('/home/qbrissaud/Documents/FD/FD_14/');
% directory = strcat('/home/qbrissaud/Documents/Results/GJI_PAPER/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4_3D/');
% directory = strcat('/home/qbrissaud/Documents/Results/GJI_PAPER/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4/') ; 
% directory = strcat('/home/qbrissaud/Documents/Results/GJI_PAPER/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beginning of the program %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading data from the atmospherical model
atmos = load(strcat(directory,'1D_iso_enhanced.txt')) ;
% atmos = load('1D_iso_enhanced.txt') ;

%%%%%%%%%%%%%%%%%%%%%%%%%
% File parameters reading
alti_atm = atmos(:,1) ;
% % rho      = 1*size(atmos(:,2)) ;
% velocity = 342*ones(size(atmos(:,3),1)) ;
% gravity  = -9.81*ones(size(atmos(:,4),1)) ;
% % Nsq      = atmos(:,5) ;
% Nsq      = (-1/(342*342))*gravity(1)*ones(size(atmos(:,5),1));
% Pressure = atmos(:,6) ;
% % Profile of density with a constant scaling H
% H        = (342*342)/(9.81+1)%8.6611888899e+03;
% rho_0    = 1.2;
% rho     = rho_0 * exp(-alti_atm/H) ;

rho      = atmos(:,2) ;
velocity = atmos(:,3) ;
gravity  = atmos(:,4) ;
Nsq      = atmos(:,5) ;
Pressure = atmos(:,6) ;
% Profile of density with a constant scaling H
% H        = (alti_atm(2)-alti_atm(1))/log(rho(1)/rho(2))
H        = (-alti_atm(2))/log(rho(2)/rho(1))
% Nsqt      = -(gravity(1)*( -1/H - gravity(1)/(velocity(1)^2) ));
rho_0    = rho(1);

nsyn=10000;
xstattab = [5000];
    zstattab = [10000];
%      xo = xmax/2;
    nstat = 1;
    Ztime(1,:) = linspace(0,nsyn,dt);
READ_DATA = 0;
    
if(READ_DATA == 1)
%%%%%%%%%%%%%%%%%%%%%%
% Stations data
pos_stat = load(strcat(directory,'stations.txt')) ;

% istattab = 1:size(pos_stat(:,1),1);
% istattab = [1:18:81] ;
% istattab = [37 38 39 40 41 42 43 44 45];
% istattab = [36 37 44 45]+9;
coef = 9;
alti = 4;
istattab = [alti+coef alti+2*coef alti+7*coef alti+8*coef]

% istattab = [1:5] + 1*9;
istattab = [2 4 6 8] + 2*9;
% istattab = 5;
% istattab = [1 2 8 9] + 4*9
% istattab = [1 10 19 28] + 8
% istattab = [10 19 55 64] + 0
% istattab = [19 28 46 55] + 0
% istattab = [13 22 49 58]

% istattab = [1 19 55 73]
% istattab = [1 10 64 73]+1
%istattab = 10:18:81;

%%%%%%%%%%%
% TEST GRAVI NO WIND
% ystattab = pos_stat(istattab,1)+250;
% zstattab = pos_stat(istattab,2)+250 ;
% zstattab(1) = zstattab(1) + 20000;
% zstattab(2) = zstattab(2) + 5000;
% zstattab(3) = zstattab(3) + 20000;
% zstattab(4) = zstattab(4) + 20000;
% xstattab = pos_stat(istattab,3)+250 ;
% xstattab(3:4) = xstattab(3:4) + 0;
% xstattab(1) = xstattab(1) - 5000;
% xstattab(2) = xstattab(2) - 5000;

%%%%%%%%%%%
% TEST GRAVI  WIND
% ystattab = pos_stat(istattab,1);
% zstattab = pos_stat(istattab,2) ;
% xstattab = pos_stat(istattab,3) ;
% 
% zstattab(1) = zstattab(1) + 7200;
% xstattab(1) = xstattab(1) - 0;
% 
% zstattab(2) = zstattab(2) + 5000;
% xstattab(2) = xstattab(2) + 3000;
% 
% zstattab(3) = zstattab(3) + 11000;
% xstattab(3) = xstattab(3) - 2500;
% 
% zstattab(4) = zstattab(4) + 13500;
% xstattab(4) = xstattab(4) + 1500;

%%%%%%%%%%%
% TEST GRAVI  WIND
ystattab = pos_stat(istattab,1);
zstattab = pos_stat(istattab,2)+0 ;
xstattab = pos_stat(istattab,3) ;

% zstattab(1) = zstattab(1) + 10000;
xstattab = xstattab + 000;
% xstattab(4) = xstattab(4) + 0;
ystattab = ystattab + 000;

xo = (pos_stat(end,1) + pos_stat(1,1))/2;
xo = xo + 0;
yo = xo;
% yo = (pos_stat(end,3) + pos_stat(1,3))/2;

nstat    = length(istattab) ;

N             = sqrt(Nsq(1))*10;
beta          = acos( 2*pi*(1/perio)/N );
beta_stat_tot = atan(pos_stat(:,2)./(pos_stat(:,1)-xo));
beta_stat     = beta_stat_tot(istattab);
pos_ok        = istattab.*(beta_stat<beta)';

%%%%%%%%%%%%%%%%%%%%%%
% title for simulation
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First loop over the stations to read semd files
% subsample the results to save memory space
nsub=100
for jstat = 1 : nstat
    istat=istattab(jstat)
    if(istat < 10)
       prefix = 'S000';
    elseif(istat < 100)
       prefix = 'S00'; 
    else
       prefix = 'S0';
    end
    % reading of the horizontal component
    file = strcat(directory,prefix,num2str(istat),'.AA.BXX.sem',signal_type) ;
    data = load(file) ;
    %if (istat == 1)
    nt = floor(max(size(data))/nsub) ;
    %end
    Xtime(jstat,1:nt) = data(1:nsub:max(size(data)),1)' ;
    Xamp(jstat,1:nt) = data(1:nsub:max(size(data)),2)' ;
    
    % reading of the vertical component
    file = strcat(directory,prefix,num2str(istat),'.AA.BXZ.sem',signal_type) ;
    data = load(file) ;
    nt=floor(max(size(data))/nsub);
    Ztime(jstat,1:nt) = data(1:nsub:max(size(data)),1)' ;
    Zamp(jstat,1:nt) = data(1:nsub:max(size(data)),2)' ;
end
clear data
dt = dt;
nsyn=nt;

end

% the square root of -1
i = complex(0,1);

length(Ztime(nstat,:))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of the analytical solution
% syn = zeros(nstat,length(Ztime(nstat,:))) ;
% synX = zeros(nstat,length(Ztime(nstat,:))) ;
% synZ = zeros(nstat,length(Ztime(nstat,:))) ;
% NFFT2 = 2^nextpow2(length(Ztime(nstat,:))); 
syn = zeros(nstat,nsyn) ;
synX = zeros(nstat,nsyn) ;
synZ = zeros(nstat,nsyn) ;
NFFT2 = 2^nextpow2(nsyn)*extent_time
% dt=256*dt; % ?????
% Fourier space
NFFT1 = 2^nextpow2((xmax-xmin)/dx)*extent_x
NFFT3 = 2^nextpow2((ymax-ymin)/dy)*extent_y
display('FFT3D matrix size :')
display(NFFT1*NFFT2*NFFT3)
k = zeros(NFFT1,NFFT2,NFFT3);
t = zeros(1,NFFT2) ;
t = dt * (0:1:NFFT2-1);
maxwind = max(t);
x = dx * (0:1:NFFT1-1) + xmin;
y = dy * (0:1:NFFT3-1) + ymin;
omega = 2.0*pi()*(1.0/(dt*NFFT2))*[[0:1:NFFT2/2] [0.0-[NFFT2/2-1:-1:1]]];
kx = 2.0*pi()*(1.0/(dx*NFFT1))*[[0:1:NFFT1/2] [0.0-[NFFT1/2-1:-1:1]]];
ky = 2.0*pi()*(1.0/(dy*NFFT3))*[[0:1:NFFT3/2] [0.0-[NFFT3/2-1:-1:1]]];

% IN THIS SENSE? -> NO!
%omega = 2.0*pi()*(1.0/(dt*NFFT2))*[ [0.0-[NFFT2/2-1:-1:1]] [0:1:NFFT2/2]];
%kx = 2.0*pi()*(1.0/(dx*NFFT1))*[[0.0-[NFFT1/2-1:-1:1]] [0:1:NFFT1/2] ];

% Define the starting Matrix
%     tendsig=1500.0;
%     [X,Y,T]=meshgrid(x,y,t);
    [X,T]=meshgrid(x,t);
%     Xp = X+2500;
%     Yp = Y+2500;
    Mo=  2 * ( exp(-((X-(xo-lambdo/4))/(lambdo/4)).^2) ...
             - exp(-((X-(xo+lambdo/4))/(lambdo/4)).^2) ).*...
             ( exp(-((T-(to-perio/4)-t0)/(perio/4)).^2)...
             - exp(-((T-(to+perio/4)-t0)/(perio/4)).^2) );%...
%              .* exp(-(sqrt((X-xo).^2+(Y-yo).^2)/(lambdo/4)).^2);
    TFMo = fftn(Mo);
%%%%%%    IMPORTANT CORRECTION   %%%%
% fft2 of matlab perform the projection on the function basis of the type:
% exp[i(Kx*X+Ky*Y)] whereas we use exp[i(Kx*X-omega*t)] !!!
% so:
% * Ky = -omega
% * and positive frequencies are in the second part of the fft table
% this is corrected by inserting minus sign in the expression below and by
% using positive frequencies of the second part of the fft table to ensure
% symmetry of the fft (see below)
% the formulas themselves (for dispersion relation) are not modified
%%%%%%    IMPORTANT CORRECTION   %%%%
    [KX,KY,Omega]=meshgrid(kx,ky,-omega);
% Assuming constant Nsq
    Nsqtab=Nsq(1,1)+0.0*Omega;
    
%     KX = Omega/velocity(1);
    
% see occhipinti 2008 for analytical solution (appendix A1)
    onestab=0.0*Nsqtab+1.0;
    
    omega_intr = Omega-(uo)*KX;
    
    KZ=sqrt(Nsqtab.*(KX.*KX)./((Omega-(uo)*KX).*(Omega-(uo)*KX))-(KX.*KX) + (1 - Nsqtab./(Omega.*Omega) ).*((Omega-(uo)*KX).*(Omega-(uo)*KX))/(velocity(1)^2) );
%     KZ=sqrt( Nsqtab.*(KX.*KX + KY.*KY)./(omega_intr.*omega_intr)-(KX.*KX + KY.*KY) - onestab/(4*H*H) + (omega_intr.*omega_intr)/(velocity(1)^2) );

    ind1=find(isnan(KZ));
    size(ind1)
    ind2=find(isinf(KZ));
    size(ind2)
    KZ(ind1)=0.0;
    KZ(ind2)=0.0;
% imaginary part of KZ should be positive in order to attenuate the signal
    indimag=find(imag(KZ)<0);
    KZ(indimag)=conj(KZ(indimag));

% real(KZ) should be positive for positive frequencies and negative for
% negative frequencies in order to shift signal in positive times
% restore the sign of KZ depending on Omega-uo*KX
%     KZnew=real(KZ).*sign((Omega-uo*KX)).*sign(KX)+1i*imag(KZ);
% !!! Why KZ should have a sign opposite to Omega for GW NOT UNDERSTOOD !!!
% => because vg perpendicular to Vphi ?
      KZnew=0.0-real(KZ).*sign(omega_intr)+1i*imag(KZ);
      KZ=KZnew;

% restore negative frequencies from postive ones:
%      KZ(:,NFFT2/2+2:NFFT2)=0.0-real(KZ(:,NFFT2/2:-1:2))+1i*imag(KZ(:,NFFT2/2:-1:2));
% restore negative frequencies from postive ones:
% Corrected for "%%%%%%    IMPORTANT CORRECTION   %%%%"
%      KZ(:,NFFT2/2:-1:2)=0.0-real(KZ(:,NFFT2/2+2:NFFT2))+1i*imag(KZ(:,NFFT2/2+2:NFFT2));
% Above not taken into account for "filt", because only positive computations 

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTERING FOR THE ANALYTICAL SOLUTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(seismotype == 1)
        
    if(simu == 1)   % forcing
        
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
% but frequencies close to (uo*KX-Omega)=0 create very large KZ...
%         TFMx=TFMo.*(KZ./KX).*filt;
%         ind1=find(isnan(TFMx));
%         ind2=find(isinf(TFMx));
%         TFMx(ind1)=0.0+i*0.0;
%         TFMX(ind2)=0.0+i*0.0;
%         Mx=0.0-exp(z_station/(2*H))*ifft2(TFMx);
%        xstattab_analytic = size(xstattab);
%        ystattab_analytic = size(ystattab);
        for istat=1:nstat
            
            z_station = zstattab(istat);
           filt=exp(i*(KZ*z_station));
           filt(1,1,:)=0.0;
            Mz=1*exp(z_station/(2*H))*ifftn((filt.*TFMo));
            
%           xcoord(istat) = xstattab(istat);
          ix=round((xstattab(istat)-xmin)/dx) + 1;
%           iy=round((ystattab(istat)-ymin)/dy) + 1;
          
          % GET X LOCATION OF STATION
          ACHECK = xstattab(istat);
          POS_sat_retrieve = X(1,ix,1);
          
           POS_sat_retrieve_p1 = X(1,ix+1,1);
           POS_sat_retrieve_m1 = X(1,ix,1);
           
           if(abs(POS_sat_retrieve - ACHECK) > abs(POS_sat_retrieve_p1 - ACHECK))
               ix = ix + 1;
           elseif(abs(POS_sat_retrieve - ACHECK) > abs(POS_sat_retrieve_m1 - ACHECK))
               ix = ix - 1;
           end
           
%            
           
           xstattab_analytic(istat) = X(1,ix,1)
               
           % GET Y LOCATION OF STATION
          ACHECK = ystattab(istat);
          POS_sat_retrieve = Y(iy,1,1);
          
           POS_sat_retrieve_p1 = Y(iy+1,1,1);
           POS_sat_retrieve_m1 = Y(iy,1,1);
           
           if(abs(POS_sat_retrieve - ACHECK) > abs(POS_sat_retrieve_p1 - ACHECK))
               iy = iy + 1;
           elseif(abs(POS_sat_retrieve - ACHECK) > abs(POS_sat_retrieve_m1 - ACHECK))
               iy = iy - 1;
           end
           
            ystattab_analytic(istat) = Y(iy,1,1)
           
           % IF NOT PART OF THE MESH
          if (ix<=0)
              ix=1
          end
          
          if (iy<=0)
              iy=1
          end
%           if(istat < 3)
%               ix = ix-1;
%           else
%               ix = ix+1;
%           end
%           synf(istat,:)=real(Mz(iy,ix,:));
          synf(istat,:)=real(Mz(ix,:));
   % remove offset a the beginning may be due to aliasing
%           synf(istat,:)=synf(istat,:)-synf(istat,1);
%           synf_X(istat,:)=synf_X(istat,:)-synf_X(istat,1);
        end
    end
end

% test only !!!!
% synf_X=synf;
%

% clear TFMx;
% clear Mx ;
% clear Mz ;

timef = t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure of the nstat vertical components and synthetic signals against 
% time
figure
for istat = 1 : nstat
    
    ax(istat)=subplot(nstat,1,nstat-istat+1) ; 
    hold on
    
    if(seismotype ==1)
         timef_final = find(timef<=Ztime(istat,end));
         timef_final = timef_final(end);
%          plot(timef(1:timef_final),real(synf(istat,1:timef_final)),'Color',[0 0 0],'LineWidth',1)
%        plot(timef,real(synf(istat,:))/max(real(synf(istat,:))),'Color',[0 0 0],'LineWidth',1)

        %        plot(Ztime(istat,:),real(synf(istat,1:length(Ztime(istat,:))))/max(real(synf(istat,1:length(Ztime(istat,:))))),'Color',[0 0 0],'LineWidth',1)
    else
    plot(Ztime(istat,:),syn(istat,:),'Color',[0 0 0],'LineWidth',1)
    end
  
    Ztime_temp = Ztime(istat,1:dt/(Ztime(istat,2)-Ztime(istat,1)):end);
      Zamp_temp = Zamp(istat,1:dt/(Ztime(istat,2)-Ztime(istat,1)):end);
      SIZE_R = size(Ztime_temp)
%       size((Zamp_temp-real(synf(istat,1:length(Ztime_temp))))')
      DIFF(istat,:) = (Zamp_temp-real(synf(istat,1:length(Ztime_temp))))';

    plot(Ztime(istat,1:nt),Zamp(istat,1:nt),'-.k','LineWidth',2)
    
    if (istat == round(nstat/2))
        ylabel([variable,' z-axis (m)'] ,'FontSize',14.3)
    end
    if (istat == 1)
        xlabel('time (s)','FontSize',14.3)
    end
    
    text(00.943*max(Ztime(istat,:)),0.95,['X position : $',num2str(xstattab(istat)),'km$ (along x)']) 
    
    legend('Analytical','Modeled','10*(Mo-An)','Location','East')

end

linkaxes(ax,'x')
% xlim([0 20000])

title(strcat('Gravito-acoustic wave propagation. Stations at altitude : ', num2str(zstattab(istat)),'km (along z)'),'FontSize',24)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure of the nstat horizontal components and synthetic signals against 
% time
% figure
% for istat = 1 : nstat
%     
%     ax(istat)=subplot(nstat,1,nstat-istat+1) ; 
%     hold on
%     
%     if(seismotype ==1)
%          plot(timef,real(synf_X(istat,:)),'Color',[0 0 0],'LineWidth',1)
% %        plot(timef,real(synf_X(istat,:))/max(real(synf_X(istat,:))),'Color',[0 0 0],'LineWidth',1)
% 
% %                plot(Ztime(istat,:),real(synf_X(istat,1:length(Xtime(istat,:))))/max(real(synf_X(istat,1:length(Xtime(istat,:))))),'Color',[0 0 0],'LineWidth',1)
%     else
%     plot(Xtime(istat,:),syn(istat,:),'Color',[0 0 0],'LineWidth',1)
%     end
%   
%      plot(Xtime(istat,1:nt),Xamp(istat,1:nt),'-.k','LineWidth',2)
% %    plot(Xtime(istat,:),Xamp(istat,:)/max(Xamp(istat,:)),'-.k','LineWidth',2)
%     
%     if (istat == round(nstat/2))
%         ylabel([variable,' x-axis (m)'] ,'FontSize',14.3)
%     end
%     if (istat == 1)
%         xlabel('time (s)','FontSize',14.3)
%     end
%     
%     text(00.943*max(Xtime(istat,:)),0.95,['X position : $',num2str((istat-1)*dx_station + x_0),'km$ (along x)']) 
%     
%     legend('Analytical','Modeled','10*(Mo-An)','Location','West')
% 
% end
% 
% linkaxes(ax,'x')
% % xlim([0 20000])
% 
% title(strcat('Gravito-acoustic wave propagation. Stations at altitude : ', num2str(z_station/1000),'km (along z)'),'FontSize',24)



