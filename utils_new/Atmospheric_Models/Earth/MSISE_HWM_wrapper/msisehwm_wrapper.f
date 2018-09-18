C Authors:       X. Bonnin, R. F. Garcia, L. Martire.
C Description:   TODO.
C Last modified: See file metadata.
C Usage:         N/A.
C Notes:         N/A.

C Compilation line:
C gfortran -o msisehwm msisehwm_wrapper.f msisehwm_wrapper_sub.for

C ****************************************************************
C * compute_msise_winds_v3.                                      *
C ****************************************************************
      program compute_msise_winds_v3

      implicit none

      integer imax
      parameter (imax=5000)

      real*4 latmod, lonmod

      real*4 z(imax),rhoat(imax),v(imax),T(imax),P(imax)
      real*4 Wind(3,imax),Windshift(2)
      real*4 MUtab(imax),MUvoltab(imax),Kappatab(imax)
      real*4 Gravity(imax),Nsqtab(imax),Htab(imax)
      real*4 Cpcoefs(5,9,3),tminrange(9,3),tmaxrange(9,3)
      real*4 molarmass(9),Ti,densmol
      integer nbspec,ntrange(9),indspec(9)
      real*4 Cp(imax),Cv(imax),gammatab(imax)
      real*4 xdum1,xdum2,PI,dt,zsat,tt(imax),xdum3,xdum4
      real*4 ttz
      real*4 MU,MUvolclass,Kappa
      real grav,Nsq,Nac,Fbr,Fac,Nsq1,Nsq2
      real GMass,radius,Hscale

      complex Ampz
      complex Amp(imax),Amp2(imax)

      integer i,j,izs,zsmin,zsmax,iz,k,k1,k2
      
      integer ispec
      integer indtemp
      
      integer NSAMPLES
      real*4 ALTMIN, ALTMAX
      
      integer itop
      parameter (itop=28033)

      character*14 evename
      character*37 directory
      character*33 folder
      character*70 location
C     The sizes of folder and evename are depending on the event name,
C     that's why they have to be changed for each event.
      
      real*4 time(itop),D(itop),lati(itop),lon(itop),R(itop)
      real*4 WS(2,itop)
      real*4 zint(itop),tmp,h(itop),t1,t2,f
      real*4 SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP(7),Dms(9),Tms(2)
      real*4 ap3(10),A(9),Vms,Pms,W(2)
      
      integer IYD,Z1,Z2,Zmin,Zmax,dh1,dh2,Iquake,day,year
      
      integer numargs, iarg
      character(len=500), dimension(:), allocatable :: args

      character*500 filename,filewind
      real*4 wind_projection

C     ****************************************************************
C     * Parameters.                                                  *
C     ****************************************************************
      AP = 0
      numargs = command_argument_count()
      if(numargs>0) then
        if(.not. numargs==13) then
          write(*,*) '13 and only 13 arguments should be given ',
     &         '(ALTMIN ALTMAX NSAMPLES latmod lonmod year day SEC ',
     1         'F107A F107 AP filename wind_projection).'
          stop
        endif
        allocate(args(numargs))
        do iarg = 1, numargs
          call get_command_argument(iarg,args(iarg))
        end do
        read(args(1), *) ALTMIN
        read(args(2), *) ALTMAX
        read(args(3), *) NSAMPLES
        read(args(4), *) latmod
        read(args(5), *) lonmod
        read(args(6), *) year
        read(args(7), *) day
        read(args(8), *) SEC
        read(args(9), *) F107A
        read(args(10), *) F107
        read(args(11), *) AP(1)
        read(args(12), *) filename
        read(args(13), *) wind_projection
        NSAMPLES = NSAMPLES - 1
        write(*, *) 'Arguments assigned.'
        if(NSAMPLES>imax) then
          write(*,*) "NSAMPLES > imax (", imax, "): Reduce NSAMPLES or",
     &               " hard-code a higher imax limit."
          stop
        endif
      else
C       Minimum and maximum wanted altitudes (in [m]).
        ALTMIN=0.
        ALTMAX=500000.
C       Number of samples.
        NSAMPLES=500
C       Latitude.
        latmod=36.529998
C       Longitude.
        lonmod=158.690000
C       Years since 2000 (because of msise00)
        year=11
C       Days since the beginning of the year 2000.
        day=69
C       Seconds since the beginning of day (in UTC).
        SEC=28059.996
C       81-day average of F10.7 radio flux (centered on the day of the event).
        F107A=106.7160
C       Daily F10.7 flux (for previous day).
        F107=131.
C       Daily magnetic index.
        AP=37
C       Output file name.
        filename='msisehwm_model_output'
C       Wind projection angle (in degrees). 0 is forward zonal (positive eastward). 90 is forward meridional (positive northward). 180 is backwards zonal (positive westward). 270 is backward meriodional (positive southward).
        wind_projection=0.
        write(*, *) 'No arguments passed, those in the source code',
     &              ' will be used. Output file is "', trim(filename),
     &              '".'
      endif
C     ****************************************************************
C     * Constants.                                                   *
C     ****************************************************************
C     Pi.
      PI=ACOS(-1.)
C     Earth and Atmosphere parameters.
      GMass=(6.67e-11)*(5.97219e+24)
      radius=6371010.0
C     Flattening of Earth deduced from the WGS84 value of 1/f. Used to compute the geodetic latitude (GLAT).
      f=1./298.257223563
C     81-day average of F10.7 flux (centered on the day of the event).
C     F107A=106.7160
C     Daily F10.7 flux for previous day.
C     F107=131.
C     Magnetic index (daily, for the last 24 hours).
C   - or when SW(9)=-1. it is an array containing :
C     (1) Daily AP
C     (2) 3 HR AP index for current time
C     (3) 3 HR AP index for 3 HRS before current time
C     (4) 3 HR AP index for 6 HRS before current time
C     (5) 3 HR AP index for 9 HRS before current time
C     (6) Average of eight 3 HR AP indexes from 12 to 33 HRS prior 
C                 to current time
C     (7) Average of eight 3 HR AP indexes from 36 to 57 HRS prior
C                 to current time
C     AP=37
C     To change the values of F107, F107A and Ap, watch the website:
C     http://www.swpc.noaa.gov/alerts/solar_indices.html
      
C     ****************************************************************
C     * Thermodynamic coefficients for species referenced in MSISE   *
C     ****************************************************************
C     He (http://webbook.nist.gov/cgi/cbook.cgi?ID=C7440597&Units=SI&Mask=1#Thermo-Gas)
      ispec=1
      molarmass(ispec)=4.0
      ntrange(ispec)=1
      tminrange(ispec,1)=290
      tmaxrange(ispec,1)=6000
      Cpcoefs(1,ispec,1)=20.78603
      Cpcoefs(2,ispec,1)=4.850638e-10
      Cpcoefs(3,ispec,1)=-1.582916e-10
      Cpcoefs(4,ispec,1)=1.525102e-10
      Cpcoefs(5,ispec,1)=3.196347e-11

C     N2 (http://webbook.nist.gov/cgi/cbook.cgi?Name=N2&Units=SI&cTG=on#Thermo-Gas)
      ispec=3
      molarmass(ispec)=28.0
      ntrange(ispec)=3
      tminrange(ispec,1)=100
      tmaxrange(ispec,1)=500
      tminrange(ispec,2)=500
      tmaxrange(ispec,2)=2000
      tminrange(ispec,3)=2000
      tmaxrange(ispec,3)=6000
      Cpcoefs(1,ispec,1)=28.98641
      Cpcoefs(2,ispec,1)=1.853978
      Cpcoefs(3,ispec,1)=-9.647459
      Cpcoefs(4,ispec,1)=16.63537
      Cpcoefs(5,ispec,1)=0.000117
      Cpcoefs(1,ispec,2)=19.50583
      Cpcoefs(2,ispec,2)=19.88705
      Cpcoefs(3,ispec,2)=-8.598535
      Cpcoefs(4,ispec,2)=1.369784
      Cpcoefs(5,ispec,2)=0.527601
      Cpcoefs(1,ispec,3)=35.51872
      Cpcoefs(2,ispec,3)=1.128728
      Cpcoefs(3,ispec,3)=-0.196103
      Cpcoefs(4,ispec,3)=0.014662
      Cpcoefs(5,ispec,3)=-4.553760

C     O2 (http://webbook.nist.gov/cgi/cbook.cgi?Name=O2&Units=SI&cTG=on#Thermo-Gas)
      ispec=4
      molarmass(ispec)=32.0
      ntrange(ispec)=3
      tminrange(ispec,1)=100
      tmaxrange(ispec,1)=700
      tminrange(ispec,2)=700
      tmaxrange(ispec,2)=2000
      tminrange(ispec,3)=2000
      tmaxrange(ispec,3)=6000
      Cpcoefs(1,ispec,1)=31.32234
      Cpcoefs(2,ispec,1)=-20.23531
      Cpcoefs(3,ispec,1)=57.86644
      Cpcoefs(4,ispec,1)=-36.50624
      Cpcoefs(5,ispec,1)=-0.007374
      Cpcoefs(1,ispec,2)=30.03235
      Cpcoefs(2,ispec,2)=8.772972
      Cpcoefs(3,ispec,2)=-3.988133
      Cpcoefs(4,ispec,2)=0.788313
      Cpcoefs(5,ispec,2)=-0.741599
      Cpcoefs(1,ispec,3)=20.91111
      Cpcoefs(2,ispec,3)=10.72071
      Cpcoefs(3,ispec,3)=-2.020498
      Cpcoefs(4,ispec,3)=0.146449
      Cpcoefs(5,ispec,3)=9.245722

C     Ar (http://webbook.nist.gov/cgi/cbook.cgi?Name=Ar&Units=SI&cTG=on#Thermo-Gas)
      ispec=5
      molarmass(ispec)=40.0
      ntrange(ispec)=1
      tminrange(ispec,1)=290
      tmaxrange(ispec,1)=6000
      Cpcoefs(1,ispec,1)=20.78600
      Cpcoefs(2,ispec,1)=2.825911e-7
      Cpcoefs(3,ispec,1)=-1.464191e-7
      Cpcoefs(4,ispec,1)=1.092131e-8
      Cpcoefs(5,ispec,1)=-3.661371e-8

C     N (http://webbook.nist.gov/cgi/cbook.cgi?ID=C17778880&Units=SI&Mask=1#Thermo-Gas)
      ispec=8
      molarmass(ispec)=14.0
      ntrange(ispec)=1
      tminrange(ispec,1)=290
      tmaxrange(ispec,1)=6000
      Cpcoefs(1,ispec,1)=21.13581
      Cpcoefs(2,ispec,1)=-0.388842
      Cpcoefs(3,ispec,1)=0.043545
      Cpcoefs(4,ispec,1)=0.024685
      Cpcoefs(5,ispec,1)=-0.025678

C     Reordering.
      nbspec=5
      indspec(1)=1
      indspec(2)=3
      indspec(3)=4
      indspec(4)=5
      indspec(5)=8

C     ****************************************************************
C     * Creation of the atmosphere and wind model.                   *
C     ****************************************************************
C     IYD like this because format needed in calls to GTD7 and to GWS5 (YYDDD).
      IYD=year*1000.0+day
      GLAT=atan(tan(latmod*PI/180.0)/((1.-f)*(1.-f)))*180./PI
      if(lonmod.lt.0.0) then
        GLONG=lonmod+360.0
      else
        GLONG=lonmod
      endif
C     Solar local time.
      STL=SEC/3600.+GLONG/15.
      if(STL.gt.24.0) then
        STL=STL-24.0
      endif

C     Activate SI units output for GTD7.
C     CALL METERS(.TRUE.)

C     Open output file.
      open(2012, file=filename, form='formatted')
C     Write information about model.
      write(2012,*) 'year', (2000+year), 'day', day, 'seconds', SEC,
     &              'lat', latmod, 'lon', lonmod
      write(2012,*) 'IYD', IYD, 'SEC', SEC, 'GLAT', GLAT, 'GLON',
     &              GLONG, 'STL', STL, 'F107A', F107A, 'F107', F107,
     &              'AP', AP
C     First loop on altitude samples to get models.
      do iz=1, NSAMPLES+1
C       Altitude of calculation for model calls. [km].
        ALT = real(ALTMIN + (iz-1)*(ALTMAX-ALTMIN)/NSAMPLES)/1000.
C       Call models.
        call GTD7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,48,Dms,Tms)
        call GWS5(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,W)
C       Altitude of calculation for model saving. [m].
        z(iz)=ALT*1000.0
C       Global atmospheric density.
        rhoat(iz)=Dms(6)*1000.0
C       rhoat(iz)=Dms(6)
C       Atmospheric temperature. [K].
        T(iz)=Tms(2)

C       Beginning of the computation of c_p ([J/mol*K]).
        Ti=T(iz)/1000.0
C       Compute density of the gas mixture.
        densmol=0.0
        do i=1,nbspec
	        ispec=indspec(i)
	        densmol=densmol+Dms(ispec)
        enddo
C       Add O (atomic oxygen) and H (atomic hydrogen) densities to the molar density.
        densmol=densmol+Dms(2)+Dms(7)
C       Compute isobaric heat capacity of the gas mixture (in [J/mol*K]).
        Cp(iz)=0.0
        do i=1,nbspec
	        ispec=indspec(i)
          indtemp=1
	        do k=1,ntrange(ispec)
	          if((T(iz).ge.tminrange(ispec,k)) .and.
     &         (T(iz).le.tmaxrange(ispec,k))) then
	            indtemp=k
	          endif
	        enddo
	Cp(iz)=Cp(iz)+(Dms(ispec)/densmol)*(Cpcoefs(1,ispec,indtemp) +
     &   Cpcoefs(2,ispec,indtemp)*Ti +
     &   Cpcoefs(3,ispec,indtemp)*Ti*Ti + 
     &   Cpcoefs(4,ispec,indtemp)*Ti*Ti*Ti +
     &   Cpcoefs(5,ispec,indtemp)/(Ti*Ti))
        enddo
C       Add O (atomic oxygen) and H (atomic hydrogen).
        Cp(iz)=Cp(iz)+(Dms(2)/densmol)*(5./2.)*8.3145
        Cp(iz)=Cp(iz)+(Dms(7)/densmol)*(5./2.)*8.3145
C       Isochoric heat capacity. [J/mol*K].
        Cv(iz)=Cp(iz)-8.3145
C       Heat capacity ratio.
        gammatab(iz)=Cp(iz)/Cv(iz)
C       Atmospheric pressure. [Pa]=[kg.s^(-2).m^(-1)].
C       P(iz)=rhoat(iz)*8.3145*T(iz)/0.0289645
        P(iz)=(densmol*1000./6.02214129e23)*(8.3145*1000.)*T(iz)
C       P(iz)=(densmol/6.02214129e23)*(8.3145)*T(iz)
C       Sound velocity. [m.s^(-1)].
        v(iz)=sqrt(gammatab(iz)*P(iz)/rhoat(iz))
C       Meridional (1) and zonal (2) winds. [m.s^(-1)].
        Wind(1,iz)=W(1)
        Wind(2,iz)=W(2)	
C       Projected wind.
        Wind(3,iz)=cos(wind_projection*PI/180.)*W(2)+
     &             sin(wind_projection*PI/180.)*W(1)
C       Thermal conductivity. [kg.m.s^(-3).K^(-1)].
        Kappatab(iz)=Kappa(T(iz))
C       Dynamic viscosity. [Pa.s]=[kg.s^(-1).m^(-1)].
        MUtab(iz)=MU(rhoat(iz),P(iz),T(iz),gammatab(iz))
C       Volumic viscosity. [Pa.s]=[kg.s^(-1).m^(-1)].
        MUvoltab(iz)=MUvolclass(MUtab(iz),Kappatab(iz),gammatab(iz),
     &                          rhoat(iz),T(iz),P(iz))
      enddo
    
C     Print column headers to file.
      write(2012,*) 'z[m] rho[kg/(m^3)] T[K] c[m/s] p[Pa] H[m] ',
     &              'g[m/(s^2)] N^2[rad^2/s^2] kappa[J/(s.m.K)] ',
     &              'mu[kg(s.m)] mu_vol[kg/(s.m)] ',
     &              'w_M[m/s] w_Z[m/s] w_P[m/s] c_p[J/(mol.K)] ',
     &              'c_v[J/(mol.K)] gamma'

C     Second loop for gravity related parameters and printing to file.
      do iz=1, NSAMPLES+1
C       Gravity. [m.s^(-2)].
        Gravity(iz)=GMass/((radius+z(iz))**2.0)
C       Scale height. [m].
        Htab(iz)=0.0
     &           -1.0/((log(rhoat(iz+1))
     &           -log(rhoat(iz)))/(z(iz+1)-z(iz)))
C       Brunt–Väisälä frequency from TODO.
        Nsq1=0.0-((gammatab(iz)-1)/gammatab(iz))*Gravity(iz)*
     &       (log(rhoat(iz+1))-log(rhoat(iz)))/(z(iz+1)-z(iz))
C       Brunt–Väisälä frequency from Lighthill, "Waves in Fluids", page 295, Eq. (49).
C       In rad^2 / s^2.
        Nsq2=(gammatab(iz)-1)*(Gravity(iz)/v(iz))**2.0
        Nsqtab(iz)=Nsq2
C       Print to file.
        write(2012,*) z(iz),rhoat(iz),T(iz),v(iz),P(iz),
     &                Htab(iz),Gravity(iz),Nsqtab(iz),
     &                Kappatab(iz),MUtab(iz),MUvoltab(iz),
     &                Wind(1,iz),Wind(2,iz),Wind(3,iz),
     &                Cp(iz),Cv(iz),gammatab(iz)
      enddo
      close(2012)
      end

C ****************************************************************
C * Subroutines for computation of attenuation and of some other *
C * parameters.                                                  *
C ****************************************************************

C     ****************************************************************
C     * alclass                                                      *
C     ****************************************************************
      function alclass(f, rho, v, MU, K, gamma)
C     - X. Bonnin - 18/07/03
C     - 'alclass' calculates the term of
C       classical relaxation absorption of sound  
C     - Henry E. Bass and James P. Chambers
C       'Absorption of sound in the Martian atmosphere',
C       J. Acoust. Soc. Am., Vol.109,No.6,p.3069-3071 (June 2001)
C     - f in Hz, rho in kg/(m**3), v in m/s, MU in kg/(m*s),
C       K in (J*kg)/(kmol*K*m*s), Cv in J/(kmol*K),
C       gamma without unit
      real*4 alclass
      real*4 f,rho,v,MU,K,gamma
      real*4 Cv,R
      PI=ACOS(-1.)
C     R in J/(kmol*K) 
      R=8.314*(1.e+3)
C     Cv in J/(kmol*K)
      Cv=R/(gamma-1.0)
C     alclass in Np/m
      alclass=((2.0*(PI*f)**2.0)/(rho*(v**3.0)))
     *   *((4.0/3.0)*MU+(gamma-1.0)*K/(gamma*Cv))
      end
      
C     ****************************************************************
C     * MUvolclass                                                   *
C     ****************************************************************
      function MUvolclass(MU, K, gamma, RHO, T, P)
C     - X. Bonnin - 18/07/03
C     - 'alclass' calculates the term of
C       classical relaxation absorption of sound  
C     - Henry E. Bass and James P. Chambers
C       'Absorption of sound in the Martian atmosphere',
C       J. Acoust. Soc. Am., Vol.109,No.6,p.3069-3071 (June 2001)
C     - f in Hz, rho in kg/(m**3), v in m/s, MU in kg/(m*s),
C       K in (J*kg)/(kmol*K*m*s), Cv in J/(kmol*K),
C       gamma without unit
      real*4 MUvolclass
      real*4 f,rho,v,MU,K,gamma,KapMol
      real*4 Cv,R
      PI=ACOS(-1.)
C     R in J/(kmol*K) 
      R=8.314*(1.e+3)
C     Kappa in (J*kg)/(kmol*K*m*s):
      KapMol=RHO*R*T*K/P
C     Cv in J/(kmol*K)
      Cv=R/(gamma-1.0)
C     MUvolclass in Pa.s
      MUvolclass=((4.0/3.0)*MU+(gamma-1.0)*KapMol/(gamma*Cv))
      end
      
C     ****************************************************************
C     * MU                                                           *
C     ****************************************************************
      function MU(rho, P, T, gamma)
C     - X. Bonnin - 18/07/03
C     - 'MU' calculates the coefficient of viscosity of Earth atmosphere (N2 approximation
C     - Following ECSS standards:
C       http://www.spenvis.oma.be/spenvis/ecss/ecss07/ecss07.html
C
C     - T in K 
C     - P in Pa
C     - rho in kg/(m**3)
C     - c speed of sound (m/s)
C     - Mu in kg/(m*s)
C       gamma without unit
      real*4 MU
      real*4 P,T,rho,gamma
      real*4 Kap,k,PI,c,L,da
C      real*4 beta,S,MU1
      PI=ACOS(-1.)
C Formula USSA76:
C      beta =1.458E-06
C      S =110.4
C      MU1=(beta*(T**1.5))/(T+S)
C Formula ECSS
C     Kap : ratio of specific heats Cp/Cv
      Kap=gamma
C     da : mean collision diameter of N2
      da=3.62e-10
C     k : Boltzman constant (J/K)
      k=1.380658e-23
C     L : mean free path (m)
      L=1.0/(sqrt(2.0)*PI*da*da*P/(k*T))
C     c : speed of sound (m/s)
      c=sqrt(Kap*P/rho)
C     MU in kg/(m*s) or Pa.s
      MU=(2.0/3.0)*L*rho*c*sqrt(2.0/(PI*Kap))
      end
      
C     ****************************************************************
C     * Kappa                                                        *
C     ****************************************************************
      function Kappa(T)
C     - X. Bonnin - 18/07/03
C     -'K' calculates the thermal conductivity in Earth atmosphere
C     - Henry E. Bass and James P. Chambers
C       'Absorption of sound in the Martian atmosphere',
C       J. Acoust. Soc. Am., Vol.109,No.6,p.3069-3071 (June 2001)
C     - USSA76 (https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539.pdf).
C     - T in K
      real*4 Kappa, T
C     Empricial formula from USSA76, in J/(s*m*K).
      Kappa = (2.64638E-03)*(T**1.5)/(T+(245.4*(10.**(-12./T))))
      end
