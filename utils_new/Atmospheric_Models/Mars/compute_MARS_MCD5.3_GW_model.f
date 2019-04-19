C Authors:       X. Bonnin, R. F. Garcia.
C Description:   TODO.
C Last modified: See file metadata.
C Usage:         N/A.
C Notes:         N/A.

C Compilation line (including links to MCD5.2 and NetCDF routines):
C gfortran compute_MARS_MCD5.3_GW_model.f ../MCD5.3/mcd/call_mcd.F ../MCD5.3/mcd/julian.F ../MCD5.3/mcd/heights.F -I../MCD5.3/mcd -I../MCD5.3/netcdf-4.0.1/include -L../MCD5.3/netcdf-4.0.1/lib -lnetcdf -o compute_MARS_MCD5.3_GW_model

C Parametrised compilation line (including links to MCD5.2 and NetCDF routines):
C MCDPATH="/home/l.martire/Documents/software/MCD5.3"
C MCDPATH="/home/l.martire/Documents/software/MCD5.2"
C NETCDFPATH="/home/l.martire/Documents/software/netcdf-c-4.6.3"
C gfortran compute_MARS_MCD5.3_GW_model.f $MCDPATH/mcd/call_mcd.F $MCDPATH/mcd/julian.F $MCDPATH/mcd/heights.F -I$MCDPATH/mcd -I$NETCDFPATH/include -L$NETCDFPATH/lib -lnetcdf -o compute_MARS_MCD5.3_GW_model

C ****************************************************************
C * compute_MARS_MCD_AGW_model.                                  *
C ****************************************************************
      program compute_MARS_MCD_AGW_model



      integer imax
      parameter (imax=30000)

      character*60 file,filewind,outfile

      real*4 z(imax),rhoat(imax),v(imax),T(imax),P(imax),Ptest(imax)
      real*4 Cptest(imax),gammatest(imax),vtest(imax)
      real*4 dust(imax),densp(imax),fr(imax),Svib(imax)
      real*4 Wind(2,imax),Windshift(2),MUvolrottab(imax)
      real*4 MUtab(imax),MUvoltab(imax),Kappatab(imax)
      real*4 Gravity(imax),Nsqtab(imax),Htab(imax),Ncuttab(imax)
      real*4 Cpcoefs(5,11,3),tminrange(11,3),tmaxrange(11,3)
      real*4 molarmass(11),Ti,densmol,Dw(11),Dms(11),molarmasstot
      integer nbspec,ntrange(11),indspec(11),jext(11)
      real*4 Cp(imax),Cv(imax),gammatab(imax)
      real*4 omega,xdum1,xdum2,PI,dt,zsat,tt(imax),xdum3,xdum4
      real*4 ttz,xpress
      real*4 MU,MUvolclass,Kappa,MUvolrot
      real grav,Nsq,Nac,Fbr,Fac,Nsq1,Nsq2
      real GMass,radius,Hscale,gamma
      real dz

      complex Amp0,Ampz
      complex Amp(imax),Amp2(imax)

      integer i,j,izs,zsmin,zsmax,iz,i1,i2,i3,k,k1,k2,ix

      integer nlimod
      real*4 xdum,ltst,ls,xtime
      integer itop
      parameter (itop=28033)
      

      character*14 evename
      character*37 directory
      character*33 folder
      character*70 location

c MCD5.2 parameters
      integer nextvar
      parameter (nextvar=100)
      integer nmeanvar
      parameter (nmeanvar=5)
      integer zkey,hireskey,datekey,perturkey,scena,extvarkeys(nextvar)
      real xz,xlon,xlat,localtime,seedin,gwlength
      character*2 LT
      real*8 xdate
      character*80 dset
c  output variables
      integer ier
      real meanvar(nmeanvar),extvar(nextvar),seedout

c The sizes of folder and evename are depending on the event name,
c that's why they have to be changed for each event
      
      real*4 time(itop),D(itop),lati(itop),lon(itop),R(itop)
      real*4 WS(2,itop)
      real*4 zint(itop),tmp,h(itop),t1,t2,Tquake,f
      real*4 SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP(7),Tms(2)
      real*4 ap3(10),A(9),Vms,Pms,W(2)
      
      integer IYD,Z1,Z2,Zmin,Zmax,dh1,dh2,Iquake,day,year
      
c thermodynamic coefficients for species referenced in GRAM

cc O2
cc http://webbook.nist.gov/cgi/cbook.cgi?Name=O2&Units=SI&cTG=on#Thermo-Gas
      ispec=4
      molarmass(ispec)=32.0
      ntrange(ispec)=3
c Assume min temp = 80, but in fact 100 in NIST data base
      tminrange(ispec,1)=80
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


cc N2
cc http://webbook.nist.gov/cgi/cbook.cgi?Name=N2&Units=SI&cTG=on#Thermo-Gas
      ispec=2
      molarmass(ispec)=28.0
      ntrange(ispec)=3
c Assume min temp = 80, but in fact 100 in NIST data base
      tminrange(ispec,1)=80
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

cc He
cc http://webbook.nist.gov/cgi/cbook.cgi?ID=C7440597&Units=SI&Mask=1#Thermo-Gas
      ispec=7
      molarmass(ispec)=4.0
      ntrange(ispec)=1
c Assume min temp = 80, but in fact 290 in NIST data base
      tminrange(ispec,1)=80
      tmaxrange(ispec,1)=6000

      Cpcoefs(1,ispec,1)=20.78603
      Cpcoefs(2,ispec,1)=4.850638e-10
      Cpcoefs(3,ispec,1)=-1.582916e-10
      Cpcoefs(4,ispec,1)=1.525102e-10
      Cpcoefs(5,ispec,1)=3.196347e-11


cc Ar
cc http://webbook.nist.gov/cgi/cbook.cgi?Name=Ar&Units=SI&cTG=on#Thermo-Gas
      ispec=3
      molarmass(ispec)=40.0
      ntrange(ispec)=1
c Assume min temp = 80, but in fact 290 in NIST data base
      tminrange(ispec,1)=80
      tmaxrange(ispec,1)=6000

      Cpcoefs(1,ispec,1)=20.78600
      Cpcoefs(2,ispec,1)=2.825911e-7
      Cpcoefs(3,ispec,1)=-1.464191e-7
      Cpcoefs(4,ispec,1)=1.092131e-8
      Cpcoefs(5,ispec,1)=-3.661371e-8

c NOT USED IN MARS GRAM
cc N
cc http://webbook.nist.gov/cgi/cbook.cgi?ID=C17778880&Units=SI&Mask=1#Thermo-Gas
c      ispec=8
c      molarmass(ispec)=14.0
c      ntrange(ispec)=1
c      tminrange(ispec,1)=290
c      tmaxrange(ispec,1)=6000

c      Cpcoefs(1,ispec,1)=21.13581
c      Cpcoefs(2,ispec,1)=-0.388842
c      Cpcoefs(3,ispec,1)=0.043545
c      Cpcoefs(4,ispec,1)=0.024685
c      Cpcoefs(5,ispec,1)=-0.025678



cc CO2
cc http://webbook.nist.gov/cgi/cbook.cgi?ID=C17778880&Units=SI&Mask=1#Thermo-Gas
      ispec=1
      molarmass(ispec)=44.0
      ntrange(ispec)=2
c Assume min temp = 80, but in fact 298 in NIST data base
      tminrange(ispec,1)=80
      tmaxrange(ispec,1)=1200
      tminrange(ispec,2)=1200
      tmaxrange(ispec,2)=6000

      Cpcoefs(1,ispec,1)=24.99735
      Cpcoefs(2,ispec,1)=55.18696
      Cpcoefs(3,ispec,1)=-33.69137
      Cpcoefs(4,ispec,1)=7.948387
      Cpcoefs(5,ispec,1)=-0.136638
      Cpcoefs(1,ispec,2)=58.16639
      Cpcoefs(2,ispec,2)=2.720074
      Cpcoefs(3,ispec,2)=-0.492289
      Cpcoefs(4,ispec,2)=0.038844
      Cpcoefs(5,ispec,2)=-6.447293

cc CO
cc http://webbook.nist.gov/cgi/cbook.cgi?Name=CO&Units=SI&cTG=on
      ispec=5
      molarmass(ispec)=28.0
      ntrange(ispec)=2
c Assume min temp = 80, but in fact 298 in NIST data base
      tminrange(ispec,1)=80
      tmaxrange(ispec,1)=1300
      tminrange(ispec,2)=1300
      tmaxrange(ispec,2)=6000

      Cpcoefs(1,ispec,1)=25.56759
      Cpcoefs(2,ispec,1)=6.096130
      Cpcoefs(3,ispec,1)=4.054656
      Cpcoefs(4,ispec,1)=-2.671301
      Cpcoefs(5,ispec,1)=0.131021
      Cpcoefs(1,ispec,2)=35.15070
      Cpcoefs(2,ispec,2)=1.300095
      Cpcoefs(3,ispec,2)=-0.205921
      Cpcoefs(4,ispec,2)=0.013550
      Cpcoefs(5,ispec,2)=-3.282780


cc H2
cc http://webbook.nist.gov/cgi/cbook.cgi?Name=H2&Units=SI&cTG=on
      ispec=8
      molarmass(ispec)=2.0
      ntrange(ispec)=3
c Assume min temp = 80, but in fact 298 in NIST data base
      tminrange(ispec,1)=80
      tmaxrange(ispec,1)=1000
      tminrange(ispec,2)=1000
      tmaxrange(ispec,2)=2500
      tminrange(ispec,3)=2500
      tmaxrange(ispec,3)=6000

      Cpcoefs(1,ispec,1)=33.066178
      Cpcoefs(2,ispec,1)=-11.363417
      Cpcoefs(3,ispec,1)=11.432816
      Cpcoefs(4,ispec,1)=-2.772874
      Cpcoefs(5,ispec,1)=-0.158558
      Cpcoefs(1,ispec,2)=18.563083
      Cpcoefs(2,ispec,2)=12.257357
      Cpcoefs(3,ispec,2)=-2.859786
      Cpcoefs(4,ispec,2)=0.268238
      Cpcoefs(5,ispec,2)=1.977990
      Cpcoefs(1,ispec,3)=43.413560
      Cpcoefs(2,ispec,3)=-4.293079
      Cpcoefs(3,ispec,3)=1.272428
      Cpcoefs(4,ispec,3)=-0.096876
      Cpcoefs(5,ispec,3)=-20.533862


cc H2O
cc http://webbook.nist.gov/cgi/cbook.cgi?Name=H2O&Units=SI&cTG=on#Thermo-Gas
      ispec=10
      molarmass(ispec)=18.0
      ntrange(ispec)=2
c Assume min temp = 80, but in fact 500 in NIST data base
      tminrange(ispec,1)=80
      tmaxrange(ispec,1)=1700
      tminrange(ispec,2)=1700
      tmaxrange(ispec,2)=6000

      Cpcoefs(1,ispec,1)=30.09200
      Cpcoefs(2,ispec,1)=6.832514
      Cpcoefs(3,ispec,1)=6.793435
      Cpcoefs(4,ispec,1)=-2.534480
      Cpcoefs(5,ispec,1)=0.082139
      Cpcoefs(1,ispec,2)=41.96426
      Cpcoefs(2,ispec,2)=8.622053
      Cpcoefs(3,ispec,2)=-1.499780
      Cpcoefs(4,ispec,2)=0.098119
      Cpcoefs(5,ispec,2)=-11.15764

cc O3
cc http://webbook.nist.gov/cgi/cbook.cgi?ID=C10028156&Units=SI&Mask=1#Thermo-Gas
      ispec=11
      molarmass(ispec)=3*16.0
      ntrange(ispec)=2
c Assume min temp = 80, but in fact 500 in NIST data base
      tminrange(ispec,1)=80
      tmaxrange(ispec,1)=1200
      tminrange(ispec,2)=1200
      tmaxrange(ispec,2)=6000

      Cpcoefs(1,ispec,1)=21.66157
      Cpcoefs(2,ispec,1)=79.86001
      Cpcoefs(3,ispec,1)=-66.02603
      Cpcoefs(4,ispec,1)=19.58363
      Cpcoefs(5,ispec,1)=-0.079251
      Cpcoefs(1,ispec,2)=57.81409
      Cpcoefs(2,ispec,2)=0.730941
      Cpcoefs(3,ispec,2)=-0.039253
      Cpcoefs(4,ispec,2)=0.002610
      Cpcoefs(5,ispec,2)=-3.560367



c Number of species for CP values are computed from formulas above
      nbspec=8
c O2
      indspec(1)=4
      jext(1)=62
c N2
      indspec(2)=2
      jext(2)=58
c Ar
      indspec(3)=3
      jext(3)=59
c CO2
      indspec(4)=1
      jext(4)=57
c CO
      indspec(5)=5
      jext(5)=60
c H2
      indspec(6)=8
      jext(6)=65
c H2O
      indspec(7)=10
      jext(7)=42
c O3
      indspec(8)=11
      jext(8)=63
         

      Amp0=(1.0,0.0)
      PI=ACOS(-1.)
      omega=2.0*PI*0.01

C Mars and Atmosphere parameters (planetary scientist companion)
      GMass=4.2828e13
      radius=3389920.0
      gamma=1.67

c Creation of the atmosphere and wind model
c
      

c create MCD5.2 atmosphere model
c  height above surface (m)
      zkey=3
c  high resolution topography (not GCM one)
      hireskey=1
c 0 + Earth date in juldian day, 1= Mars date in degrees Ls
      datekey=1
c perturbations on top of model (0=none)
      perturkey=0
c model scenario (1 = Climatology ave solar)
      scena=1
c extravariables requests (1 everywhere request everything)
      extvarkeys=1
c directory of model set:
      dset='/home/r.garcia/INSIGHT/Science/Mars_atmosphere_models/'
      dset(55:66)='MCD5.3/data/'
      write(*,*) dset
c SEED random number generator of perturbation
      seedin=0.0
c  vertical wavelength of gravity waves
      gwlength=0.0

c  INSIGHT coordinates
c      xlon=135.7
c      xlat=4.4
c INSIGHT coordinates used by Aymeric Spiga
      xlon=135.6
      xlat=4.5

c  Local time at xlon (only used if datekey=1, if not should be 0)
c  in martian hours      
c      localtime=0.0
c      localtime=12.0
c Date in Ls (solar longitude)
c http://www-mars.lmd.jussieu.fr/mars/time/martian_time.html
c of Earth date 01/10/2015
c      xdate=48.6
c of Earth date 01/10/2020
c      xdate=288.1

c read outputs of MCD at LS=295 for local time from 0 to 22      
	xdate=295
	do ilt=0,23
c	do ilt=0,0
		localtime=ilt
		if (ilt<10) then
		write(LT,'(a,i1)') '0',ilt
		else
		write(LT,'(i2)') ilt
		end if

c  Define size of verical profile:
c      nlimod=202
c      dz=1000.0
c For R2020
      nlimod=25000
      dz=10.0

c output file
       if (xdate<10) then
 101   format('MARS_MCD_AGW_model_INSIGHT_LS_',i1,'_LT_',a2,'h.dat')
       write(outfile,101) nint(xdate),LT
	else
	if (xdate<100) then
 102   format('MARS_MCD_AGW_model_INSIGHT_LS_',i2,'_LT_',a2,'h.dat')
       write(outfile,102) nint(xdate),LT
	else
 103   format('MARS_MCD_AGW_model_INSIGHT_LS_',i3,'_LT_',a2,'h.dat')
       write(outfile,103) nint(xdate),LT
	endif
       endif
	write(*,*) outfile
c      open(2013,file='MARS_MCD_AGW_model_INSIGHT_1_Oct_2016_0h.dat',
c      open(2013,file='MARS_MCD_AGW_model_R2020_1_Oct_2020_12h.dat',
      open(2013,file=outfile,
     .form='formatted')
      write(2013,*) '    z(iz),         rhoat(iz),       T(iz),       ',
     .   '        v(iz),       P(iz),',
     .   '        Htab(iz),       Gravity(iz),       Nsqtab(iz),    ',
     .   '     Ncuttab(iz),        ',
     .      'Kappatab(iz),      MUtab(iz),    MUvoltab(iz),  ',
     .      '  MUvolrottab(iz),      fr(iz),      Svib(iz),      ',
     .      '  Mer. Wind(1,iz),   Zon. Wind(2,iz),          ',
     .      'Cp(iz)    ,     Cv(iz),   gammatab(iz)'

      do iz=1,nlimod

c	write(*,*) iz,nlimod

c     altitude of the model computation
      z(iz)=real(iz-1)*dz
c  Donne to avoid zero winds and surface temperature
      if (z(iz)<1.0) then
	z(iz)=1.0
      endif
      xz=z(iz)

c	write(*,*) iz, z(iz),extvarkeys
	extvarkeys=1
c	write(*,*) zkey,xz,xlon,xlat,hireskey,
c     &      datekey,xdate,localtime,dset,scena,
c     &      perturkey,seedin,gwlength,extvarkeys,ier

      call call_mcd(zkey,xz,xlon,xlat,hireskey,
     &      datekey,xdate,localtime,dset,scena,
     &      perturkey,seedin,gwlength,extvarkeys,
     &      P(iz),rhoat(iz),T(iz),W(2),W(1),
     &      meanvar,extvar,seedout,ier)

c  MCD model provide mixing ratios in mol/mol (r=ni/(ntot-ni))




cc computation of Cp in J/mol*K
      Ti=T(iz)/1000.0
cc Compute total molar mass from mixing ratios
      molarmasstot=0.0
      do i=1,nbspec
	ispec=indspec(i)
        molarmasstot=molarmasstot+extvar(jext(i))*
     .  (molarmass(ispec)/1000.0)
      end do
cc    add atomic oxygen and atomic hydrogen
        molarmasstot=molarmasstot+(extvar(61))*
     .  (16.0/1000.0)
        molarmasstot=molarmasstot+(extvar(64))*
     .  (1.0/1000.0)

cc compute relative densities in moles
      densmol=rhoat(iz)/molarmasstot
      xdum=0.0
      do i=1,nbspec
	ispec=indspec(i)
        Dms(ispec)=(extvar(jext(i)))*densmol
        xdum=xdum+Dms(ispec)/densmol
      end do
cc    add atomic oxygen and atomic hydrogen to the molar density
      Dms(6)=(extvar(61))*densmol
      Dms(9)=(extvar(64))*densmol
      xdum=xdum+(Dms(6)+Dms(9))/densmol
c      write(*,*) 'test densities:',xdum

cc compute Heat capacity at constant pressure of the gas mixture
      Cp(iz)=0.0
      do i=1,nbspec
	ispec=indspec(i)
        indtemp=1
	do k=1,ntrange(ispec)
	  if ((T(iz).ge.tminrange(ispec,k)).and.
     .        (T(iz).le.tmaxrange(ispec,k))) then
	    indtemp=k
	  end if
	end do
cc in J/mol*K
	xdum=Cp(iz)
	Cp(iz)=Cp(iz)+(Dms(ispec)/densmol)*(Cpcoefs(1,ispec,indtemp) +
     .   Cpcoefs(2,ispec,indtemp)*Ti +
     .   Cpcoefs(3,ispec,indtemp)*Ti*Ti + 
     .   Cpcoefs(4,ispec,indtemp)*Ti*Ti*Ti +
     .   Cpcoefs(5,ispec,indtemp)/(Ti*Ti))
c      write(*,*) iz,ispec,(Dms(ispec)/densmol),Cp(iz),
c     . (Cp(iz)-xdum)/(Dms(ispec)/densmol),indtemp
      end do
cc    add atomic oxygen and atomic hydrogen
      Cp(iz)=Cp(iz)+(Dms(6)/densmol)*(5/2)*8.3145
      Cp(iz)=Cp(iz)+(Dms(9)/densmol)*(5/2)*8.3145
c      write(*,*) iz,(Dms(6)/densmol),(Dms(9)/densmol),Cp(iz),densmol
      Cptest(iz)=Cp(iz)
      Cp(iz)=extvar(8)*molarmasstot

cc computation of Cv in J/mol*K
      Cv(iz)=Cp(iz)-8.3145
      gammatest(iz)=Cp(iz)/Cv(iz)
      gammatab(iz)=extvar(9)

cc atmospheric pressure from perfect gas hypothesis (Pa)
cc      P(iz)=rhoat(iz)*8.3145*T(iz)/0.0289645
cc      P(iz)=(densmol/6.02214129e23)*(8.3145*1000)*T(iz)
      Ptest(iz)=(densmol)*(8.3145)*T(iz)
cc sound velocity (m/s)
      vtest(iz)=sqrt(gammatest(iz)*Ptest(iz)/rhoat(iz))
      v(iz)=sqrt(gammatab(iz)*P(iz)/rhoat(iz))
cc Meridional (1) and zonal (2) winds in m/s
      Wind(1,iz)=W(1)
      Wind(2,iz)=W(2)	
cc Dynamic viscosity in Pa.s
      MUtab(iz)=MU(rhoat(iz),P(iz),T(iz),gammatab(iz))
cc atmosphere thermal conductivity in (J*kg)/(mol*K*m*s)
      Kappatab(iz)=Kappa(MUtab(iz),gammatab(iz),T(iz))
cc Volumic viscosity in Pa.s
      MUvoltab(iz)=MUvolclass(MUtab(iz),Kappatab(iz),gammatab(iz))
cc Rotational volumic viscosity in Pa.s
cc to be added to volumic viscosity to compute attenuation effects
      MUvolrottab(iz)=MUvolrot(MUtab(iz),T(iz))
cc Vibrational attenuation parameters
      call Attenuation_vibration(Dms(1)/densmol,Dms(2)/densmol,
     .     Dms(3)/densmol,Dms(10)/densmol,
     .     MUtab(iz),P(iz),T(iz),gammatab(iz),fr(iz),Svib(iz))


c Compare different parameters both computed by program and provided by model
c Pressure
c	write(*,*) iz,z(iz),'Pressure',Ptest(iz),P(iz)
c Air vissosity
c	write(*,*) iz,'Viscosity',MUtab(iz),extvar(54)
c Cp
c	write(*,*) iz,z(iz),molarmasstot,'CP',Cptest(iz),
c     .  extvar(8),extvar(8)*molarmasstot
c Gamma
c	write(*,*) iz,extvar(53),'Gamma',gammatest(iz),extvar(9)
c sound velocity
c	write(*,*) iz,'Sound speed',vtest(iz),
c     .  v(iz)

        end do

c second loop for gravity related parameters
        do iz=1,(nlimod-1)
cc gravity
c      write(*,*) GMass,radius,z(iz)
      Gravity(iz)=GMass/((radius+z(iz))**2.0)
c      write(*,*) grav,rhoat(iz),rhoat(iz-1),z(iz),z(iz-1)
      Htab(iz)=0.0
     .-1.0/((log(rhoat(iz+1))-log(rhoat(iz)))/(z(iz+1)-z(iz)))
c      Nsq1=0.0-((gammatab(iz)-1)/gammatab(iz))*Gravity(iz)*
c     .(log(rhoat(iz+1))-log(rhoat(iz)))/(z(iz+1)-z(iz))

c General formulation formulation of Nsquare 
c From lighthill (Waves in fluids) eq(12) p287
      Nsq1=0.0
     . -Gravity(iz)*(log(rhoat(iz+1))-log(rhoat(iz)))/(z(iz+1)-z(iz))
     . -Gravity(iz)*Gravity(iz)/(v(iz)*v(iz))

c Nsquare From lighthill (Waves in fluids) eq(49) p295
      Nsq2=(gammatab(iz)-1)*(Gravity(iz)/v(iz))**2.0
      Nsqtab(iz)=Nsq1
c Acoustic cut off frequency
      Ncuttab(iz)=v(iz)/(2*Htab(iz))
c      write(*,*) z(iz),Gravity(iz),Nsq1,Nsq2,Htab(iz),Ncuttab(iz)
    
      write(2013,*) z(iz),rhoat(iz),T(iz),v(iz),P(iz),
     .              Htab(iz),Gravity(iz),Nsqtab(iz),Ncuttab(iz),
     .              Kappatab(iz),MUtab(iz),MUvoltab(iz),
     .              MUvolrottab(iz),fr(iz),Svib(iz),
     .              Wind(1,iz),Wind(2,iz),
     .              Cp(iz),Cv(iz),gammatab(iz)

        end do
      
      close(2012)
      close(2013)

	write(*,*) outfile
	write(*,*) LT

c END Loop on local times:
	end do


      end

C******* Subroutines for computation of attenuation + other parameters
c **********

      function alclass(f,rho,v,MU,K,gamma)
C Modified for Mars assuming formulas of Bass and Chambers (2001)
c     - X. Bonnin - 18/07/03
c     - 'alclass' calculates the term of
c       classical relaxation absorption of sound  
c     - Henry E. Bass and James P. Chambers
c       'Absorption of sound in the Martian atmosphere',
c       J. Acoust. Soc. Am., Vol.109,No.6,p.3069-3071 (June 2001)
c     - f in Hz, rho in kg/(m**3), v in m/s, MU in kg/(m*s),
c       K in (J*kg)/(kmol*K*m*s), Cv in J/(kmol*K),
c       gamma without unit

      real*4 alclass
      real*4 f,rho,v,MU,K,gamma
      real*4 Cv,R

      PI=ACOS(-1.)

c     R in J/(kmol*K) 
c      R=8.314*(1.e+3)
c     Cv in J/(kmol*K)
c      Cv=R/(gamma-1.0)
c     alclass in Np/m
c      alclass=((2.0*(PI*f)**2.0)/(rho*(v**3.0)))
c     *   *((4.0/3.0)*MU+(gamma-1.0)*K/(gamma*Cv))

c     R in J/(mol*K) 
      R=8.314
c     Cv in J/(mol*K)
      Cv=R/(gamma-1.0)
c     alclass in Np/m
      alclass=((2.0*(PI*f)**2.0)/(rho*(v**3.0)))
     *   *((4.0/3.0)*MU+(gamma-1.0)*K/(gamma*Cv))

      end
      
c **********
      subroutine Attenuation_vibration(Xco2,Xn2,Xar,Xh2o,
     .           Mu,P,T,gamma,fr,S)

C Modified for Mars assuming formulas of Bass and Chambers (2001)
c     - Henry E. Bass and James P. Chambers
c       'Absorption of sound in the Martian atmosphere',
c       J. Acoust. Soc. Am., Vol.109,No.6,p.3069-3071 (June 2001)
c     - Xco2,Xn2,Xar,Xh2o molar fractions (no units)
c     - f in Hz, rho in kg/(m**3), v in m/s, Mu in kg/(m*s),
c       K in (J*kg)/(kmol*K*m*s), Cv in J/(kmol*K),
c       gamma without unit

      real*4 Mu,P,T
      real*4 f,rho,v,gamma
      real*4 Cv,R
      real*4 Xco2,Xn2,Xar,Xh2o
      real*4 kco2,kn2,kar,kh2o
      real*4 kk,tauvt,tauvs
      real*4 Cprime,theta
      real*4 Cpinf,Cp0,Cvinf,Cv0,Cp
      real*4 fr,S

      PI=ACOS(-1.)
      theta=960.0
c     R in J/(mol*K) 
      R=8.314
c not sure of 4 lines below !!! => Cpinf = Cp0
      Cvinf=(5./2.)*R
      Cpinf=Cvinf*gamma
      Cv0=(5./2.)*R
      Cp0=Cv0*gamma
      Cp=(gamma/(gamma-1))*R


      kco2=((0.219*P)/(Mu))*(exp(-60.75/((T**(1./3.)))))
      kn2=((1.44*P)/(Mu))*exp(-78.29/((T**(1./3.))))
      kar=kn2
      kh2o=((6.e-2)*P)/Mu

      kk=(Xco2*kco2)+(Xn2*kn2)+(Xar*kar)+(Xh2o*kh2o)
      tauvt=1.0/(kk*(1.0-exp(-theta/T)))
      tauvs=(Cpinf/Cp)*tauvt

c	write(*,*) Cpinf,Cp,tauvt,kk,theta,T
c        write(*,*) Mu,P,T

      fr=1.0/(2.0*PI*tauvs)

      Cprime=R*(Xco2*((theta/T)**2))
     .   *exp(-theta/T)/((1-exp(-theta/T))**2)

      S=(Cprime*R)/(Cpinf*(Cvinf+Cprime))

      end

c **********   

      function MUvolclass(MU,K,gamma)

C Modified for Mars assuming formulas of Bass and Chambers (2001)
c     - X. Bonnin - 18/07/03
c     - 'alclass' calculates the term of
c       classical relaxation absorption of sound  
c     - Henry E. Bass and James P. Chambers
c       'Absorption of sound in the Martian atmosphere',
c       J. Acoust. Soc. Am., Vol.109,No.6,p.3069-3071 (June 2001)
c     - f in Hz, rho in kg/(m**3), v in m/s, MU in kg/(m*s),
c       K in (J*kg)/(kmol*K*m*s), Cv in J/(kmol*K),
c       gamma without unit

      real*4 MUvolclass
      real*4 f,rho,v,MU,K,gamma
      real*4 Cv,R

      PI=ACOS(-1.)

c     R in J/(kmol*K) 
c      R=8.314*(1.e+3)
c     Cv in J/(kmol*K)
c      Cv=R/(gamma-1.0)
c     MUvolclass in Pa.s
c      MUvolclass=((4.0/3.0)*MU+(gamma-1.0)*K/(gamma*Cv))

c     R in J/(mol*K) 
      R=8.314
c     Cv in J/(mol*K)
      Cv=R/(gamma-1.0)
c     MUvolclass in Pa.s
      MUvolclass=((4.0/3.0)*MU+(gamma-1.0)*K/(gamma*Cv))

      end

c **********   
      function MUvolrot(MU,T)

C computed for Mars assuming formulas of Bass and Chambers (2001)
c     alpha_rot=((2*pi*pi*f*f)/(gamma*P*c))*MUvolrot
c     So, MUvolrot can simply be added to MUvolclass before computing
c     Attenuation parameters
c     - Henry E. Bass and James P. Chambers
c       'Absorption of sound in the Martian atmosphere',
c       J. Acoust. Soc. Am., Vol.109,No.6,p.3069-3071 (June 2001)
c     - f in Hz, rho in kg/(m**3), v in m/s, MU in kg/(m*s),
c       Cp in J/(mol*K),
c       gamma without unit

      real*4 MUvolrot,T
      real*4 f,rho,v,MU,K,gamma
      real*4 Cp,R,Zrot

      PI=ACOS(-1.)

c     rotational collision number:
      Zrot=61.1*exp(-16.8/((T**(1./3.))))
c     exclude vibrational effects in these formulas
      gamma=7.0/5.0
c     R in J/(mol*K) 
      R=8.314
c     Cv in J/(mol*K)
      Cp=(7.0/2.0)*R
c     MUvolrot in Pa.s
      MUvolrot=MU*(gamma*(gamma-1.0)*R/(1.25*Cp))*Zrot

c	write(*,*) MUvolrot,gamma,Cp,Zrot,MU
      end
      


c **********
      function MU(rho,P,T,gamma)
C Modified for Mars assuming formulas of Bass and Chambers (2001)
c     - X. Bonnin - 18/07/03
c     - 'MU' calculates the coefficient of viscosity of Earth atmosphere (N2 approximation
c     - Following ECSS standards:
c       http://www.spenvis.oma.be/spenvis/ecss/ecss07/ecss07.html
c
c     - T in K 
c     - P in Pa
c     - rho in kg/(m**3)
c     - c : speed of sound (m/s)
c     - Mu in kg/(m*s)
c       gamma without unit

      real*4 MU
      real*4 P,T,rho,gamma
      real*4 Kap,k,PI,c,L,da
      real*4 beta,S,MU1

      PI=ACOS(-1.)

c Formula USSA76:
c      beta =1.458E-06
c      S =110.4
c      MU1=(beta*(T**1.5))/(T+S)
      

c Formula ECSS
c     Kap : ratio of specific heats Cp/Cv
c      Kap=gamma
c     da : mean collision diameter of N2
c      da=3.62e-10
c     k : Boltzman constant (J/K)
c      k=1.380658e-23
c     L : mean free path (m)
c      L=1.0/(sqrt(2.0)*PI*da*da*P/(k*T))
c     c : speed of sound (m/s)
c      c=sqrt(Kap*P/rho)
c     MU in kg/(m*s) or Pa.s
c      MU=(2.0/3.0)*L*rho*c*sqrt(2.0/(PI*Kap))

c Formula of Sutherland equation for CO2 (Bass and chambers, 2001)
c in Pa.s
      beta=1.49e-6
      S=217.0
      MU=(beta*(T**(1.0/2.0)))/(1.0+(S/T))

      end
      
c **********

      function Kappa(Mu,gamma,T)
C Modified for Mars assuming formulas of Bass and Chambers (2001)
c      function Kappa(rho,P,T)
c     - X. Bonnin - 18/07/03
c     -'K' calculates the thermal conductivity in Earth atmosphere
c     - Henry E. Bass and James P. Chambers
c       'Absorption of sound in the Martian atmosphere',
c       J. Acoust. Soc. Am., Vol.109,No.6,p.3069-3071 (June 2001)
c     - T in K, Cv in J/(kmol*K) 

      real*4 Kappa,PI
      real*4 rho,T
      real*4 Kt,R,M,Mu,Cv

      PI=ACOS(-1.)

c Formula USSA76
c     Kt (J/(m*s*K))
c      Kt=(2.64638E-03)*(T**1.5)/(T+(245.4*(10**(-12./T))))
c Molar mass in (kg/Kmol):
c      M=rho*(8.314*1000.0)*T/P
c     K in (J*kg)/(kmol*K*m*s) 
c      Kappa=M*Kt

c Formula of Eucken expression for CO2 (Bass and chambers, 2001)
c     R in J/(mol*K) 
      R=8.314
c     Cv in J/(kmol*K)
      Cv=R/(gamma-1.0)
c in (J*kg)/(mol*K*m*s)
      Kappa=(15.0*R*Mu/4.0)
     .   *((4.0*Cv)/(15.0*R)+3.0/5.0)


      end

c **********

         


