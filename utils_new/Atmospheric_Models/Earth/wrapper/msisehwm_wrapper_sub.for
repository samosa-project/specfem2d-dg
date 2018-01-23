C ****************************************************************
C *****************            HWM93             *****************
C ****************************************************************
       SUBROUTINE GWS5(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,W)
C      Horizontal wind model HWM93 covering all altitude regions
C      A. E. HEDIN  (1/25/93) (4/9/93)
C      Calling argument list made similar to GTS5 subroutine for
C       MSIS-86 density model and GWS4 for thermospheric winds.
C        IYD - YEAR AND DAY AS YYDDD
C        SEC - UT(SEC)  (Not important in lower atmosphere)
C        ALT - ALTITUDE(KM) 
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS)
C        F107A - 3 MONTH AVERAGE OF F10.7 FLUX (Use 150 in lower atmos.)
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY ( " )
C        AP - Two element array with
C             AP(1) = MAGNETIC INDEX(DAILY) (use 4 in lower atmos.)
C             AP(2)=CURRENT 3HR ap INDEX (used only when SW(9)=-1.)
C     Note:  Ut, Local Time, and Longitude are used independently in the
C            model and are not of equal importance for every situation.  
C            For the most physically realistic calculation these three
C            variables should be consistent.
C      OUTPUT
C        W(1) = MERIDIONAL (m/sec + Northward)
C        W(2) = ZONAL (m/sec + Eastward)
C          ADDITIONAL COMMENTS
C               TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SW)
C               WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1. 
C               FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
C               FOR THE FOLLOWING VARIATIONS
C               1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
C               3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
C               5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
C               7 - DIURNAL               8 - SEMIDIURNAL
C               9 - DAILY AP             10 - ALL UT/LONG EFFECTS
C              11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
C              13 - MIXED AP/UT/LONG     14 - TERDIURNAL
C              16 - ALL WINDF VAR        17 - ALL WZL VAR
C              18 - ALL UN1 VAR          19 - ALL WDZL VAR
C              24 - ALL B FIELDS (DIV)   25 - ALL C FIELDS (CURL)
C
C              To get current values of SW: CALL TRETRV(SW)
C
C             For example, to get zonal averages (no diurnal or
C             longitudinal variations) set SW(7),SW(8), SW(14),
C             and SW(10) equal to 0.  To just remove tidal variations 
C             set SW(7),SW(8), and SW(14) equal to 0.
      PARAMETER (MN1=5,MN2=14)
      DIMENSION AP(1),W(2),WINDF(2),WW(2),SV(25)
      DIMENSION WZL(2),WDZL(2)
      DIMENSION ZN1(MN1),UN1(MN1,2),UGN1(2,2)
      DIMENSION ZN2(MN2),UN2(MN2,2),UGN2(2,2)
      COMMON/PARMW5/PWB(200),PWC(200),PWBL(150),PWCL(150),PWBLD(150),
     $ PWCLD(150),PB12(150),PC12(150),PB13(150),PC13(150),
     $ PB14(150),PC14(150),PB15(150),PC15(150),
     $ PB15D(150),PC15D(150),PWP(100,26)
      COMMON/CSW/SW(25),ISW,SWC(25)               
      COMMON/HWMC/WBT(2),WCT(2)
      COMMON/DATW/ISD(3),IST(2),NAM(2)
      COMMON/DATIME/ISDATE(3),ISTIME(2),NAME(2)
      SAVE
      EXTERNAL INITW5,GWSBK5
      DATA S/.016/,ZL/200./,SV/25*1./,NNN/3/,MN2S/1/,MN2M/1/
      DATA ZN1/200.,150.,130.,115.,100./
      DATA ZN2/100.,90.,82.5,75.,67.5,60.,52.5,45.,37.5,30.,22.5,
     $ 15.,7.5,0/
C      Put identification data into common/datime/
      DO 1 I=1,3
        ISDATE(I)=ISD(I)
    1 CONTINUE
      DO 2 I=1,2
        ISTIME(I)=IST(I)
        NAME(I)=NAM(I)
    2 CONTINUE
      IF(ISW.NE.64999) CALL TSELEC(SV)
      YRD=IYD
      WW(1)=W(1)
      WW(2)=W(2)
C
      IF(ALT.LE.ZN1(MN1)) GOTO 50
C
C       EXOSPHERE WIND
      CALL GLBW5E(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PWB,PWC,WINDF)
      WINDF(1)=SW(16)*WINDF(1)
      WINDF(2)=SW(16)*WINDF(2)
C       WIND  AT ZL
      CALL GLBW5M(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PWBL,PWCL,WW)
      WZL(1)=(PWBL(1)*WINDF(1)+WW(1))*SW(17)*SW(18)
      WZL(2)=(PWBL(1)*WINDF(2)+WW(2))*SW(17)*SW(18)
      UN1(1,1)=WZL(1)
      UN1(1,2)=WZL(2)
C       WIND DERIVATIVE AT ZL
      WW(1)=0
      WW(2)=0
      CALL GLBW5M(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PWBLD,PWCLD,WW)
      WDZL(1)=(PWBLD(1)*WINDF(1)+WW(1))*SW(19)*SW(18)
      WDZL(2)=(PWBLD(1)*WINDF(2)+WW(2))*SW(19)*SW(18)
      UGN1(1,1)=WDZL(1)*S
      UGN1(1,2)=WDZL(2)*S
C
      IF(ALT.GE.ZL) GOTO 90
C
C        WIND AT ZN1(2) (150)
      CALL GLBW5M(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PB12,PC12,WW)
      UN1(2,1)=(PB12(1)*WINDF(1)+WW(1))*SW(18)
      UN1(2,2)=(PB12(1)*WINDF(2)+WW(2))*SW(18)
C        WIND AT ZN1(3) (130)
      CALL GLBW5M(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PB13,PC13,WW)
      UN1(3,1)=WW(1)*SW(18)
      UN1(3,2)=WW(2)*SW(18)
C        WIND AT ZN1(4) (115)
      CALL GLBW5M(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PB14,PC14,WW)
      UN1(4,1)=WW(1)*SW(18)
      UN1(4,2)=WW(2)*SW(18)
C
   50 CONTINUE
      MNN=MAX(1,MIN(MN2,NNN+1))
      IF(ALT.LT.ZN2(MNN)) GOTO 40
C
C        WIND AT ZN1(5) (100)
      CALL GLBW5M(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PB15,PC15,WW)
      UN1(5,1)=WW(1)*SW(18)
      UN1(5,2)=WW(2)*SW(18)
C         WIND DERIVATIVE AT ZN1(5) (100)
      CALL GLBW5M(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PB15D,PC15D,WW)
      UGN1(2,1)=WW(1)*SW(18)
      UGN1(2,2)=WW(2)*SW(18)
C
      IF(ALT.GE.ZN1(MN1)) GOTO 90
C
      UGN2(1,1)=UGN1(2,1)
      UGN2(1,2)=UGN1(2,2)
      UN2(1,1)=UN1(5,1)
      UN2(1,2)=UN1(5,2)
      GOTO 45
   40 CONTINUE
      UGN2(1,1)=1.E30
      UGN2(1,2)=1.E30
      UN2(1,1)=0
      UN2(1,2)=0
   45 CONTINUE
C
      DO 10 I=1,MN2
        IF(ALT.GT.ZN2(I)) GOTO 12
   10 CONTINUE
      I=MN2
   12 IZ=I
      MN2S=MAX(1,MIN(IZ-1,IZ-NNN))
      MN2E=MIN(MN2,MAX(MN2S+1,IZ-1+NNN))
      DO 20 I=MN2S,MN2E
        II=2*(I-2)+1
        IF(I.GT.1) THEN
          CALL GLBW5S(IYD,GLAT,GLONG,STL,PWP(1,II),PWP(1,II+1),WW)
          UN2(I,1)=WW(1)*SW(20)
          UN2(I,2)=WW(2)*SW(20)
        ENDIF
   20 CONTINUE
      MN2M=MN2E-MN2S+1
      UGN2(2,1)=1.E30
      UGN2(2,2)=1.E30
   90 CONTINUE
C       WIND AT ALTITUDE
      IF(W(1).NE.9898)
     $ W(1)= WPROF(ALT,ZL,S,WINDF(1),WZL(1),WDZL(1),
     $  MN1,ZN1,UN1(1,1),UGN1(1,1),MN2M,ZN2(MN2S),UN2(MN2S,1),UGN2(1,1))
      IF(W(2).NE.9898)
     $ W(2)= WPROF(ALT,ZL,S,WINDF(2),WZL(2),WDZL(2),
     $  MN1,ZN1,UN1(1,2),UGN1(1,2),MN2M,ZN2(MN2S),UN2(MN2S,2),UGN2(1,2))
      RETURN
C       Set number of nodes calculated each side of required altitude
C         to adjust profile accuracy vs efficiency
      ENTRY SETNW5(NNW)
      NNN=NNW
      END
C-----------------------------------------------------------------------
      FUNCTION WPROF(Z,ZL,S,UINF,ULB,ULBD,MN1,ZN1,UN1,UGN1,
     $   MN2,ZN2,UN2,UGN2)
      DIMENSION ZN1(MN1),UN1(MN1),UGN1(2),XS(15),YS(15),Y2OUT(15)
      DIMENSION ZN2(MN2),UN2(MN2),UGN2(2)
      SAVE
      IF(Z.GE.ZL) THEN
        X=S*(Z-ZL)
        F=EXP(-X)
        WPROF=UINF+(ULB-UINF)*F+(ULB-UINF+ULBD)*X*F
        RETURN
      ENDIF
      IF(Z.GE.ZN1(MN1).AND.Z.LT.ZN1(1)) THEN
        MN=MN1
        Z1=ZN1(1)
        Z2=ZN1(MN)
        ZDIF=Z2-Z1
        DO 10 K=1,MN
          XS(K)=(ZN1(K)-Z1)/ZDIF
          YS(K)=UN1(K)
   10   CONTINUE
        YD1=UGN1(1)*ZDIF
        YD2=UGN1(2)*ZDIF
        CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
C      Eq.
        X=(Z-Z1)/ZDIF
C      Eq. 
        CALL SPLINT(XS,YS,Y2OUT,MN,X,Y)
        WPROF=Y
        RETURN
      ENDIF
      IF(Z.LT.ZN2(1)) THEN
        MN=MN2
        Z1=ZN2(1)
        Z2=ZN2(MN)
        ZDIF=Z2-Z1
        DO 20 K=1,MN
          XS(K)=(ZN2(K)-Z1)/ZDIF
          YS(K)=UN2(K)
   20   CONTINUE
        YD1=UGN2(1)
        IF(UGN2(1).LT.1.E30) YD1=UGN2(1)*ZDIF
        YD2=UGN2(2)
        IF(UGN2(2).LT.1.E30) YD2=UGN2(2)*ZDIF
        CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
C      Eq.
        X=(Z-Z1)/ZDIF
C      Eq. 
        CALL SPLINT(XS,YS,Y2OUT,MN,X,Y)
        WPROF=Y
        RETURN
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GLBW5E(YRD,SEC,LAT,LONG,STL,F107A,F107,AP,PB,PC,WW)
      REAL LAT,LONG
      DIMENSION WB(2,15),WC(2,15),PB(200),PC(200),WW(2)
      DIMENSION AP(1)
      COMMON/CSW/SW(25),ISW,SWC(25)
      COMMON/HWMC/WBT(2),WCT(2)
C      COMMON/VPOLY/BT(20,20),BP(20,20),CSTL,SSTL,C2STL,S2STL,
C     $ C3STL,S3STL,IYR,DAY,DF,DFA,DFC,APD,APDF,APDFC,APT,SLT
      COMMON/VPOLY2/XVL,LVL,MVL,CLAT,SLAT,BT(20,20),BP(20,20)
      COMMON/LTCOMP/TLL,NSVL,CSTL,SSTL,C2STL,S2STL,C3STL,S3STL
      COMMON/LGCOMP/XLL,NGVL,CLONG,SLONG,C2LONG,S2LONG
      SAVE
      DATA DGTR/.017453/,SR/7.2722E-5/,HR/.2618/,DR/1.72142E-2/
      DATA NSW/14/,WB/30*0/,WC/30*0/
      DATA PB14/-1./,PB18/-1./
      DATA SW9/1./,LV/12/,MV/3/,NSV/3/,NGV/2/,PSET/3./
      G0(A)=(A-4.+(PB(26)-1.)*(A-4.+(EXP(-ABS(PB(25))*(A-4.))-1.)/
     * ABS(PB(25))))
C       CONFIRM PARAMETER SET
      IF(PB(100).EQ.0) PB(100)=PSET
      IF(PB(100).NE.PSET) THEN
        WRITE(6,900) PB(100),PC(100)
  900   FORMAT(1X,'WRONG PARAMETER SET FOR GLBW5E',3F10.1)
        STOP
      ENDIF
C
      DO 10 J=1,NSW
        WB(1,J)=0
        WB(2,J)=0
        WC(1,J)=0
        WC(2,J)=0
   10 CONTINUE
      IF(SW(9).GT.0) SW9=1.
      IF(SW(9).LT.0) SW9=-1.
      IYR = YRD/1000.
      DAY = YRD - IYR*1000.
      IF(XVL.NE.LAT.OR.LV.GT.LVL.OR.MV.GT.MVL) THEN
        SLAT=SIN(DGTR*LAT)
        CLAT=COS(DGTR*LAT)
        CALL VSPHR1(SLAT,CLAT,LV,MV,BT,BP,20)
        XVL=LAT
        LVL=LV
        MVL=MV
      ENDIF
      IF(TLL.NE.STL.OR.NSV.GT.NSVL)  THEN
        SSTL = SIN(HR*STL)
        CSTL = COS(HR*STL)
        S2STL = SIN(2.*HR*STL)
        C2STL = COS(2.*HR*STL)
        S3STL = SIN(3.*HR*STL)
        C3STL = COS(3.*HR*STL)
        TLL = STL
        NSVL=NSV
      ENDIF
      IF(DAY.NE.DAYL.OR.PB(14).NE.PB14) THEN
        CD14=COS(DR*(DAY-PB(14)))
C        SD14=SIN(DR*(DAY-PB(14)))
      ENDIF
      IF(DAY.NE.DAYL.OR.PB(18).NE.PB18) CD18=COS(2.*DR*(DAY-PB(18)))
      DAYL=DAY
      PB14=PB(14)
      PB18=PB(18)
      IF(XLL.NE.LONG) THEN
        SLONG=SIN(DGTR*LONG)
        CLONG=COS(DGTR*LONG)
        S2LONG=SIN(2.*DGTR*LONG)
        C2LONG=COS(2.*DGTR*LONG)
        XLL=LONG
        NGVL=2
      ENDIF
C       F10.7 EFFECT
      DF=F107-F107A
      DFA=F107A-150.
      DFC=DFA+PB(20)*DF
C       TIME INDEPENDENT
      F1B=1.+PB(22)*DFC*SWC(1)
      IF(WW(1).NE.9898) THEN
       WB(1,2)=(PB(2)*BT(3,1)+PB(3)*BT(5,1)+PB(23)*BT(7,1))*F1B
      ENDIF
      WB(2,2)=0.
      F1C=1.+PC(22)*DFC*SWC(1)
      WC(1,2)=0.
      IF(WW(2).NE.9898) THEN
       WC(2,2)=-(PC(2)*BT(2,1)+PC(3)*BT(4,1)+PC(23)*BT(6,1))*F1C
     $ -(PC(27)*BT(3,1)+PC(15)*BT(5,1)+PC(60)*BT(7,1)
     $ +PC(161)*BT(9,1)+PC(162)*BT(11,1)+PC(163)*BT(13,1))*F1C
      ENDIF
C       SYMMETRICAL ANNUAL
C       SYMMETRICAL SEMIANNUAL
      IF(WW(1).NE.9898) THEN
       WB(1,4)=(PB(17)*BT(3,1)+PB(31)*BT(5,1))*CD18
      ENDIF
      WB(2,4)=0
      WC(1,4)=0
      IF(WW(2).NE.9898) THEN
       WC(2,4)=-(PC(17)*BT(2,1)+PC(31)*BT(4,1))*CD18
      ENDIF
C       ASYMMETRICAL ANNUAL
      F5B=1.+PB(48)*DFC*SWC(1)
      IF(WW(1).NE.9898) THEN
       WB(1,5)=(PB(10)*BT(2,1)+PB(11)*BT(4,1))*CD14*F5B
      ENDIF
      WB(2,5)=0
      F5C=1.+PC(48)*DFC*SWC(1)
      WC(1,5)=0
      IF(WW(2).NE.9898) THEN
       WC(2,5)=-(PC(10)*BT(3,1)+PC(11)*BT(5,1))*CD14*F5C
      ENDIF
C       ASYMMETRICAL SEMIANNUAL
C         none
C       DIURNAL      
      IF(SW(7).EQ.0) GOTO 200
      F7B=1.+PB(50)*DFC*SWC(1)
      F75B=1.+PB(83)*DFC*SWC(1)
      IF(WW(1).NE.9898) THEN
       WB(1,7)=(PB(7)*BT(2,2)+PB(8)*BT(4,2)+PB(29)*BT(6,2)
     $ +PB(142)*BT(8,2)+PB(144)*BT(10,2)
     $  +PB(182)*BT(3,2)+PB(184)*BT(5,2)
     $  )*SSTL*F7B
     $ +(PB(13)*BT(3,2)+PB(146)*BT(5,2))
     $    *CD14*SSTL*F75B*SWC(5)
     $ +(PB(171)*BT(2,2)+PB(173)*BT(4,2))
     $    *CD18*SSTL*F75B*SWC(4)
     $ + (PB(4)*BT(2,2)+PB(5)*BT(4,2)+PB(28)*BT(6,2)
     $ +PB(141)*BT(8,2)+PB(143)*BT(10,2)
     $  +PB(181)*BT(3,2)+PB(183)*BT(5,2)
     $  )*CSTL*F7B
     $ +(PB(12)*BT(3,2)+PB(145)*BT(5,2))
     $      *CD14*CSTL*F75B*SWC(5)
     $ +(PB(170)*BT(2,2)+PB(172)*BT(4,2))
     $    *CD18*CSTL*F75B*SWC(4)
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,7)=-(PB(4)*BP(2,2)+PB(5)*BP(4,2)+PB(28)*BP(6,2)
     $   +PB(141)*BP(8,2)+PB(143)*BP(10,2)
     $   +PB(181)*BP(3,2)+PB(183)*BP(5,2)
     $  )*SSTL*F7B
     $ -(PB(12)*BP(3,2)+PB(145)*BP(5,2))
     $    *CD14*SSTL*F75B*SWC(5)
     $ -(PB(170)*BP(2,2)+PB(172)*BP(4,2))
     $    *CD18*SSTL*F75B*SWC(4)
     $ + (PB(7)*BP(2,2)+PB(8)*BP(4,2)+PB(29)*BP(6,2)
     $   +PB(142)*BP(8,2)+PB(144)*BP(10,2)
     $   +PB(182)*BP(3,2)+PB(184)*BP(5,2)
     $  )*CSTL*F7B
     $ +(PB(13)*BP(3,2)+PB(146)*BP(5,2))
     $    *CD14*CSTL*F75B*SWC(5)
     $ +(PB(171)*BP(2,2)+PB(173)*BP(4,2))
     $    *CD18*CSTL*F75B*SWC(4)
      ENDIF
      F7C=1.+PC(50)*DFC*SWC(1)
      F75C=1.+PC(83)*DFC*SWC(1)
      IF(WW(1).NE.9898) THEN
       WC(1,7)=-(PC(4)*BP(3,2)+PC(5)*BP(5,2)+PC(28)*BP(7,2)
     $   +PC(141)*BP(9,2)+PC(143)*BP(11,2)
     $   +PC(181)*BP(2,2)+PC(183)*BP(4,2)+PC(185)*BP(6,2)
     $   +PC(187)*BP(8,2)+PC(189)*BP(10,2)
     $  )*SSTL*F7C
     $ -(PC(12)*BP(2,2)+PC(145)*BP(4,2))
     $    *CD14*SSTL*F75C*SWC(5)
     $ -(PC(170)*BP(3,2)+PC(172)*BP(5,2))
     $    *CD18*SSTL*F75C*SWC(4)
     $ +(PC(7)*BP(3,2)+PC(8)*BP(5,2)+PC(29)*BP(7,2)
     $ +PC(142)*BP(9,2)+PC(144)*BP(11,2)
     $ +PC(182)*BP(2,2)+PC(184)*BP(4,2)+PC(186)*BP(6,2)
     $ +PC(188)*BP(8,2)+PC(190)*BP(10,2)
     $  )*CSTL*F7C
     $ +(PC(13)*BP(2,2)+PC(146)*BP(4,2))
     $     *CD14*CSTL*F75C*SWC(5)
     $ +(PC(171)*BP(3,2)+PC(173)*BP(5,2))
     $    *CD18*CSTL*F75C*SWC(4)
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,7)=-(PC(7)*BT(3,2)+PC(8)*BT(5,2)+PC(29)*BT(7,2)
     $ +PC(142)*BT(9,2)+PC(144)*BT(11,2)
     $ +PC(182)*BT(2,2)+PC(184)*BT(4,2)+PC(186)*BT(6,2)
     $ +PC(188)*BT(8,2)+PC(190)*BT(10,2)
     $  )*SSTL*F7C
     $ -(PC(13)*BT(2,2)+PC(146)*BT(4,2))
     $    *CD14*SSTL*F75C*SWC(5)
     $ -(PC(171)*BT(3,2)+PC(173)*BT(5,2))
     $    *CD18*SSTL*F75C*SWC(4)
     $ -(PC(4)*BT(3,2)+PC(5)*BT(5,2)+PC(28)*BT(7,2)
     $ +PC(141)*BT(9,2)+PC(143)*BT(11,2)
     $ +PC(181)*BT(2,2)+PC(183)*BT(4,2)+PC(185)*BT(6,2)
     $ +PC(187)*BT(8,2)+PC(189)*BT(10,2)
     $  )*CSTL*F7C
     $ -(PC(12)*BT(2,2)+PC(145)*BT(4,2))
     $    *CD14*CSTL*F75C*SWC(5)
     $ -(PC(170)*BT(3,2)+PC(172)*BT(5,2))
     $    *CD18*CSTL*F75C*SWC(4)
      ENDIF
  200 CONTINUE
C       SEMIDIURNAL
      IF(SW(8).EQ.0) GOTO 210
      F8B=1.+PB(90)*DFC*SWC(1)
      IF(WW(1).NE.9898) THEN
       WB(1,8)=(PB(9)*BT(3,3)+PB(43)*BT(5,3)
     $   +PB(111)*BT(7,3)
     $   +(PB(34)*BT(4,3)+PB(148)*BT(6,3))*CD14*SWC(5)
     $   +(PB(134)*BT(3,3))*CD18*SWC(4) 
     $   +PB(152)*BT(4,3)+PB(154)*BT(6,3)+PB(156)*BT(8,3)
     $   +PB(158)*BT(10,3)
     $  )*S2STL*F8B
     $ +(PB(6)*BT(3,3)+PB(42)*BT(5,3)
     $   +PB(110)*BT(7,3)
     $   +(PB(24)*BT(4,3)+PB(147)*BT(6,3))*CD14*SWC(5)
     $   +(PB(135)*BT(3,3))*CD18*SWC(4)
     $   +PB(151)*BT(4,3)+PB(153)*BT(6,3)+PB(155)*BT(8,3)
     $   +PB(157)*BT(10,3)
     $  )*C2STL*F8B
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,8)=-(PB(6)*BP(3,3)+PB(42)*BP(5,3)
     $   +PB(110)*BP(7,3)
     $   +(PB(24)*BP(4,3)+PB(147)*BP(6,3))*CD14*SWC(5)
     $   +(PB(135)*BP(3,3))*CD18*SWC(4)
     $   +PB(151)*BP(4,3)+PB(153)*BP(6,3)+PB(155)*BP(8,3)
     $   +PB(157)*BP(10,3)
     $  )*S2STL*F8B
     $   + (PB(9)*BP(3,3)+PB(43)*BP(5,3)
     $   +PB(111)*BP(7,3)
     $   +(PB(34)*BP(4,3)+PB(148)*BP(6,3))*CD14*SWC(5)
     $   +(PB(134)*BP(3,3))*CD18*SWC(4)
     $   +PB(152)*BP(4,3)+PB(154)*BP(6,3)+PB(156)*BP(8,3)
     $   +PB(158)*BP(10,3)
     $  )*C2STL*F8B
      ENDIF
      F8C=1.+PC(90)*DFC*SWC(1)
      IF(WW(1).NE.9898) THEN
       WC(1,8)=-(PC(6)*BP(4,3)+PC(42)*BP(6,3)
     $   +PC(110)*BP(8,3)
     $   +(PC(24)*BP(3,3)+PC(147)*BP(5,3))*CD14*SWC(5)
     $   +(PC(135)*BP(4,3))*CD18*SWC(4)
     $   +PC(151)*BP(3,3)+PC(153)*BP(5,3)+PC(155)*BP(7,3)
     $   +PC(157)*BP(9,3)
     $  )*S2STL*F8C
     $ +(PC(9)*BP(4,3)+PC(43)*BP(6,3)
     $   +PC(111)*BP(8,3)
     $   +(PC(34)*BP(3,3)+PC(148)*BP(5,3))*CD14*SWC(5)
     $   +(PC(134)*BP(4,3))*CD18*SWC(4)
     $   +PC(152)*BP(3,3)+PC(154)*BP(5,3)+PC(156)*BP(7,3)
     $   +PC(158)*BP(9,3)
     $  )*C2STL*F8C
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,8)=-(PC(9)*BT(4,3)+PC(43)*BT(6,3)
     $   +PC(111)*BT(8,3)
     $   +(PC(34)*BT(3,3)+PC(148)*BT(5,3))*CD14*SWC(5)
     $   +(PC(134)*BT(4,3))*CD18*SWC(4)
     $   +PC(152)*BT(3,3)+PC(154)*BT(5,3)+PC(156)*BT(7,3)
     $   +PC(158)*BT(9,3)
     $  )*S2STL*F8C
     $ - (PC(6)*BT(4,3)+PC(42)*BT(6,3)
     $   +PC(110)*BT(8,3)
     $   +(PC(24)*BT(3,3)+PC(147)*BT(5,3))*CD14*SWC(5)
     $   +(PC(135)*BT(4,3))*CD18*SWC(4)
     $   +PC(151)*BT(3,3)+PC(153)*BT(5,3)+PC(155)*BT(7,3)
     $   +PC(157)*BT(9,3)
     $  )*C2STL*F8C
      ENDIF
  210 CONTINUE
C        TERDIURNAL
      IF(SW(14).EQ.0) GOTO 220
      F14B=1.
      IF(WW(1).NE.9898) THEN
       WB(1,14)=(PB(40)*BT(4,4)+PB(149)*BT(6,4)
     $   +PB(114)*BT(8,4)
     $   +(PB(94)*BT(5,4)+PB(47)*BT(7,4))*CD14*SWC(5)
     $  )*S3STL*F14B
     $ + (PB(41)*BT(4,4)+PB(150)*BT(6,4)
     $   +PB(115)*BT(8,4)
     $   +(PB(95)*BT(5,4)+PB(49)*BT(7,4))*CD14*SWC(5)
     $  )*C3STL*F14B
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,14)=-(PB(41)*BP(4,4)+PB(150)*BP(6,4)
     $   +PB(115)*BP(8,4)
     $   +(PB(95)*BP(5,4)+PB(49)*BP(7,4))*CD14*SWC(5)
     $  )*S3STL*F14B
     $ + (PB(40)*BP(4,4)+PB(149)*BP(6,4)
     $   +PB(114)*BP(8,4)
     $   +(PB(94)*BP(5,4)+PB(47)*BP(7,4))*CD14*SWC(5)
     $  )*C3STL*F14B
      ENDIF
      F14C=1.
      IF(WW(1).NE.9898) THEN
       WC(1,14)=-(PC(41)*BP(5,4)+PC(150)*BP(7,4)
     $   +PC(115)*BP(9,4)
     $   +(PC(95)*BP(4,4)+PC(49)*BP(6,4))*CD14*SWC(5)
     $  )*S3STL*F14C
     $ + (PC(40)*BP(5,4)+PC(149)*BP(7,4)
     $   +PC(114)*BP(9,4)
     $   +(PC(94)*BP(4,4)+PC(47)*BP(6,4))*CD14*SWC(5)
     $  )*C3STL*F14C
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,14)=-(PC(40)*BT(5,4)+PC(149)*BT(7,4)
     $   +PC(114)*BT(9,4)
     $   +(PC(94)*BT(4,4)+PC(47)*BT(6,4))*CD14*SWC(5)
     $  )*S3STL*F14C
     $ - (PC(41)*BT(5,4)+PC(150)*BT(7,4)
     $   +PC(115)*BT(9,4)
     $   +(PC(95)*BT(4,4)+PC(49)*BT(6,4))*CD14*SWC(5)
     $  )*C3STL*F14C
      ENDIF
  220 CONTINUE
C        MAGNETIC ACTIVITY
      IF(SW(9).EQ.0.) GOTO 40
      IF(SW9.EQ.-1.) GOTO 30
C           daily AP
      APD=AP(1)-4.
      APDF=(APD+(PB(45)-1.)*(APD+(EXP(-PB(44)*APD)-1.)/PB(44)))
C      APDFC=(APD+(PC(45)-1.)*(APD+(EXP(-PC(44)*APD)-1.)/PC(44)))
      APDFC=APDF
      IF(APD.EQ.0.) GOTO 40
      IF(WW(1).NE.9898) THEN
       WB(1,9)=(PB(46)*BT(3,1)+PB(35)*BT(5,1)+PB(33)*BT(7,1))*APDF
     $  +(PB(175)*BT(3,3)+PB(177)*BT(5,3))*S2STL*APDF
     $  +(PB(174)*BT(3,3)+PB(176)*BT(5,3))*C2STL*APDF
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,9)=0                                              
     $  -(PB(174)*BP(3,3)+PB(176)*BP(5,3))*S2STL*APDF
     $  +(PB(175)*BP(3,3)+PB(177)*BP(5,3))*C2STL*APDF
      ENDIF
      IF(WW(1).NE.9898) THEN
       WC(1,9)=SWC(7)*WC(1,7)*PC(122)*APDFC
     $  -(PC(174)*BP(4,3)+PC(176)*BP(6,3))*S2STL*APDFC
     $  +(PC(175)*BP(4,3)+PC(177)*BP(6,3))*C2STL*APDFC
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,9)=-(PC(46)*BT(2,1)+PC(35)*BT(4,1)+PC(33)*BT(6,1))*APDFC
     $ +SWC(7)*WC(2,7)*PC(122)*APDFC
     $ -(PC(175)*BT(4,3)+PC(177)*BT(6,3))*S2STL*APDFC
     $ -(PC(174)*BT(4,3)+PC(176)*BT(6,3))*C2STL*APDFC
      ENDIF
      GO TO 40
   30 CONTINUE
      IF(PB(25).LT.1.E-4) PB(25)=1.E-4
c Problem not solved because AP is of dimension 1
c     APT=G0(AP(2))
       APT=0
      IF(APT.EQ.0) GOTO 40
      IF(WW(1).NE.9898) THEN
       WB(1,9)=(PB(97)*BT(3,1)+PB(55)*BT(5,1)+PB(51)*BT(7,1))*APT
     $  +(PB(160)*BT(3,3)+PB(179)*BT(5,3))*S2STL*APT
     $  +(PB(159)*BT(3,3)+PB(178)*BT(5,3))*C2STL*APT
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,9)=0
     $  -(PB(159)*BP(3,3)+PB(178)*BP(5,3))*S2STL*APT
     $  +(PB(160)*BP(3,3)+PB(179)*BP(5,3))*C2STL*APT
      ENDIF
      IF(WW(1).NE.9898) THEN
       WC(1,9)=SWC(7)*WC(1,7)*PC(129)*APT
     $  -(PC(159)*BP(4,3)+PC(178)*BP(6,3))*S2STL*APT
     $  +(PC(160)*BP(4,3)+PC(179)*BP(6,3))*C2STL*APT
      ENDIF
      IF(WW(2).NE.9898) THEN
      WC(2,9)=-(PC(97)*BT(2,1)+PC(55)*BT(4,1)+PC(51)*BT(6,1))*APT
     $ +SWC(7)*WC(2,7)*PC(129)*APT
     $ -(PC(160)*BT(4,3)+PC(179)*BT(6,3))*S2STL*APT
     $ -(PC(159)*BT(4,3)+PC(178)*BT(6,3))*C2STL*APT
      ENDIF
  40  CONTINUE
      IF(SW(10).EQ.0) GOTO 49
C        LONGITUDINAL
      DBASY1=1.+PB(199)*SLAT
      DBASY2=1.+PB(200)*SLAT
      F11B=1.+PB(81)*DFC*SWC(1)
      IF(SW(11).EQ.0) GOTO 230
      IF(WW(1).NE.9898) THEN
       WB(1,11)=(PB(91)*BT(3,2)+PB(92)*BT(5,2)+PB(93)*BT(7,2))
     $  *SLONG*DBASY1*F11B
     $ + (PB(65)*BT(3,2)+PB(66)*BT(5,2)+PB(67)*BT(7,2))
     $  *CLONG*DBASY1*F11B
     $  +(PB(191)*BT(3,3)+PB(193)*BT(5,3)+PB(195)*BT(7,3)
     $   +PB(197)*BT(9,3)
     $  )*S2LONG*DBASY2*F11B
     $ + (PB(192)*BT(3,3)+PB(194)*BT(5,3)+PB(196)*BT(7,3)
     $    +PB(198)*BT(9,3)
     $  )*C2LONG*DBASY2*F11B
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,11)=-(PB(65)*BP(3,2)+PB(66)*BP(5,2)+PB(67)*BP(7,2))
     $  *SLONG*DBASY1*F11B
     $ + (PB(91)*BP(3,2)+PB(92)*BP(5,2)+PB(93)*BP(7,2))
     $  *CLONG*DBASY1*F11B
     $ -(PB(192)*BP(3,3)+PB(194)*BP(5,3)+PB(196)*BP(7,3)
     $   +PB(198)*BP(9,3)
     $  )*S2LONG*DBASY2*F11B
     $ + (PB(191)*BP(3,3)+PB(193)*BP(5,3)+PB(195)*BP(7,3)
     $    +PB(197)*BP(9,3)
     $  )*C2LONG*DBASY2*F11B
      ENDIF
      DCASY1=1.+PC(199)*SLAT
      DCASY2=1.+PC(200)*SLAT
      F11C=1.+PC(81)*DFC*SWC(1)
      IF(WW(1).NE.9898) THEN
       WC(1,11)=-(PC(65)*BP(2,2)+PC(66)*BP(4,2)+PC(67)*BP(6,2)
     $ +PC(73)*BP(8,2)+PC(74)*BP(10,2)
     $  )*SLONG*DCASY1*F11C
     $ + (PC(91)*BP(2,2)+PC(92)*BP(4,2)+PC(93)*BP(6,2)
     $ +PC(87)*BP(8,2)+PC(88)*BP(10,2)
     $  )*CLONG*DCASY1*F11C
     $  -(PC(192)*BP(4,3)+PC(194)*BP(6,3)+PC(196)*BP(8,3)
     $ +PC(198)*BP(10,3)
     $  )*S2LONG*DCASY2*F11C
     $ + (PC(191)*BP(4,3)+PC(193)*BP(6,3)+PC(195)*BP(8,3)
     $ +PC(197)*BP(10,3)
     $  )*C2LONG*DCASY2*F11C
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,11)=-(PC(91)*BT(2,2)+PC(92)*BT(4,2)+PC(93)*BT(6,2)
     $ +PC(87)*BT(8,2)+PC(88)*BT(10,2)
     $  )*SLONG*DCASY1*F11C
     $ - (PC(65)*BT(2,2)+PC(66)*BT(4,2)+PC(67)*BT(6,2)
     $ +PC(73)*BT(8,2)+PC(74)*BT(10,2)
     $  )*CLONG*DCASY1*F11C
     $  -(PC(191)*BT(4,3)+PC(193)*BT(6,3)+PC(195)*BT(8,3)
     $ +PC(197)*BT(10,3)
     $  )*S2LONG*DCASY2*F11C
     $ - (PC(192)*BT(4,3)+PC(194)*BT(6,3)+PC(196)*BT(8,3)
     $ +PC(198)*BT(10,3)
     $  )*C2LONG*DCASY2*F11C
      ENDIF
  230 CONTINUE
C       UT & MIXED UT/LONG
      UTBASY=1.
      F12B=1.+PB(82)*DFC*SWC(1)
      IF(SW(12).EQ.0) GOTO 240
      IF(WW(1).NE.9898) THEN
       WB(1,12)=(PB(69)*BT(2,1)+PB(70)*BT(4,1)+PB(71)*BT(6,1)
     $ +PB(116)*BT(8,1)+PB(117)*BT(10,1)+PB(118)*BT(12,1)
     $  )*COS(SR*(SEC-PB(72)))*UTBASY*F12B
     $ + (PB(77)*BT(4,3)+PB(78)*BT(6,3)+PB(79)*BT(8,3))
     $  *COS(SR*(SEC-PB(80))+2.*DGTR*LONG)*UTBASY*F12B*SWC(11)
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,12)=(PB(77)*BP(4,3)+PB(78)*BP(6,3)+PB(79)*BP(8,3))
     $  *COS(SR*(SEC-PB(80)+21600.)+2.*DGTR*LONG)
     $    *UTBASY*F12B*SWC(11)
      ENDIF
      UTCASY=1.
      F12C=1.+PC(82)*DFC*SWC(1)
      IF(WW(1).NE.9898) THEN
       WC(1,12)=(PC(77)*BP(3,3)+PC(78)*BP(5,3)+PC(79)*BP(7,3)
     $ +PC(165)*BP(9,3)+PC(166)*BP(11,3)+PC(167)*BP(13,3)
     $  )*COS(SR*(SEC-PC(80))+2.*DGTR*LONG)*UTCASY*F12C*SWC(11)
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,12)=-(PC(69)*BT(3,1)+PC(70)*BT(5,1)+PC(71)*BT(7,1)
     $ +PC(116)*BT(9,1)+PC(117)*BT(11,1)+PC(118)*BT(13,1)
     $  )*COS(SR*(SEC-PC(72)))*UTCASY*F12C
     $ + (PC(77)*BT(3,3)+PC(78)*BT(5,3)+PC(79)*BT(7,3)
     $ +PC(165)*BT(9,3)+PC(166)*BT(11,3)+PC(167)*BT(13,3)
     $  )*COS(SR*(SEC-PC(80)+21600.)+2.*DGTR*LONG)
     $   *UTCASY*F12C*SWC(11)
      ENDIF
  240 CONTINUE
C       MIXED LONG,UT,AP
      IF(SW(13).EQ.0) GOTO 48
      IF(SW9.EQ.-1.) GO TO 45
      IF(APD.EQ.0) GOTO 48
      IF(WW(1).NE.9898) THEN
       WB(1,13)=
     $ (PB(61)*BT(3,2)+PB(62)*BT(5,2)+PB(63)*BT(7,2))
     $  *COS(DGTR*(LONG-PB(64)))*APDF*SWC(11)+
     $  (PB(84)*BT(2,1)+PB(85)*BT(4,1)+PB(86)*BT(6,1))
     $  *COS(SR*(SEC-PB(76)))*APDF*SWC(12)
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,13)=(PB(61)*BP(3,2)+PB(62)*BP(5,2)+PB(63)*BP(7,2))
     $  *COS(DGTR*(LONG-PB(64)+90.))*APDF*SWC(11)
      ENDIF
      IF(WW(1).NE.9898) THEN 
       WC(1,13)=SWC(11)*WC(1,11)*PC(61)*APDFC
     $ +SWC(12)*WC(1,12)*PC(84)*APDFC
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,13)=SWC(11)*WC(2,11)*PC(61)*APDFC
     $ +SWC(12)*WC(2,12)*PC(84)*APDFC
      ENDIF
      GOTO 48
   45 CONTINUE
      IF(APT.EQ.0) GOTO 48
      IF(WW(1).NE.9898) THEN
       WB(1,13)=
     $  (PB(53)*BT(3,2)+PB(99)*BT(5,2)+PB(68)*BT(7,2))
     $  *COS(DGTR*(LONG-PB(98)))*APT*SWC(11)+
     $  (PB(56)*BT(2,1)+PB(57)*BT(4,1)+PB(58)*BT(6,1))
     $  *COS(SR*(SEC-PB(59)))*APT*SWC(12)
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,13)=(PB(53)*BP(3,2)+PB(99)*BP(5,2)+PB(68)*BP(7,2))
     $  *COS(DGTR*(LONG-PB(98)+90.))*APT*SWC(11)
      ENDIF
      IF(WW(1).NE.9898) THEN
       WC(1,13)=SWC(11)*WC(1,11)*PC(53)*APT
     $ +SWC(12)*WC(1,12)*PC(56)*APT
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,13)=SWC(11)*WC(2,11)*PC(53)*APT
     $ +SWC(12)*WC(2,12)*PC(56)*APT
      ENDIF
   48 CONTINUE
   49 CONTINUE
      WBT(1)=0             
      WBT(2)=0
      WCT(1)=0
      WCT(2)=0                                 
C       SUM WINDS AND CHANGE MERIDIONAL SIGN TO + NORTH
      DO 50 K=1,NSW
        WBT(1)=WBT(1)-ABS(SW(K))*WB(1,K)
        WCT(1)=WCT(1)-ABS(SW(K))*WC(1,K)
        WBT(2)=WBT(2)+ABS(SW(K))*WB(2,K)
        WCT(2)=WCT(2)+ABS(SW(K))*WC(2,K)
   50 CONTINUE
      IF(WW(1).NE.9898) WW(1)=WBT(1)*SW(24)+WCT(1)*SW(25)
      IF(WW(2).NE.9898) WW(2)=WBT(2)*SW(24)+WCT(2)*SW(25)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GLBW5M(YRD,SEC,LAT,LONG,STL,F107A,F107,AP,PB,PC,WW)
      REAL LAT,LONG
      DIMENSION WB(2,15),WC(2,15),PB(150),PC(150),WW(2)
      DIMENSION AP(1)
      COMMON/CSW/SW(25),ISW,SWC(25)
      COMMON/HWMC/WBT(2),WCT(2)
C      COMMON/VPOLY/BT(20,20),BP(20,20),CSTL,SSTL,C2STL,S2STL,
C     $ C3STL,S3STL,IYR,DAY,DF,DFA,DFC,APD,APDF,APDFC,APT,STL
      COMMON/VPOLY2/XVL,LVL,MVL,CLAT,SLAT,BT(20,20),BP(20,20)
      COMMON/LTCOMP/TLL,NSVL,CSTL,SSTL,C2STL,S2STL,C3STL,S3STL
      COMMON/LGCOMP/XLL,NGVL,CLONG,SLONG,C2LONG,S2LONG
      SAVE
      DATA DGTR/.017453/,SR/7.2722E-5/,HR/.2618/,DR/1.72142E-2/
      DATA PB14/-1./,PB18/-1./
      DATA NSW/14/,WB/30*0/,WC/30*0/
      DATA SW9/1./,LV/10/,MV/2/,NSV/2/,PSET/4./
      G0(A)=(A-4.+(PB(26)-1.)*(A-4.+(EXP(-ABS(PB(25))*(A-4.))-1.)/
     * ABS(PB(25))))
C       CONFIRM PARAMETER SET
      IF(PB(100).EQ.0) PB(100)=PSET
      IF(PB(100).NE.PSET) THEN
        WRITE(6,900) PSET,PB(100),PC(100)
  900   FORMAT(1X,'WRONG PARAMETER SET FOR GLBW5M',3F10.1)
        STOP
      ENDIF
C
      DO 10 J=1,NSW
        WB(1,J)=0
        WB(2,J)=0
        WC(1,J)=0
        WC(2,J)=0
   10 CONTINUE
      IF(SW(9).GT.0) SW9=1.
      IF(SW(9).LT.0) SW9=-1.
      IYR = YRD/1000.
      DAY = YRD - IYR*1000.
      IF(XVL.NE.LAT.OR.LV.GT.LVL.OR.MV.GT.MVL) THEN
        SLAT=SIN(DGTR*LAT)
        CLAT=COS(DGTR*LAT)
        CALL VSPHR1(SLAT,CLAT,LV,MV,BT,BP,20)
        XVL=LAT
        LVL=LV
        MVL=MV
      ENDIF
      IF(TLL.NE.STL.OR.NSV.GT.NSVL)  THEN
        SSTL = SIN(HR*STL)
        CSTL = COS(HR*STL)
        S2STL = SIN(2.*HR*STL)
        C2STL = COS(2.*HR*STL)
        TLL = STL
        NSVL=NSV
      ENDIF
      IF(DAY.NE.DAYL.OR.PB(14).NE.PB14) CD14=COS(DR*(DAY-PB(14)))
      IF(DAY.NE.DAYL.OR.PB(18).NE.PB18) CD18=COS(2.*DR*(DAY-PB(18)))
      IF(DAY.NE.DAYL.OR.PB(19).NE.PB19) CD19B=COS(2.*DR*(DAY-PB(19)))
      DAYL=DAY
      PB14=PB(14)
      PB18=PB(18)
      PB19=PB(19)
C       F10.7 EFFECT
      DF=F107-F107A
      DFA=F107A-150.
      DFC=DFA+PB(20)*DF
C       TIME INDEPENDENT
      F1B=1.
      IF(WW(1).NE.9898) THEN
       WB(1,2)=(PB(2)*BT(3,1)+PB(3)*BT(5,1)+PB(23)*BT(7,1))*F1B
      ENDIF
      WB(2,2)=0.
      F1C=1.
      WC(1,2)=0.
      IF(WW(2).NE.9898) THEN
       WC(2,2)=-(PC(2)*BT(2,1)+PC(3)*BT(4,1)+PC(23)*BT(6,1))*F1C
     $ -(PC(27)*BT(3,1)+PC(15)*BT(5,1)+PC(60)*BT(7,1))*F1C
      ENDIF
C       SYMMETRICAL ANNUAL
C       SYMMETRICAL SEMIANNUAL
      IF(WW(1).NE.9898) THEN
       WB(1,4)=(PB(17)*BT(3,1)+PB(31)*BT(5,1))*CD18
      ENDIF
      WB(2,4)=0
      WC(1,4)=0
      IF(WW(2).NE.9898) THEN
       WC(2,4)=-(PC(17)*BT(2,1)+PC(31)*BT(4,1))*CD18
      ENDIF
C       ASYMMETRICAL ANNUAL
      F5B=1.
      IF(WW(1).NE.9898) THEN
       WB(1,5)=(PB(10)*BT(2,1)+PB(11)*BT(4,1))*CD14*F5B
      ENDIF
      WB(2,5)=0
      F5C=1.
      WC(1,5)=0
      IF(WW(2).NE.9898) THEN
       WC(2,5)=-(PC(10)*BT(3,1)+PC(11)*BT(5,1))*CD14*F5C
      ENDIF
C       ASYMMETRICAL SEMIANNUAL
C       DIURNAL      
      IF(SW(7).EQ.0) GOTO 200
      F7B=1.
      F75B=1.
      IF(WW(1).NE.9898) THEN
       WB(1,7)=(PB(7)*BT(2,2)+PB(8)*BT(4,2)+PB(29)*BT(6,2)
     $         +PB(89)*BT(3,2)
     $  )*SSTL*F7B
     $ +(PB(13)*BT(3,2)+PB(146)*BT(5,2))
     $    *CD14*SSTL*F75B*SWC(5)
     $ + (PB(4)*BT(2,2)+PB(5)*BT(4,2)+PB(28)*BT(6,2)
     $         +PB(88)*BT(3,2)
     $  )*CSTL*F7B
     $ +(PB(12)*BT(3,2)+PB(145)*BT(5,2))
     $      *CD14*CSTL*F75B*SWC(5)
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,7)=-(PB(4)*BP(2,2)+PB(5)*BP(4,2)+PB(28)*BP(6,2)
     $         +PB(88)*BP(3,2)
     $  )*SSTL*F7B
     $ -(PB(12)*BP(3,2)+PB(145)*BP(5,2))
     $    *CD14*SSTL*F75B*SWC(5)
     $ + (PB(7)*BP(2,2)+PB(8)*BP(4,2)+PB(29)*BP(6,2)
     $         +PB(89)*BP(3,2)
     $  )*CSTL*F7B
     $ +(PB(13)*BP(3,2)+PB(146)*BP(5,2))
     $    *CD14*CSTL*F75B*SWC(5)
      ENDIF
      F7C=1.
      F75C=1.
      IF(WW(1).NE.9898) THEN
       WC(1,7)=-(PC(4)*BP(3,2)+PC(5)*BP(5,2)+PC(28)*BP(7,2)
     $         +PC(88)*BP(2,2)
     $   +PC(141)*BP(9,2)+PC(143)*BP(11,2)
     $  )*SSTL*F7C
     $ -(PC(12)*BP(2,2)+PC(145)*BP(4,2))
     $    *CD14*SSTL
     $   *F75C*SWC(5)
     $ +(PC(7)*BP(3,2)+PC(8)*BP(5,2)+PC(29)*BP(7,2)
     $         +PC(89)*BP(2,2)
     $ +PC(142)*BP(9,2)+PC(144)*BP(11,2)
     $  )*CSTL*F7C
     $ +(PC(13)*BP(2,2)+PC(146)*BP(4,2))
     $     *CD14*CSTL
     $   *F75C*SWC(5)
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,7)=-(PC(7)*BT(3,2)+PC(8)*BT(5,2)+PC(29)*BT(7,2)
     $         +PC(89)*BT(2,2)
     $ +PC(142)*BT(9,2)+PC(144)*BT(11,2)
     $  )*SSTL*F7C
     $ -(PC(13)*BT(2,2)+PC(146)*BT(4,2))
     $    *CD14*SSTL
     $   *F75C*SWC(5)
     $ -(PC(4)*BT(3,2)+PC(5)*BT(5,2)+PC(28)*BT(7,2)
     $         +PC(88)*BT(2,2)
     $ +PC(141)*BT(9,2)+PC(143)*BT(11,2)
     $  )*CSTL*F7C
     $ -(PC(12)*BT(2,2)+PC(145)*BT(4,2))
     $    *CD14*CSTL
     $   *F75C*SWC(5)
      ENDIF
  200 CONTINUE
C       SEMIDIURNAL
      IF(SW(8).EQ.0) GOTO 210
      F8B=1.+PB(90)*DFC*SWC(1)
      IF(WW(1).NE.9898) THEN
       WB(1,8)=(PB(9)*BT(3,3)+PB(43)*BT(5,3)+PB(111)*BT(7,3)
     $         +PB(98)*BT(4,3)
     $   +(PB(34)*BT(4,3)+PB(148)*BT(6,3))*CD14*SWC(5)
     $   +(PB(37)*BT(4,3))*CD19B*SWC(6)
     $  )*S2STL*F8B
     $ +(PB(6)*BT(3,3)+PB(42)*BT(5,3)+PB(110)*BT(7,3)
     $         +PB(96)*BT(4,3)
     $   +(PB(24)*BT(4,3)+PB(147)*BT(6,3))*CD14*SWC(5)
     $   +(PB(36)*BT(4,3))*CD19B*SWC(6)
     $  )*C2STL*F8B
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,8)=-(PB(6)*BP(3,3)+PB(42)*BP(5,3)+PB(110)*BP(7,3)
     $          +PB(96)*BP(4,3)
     $   +(PB(24)*BP(4,3)+PB(147)*BP(6,3))*CD14*SWC(5)
     $   +(PB(36)*BP(4,3))*CD19B*SWC(6)
     $  )*S2STL*F8B
     $   + (PB(9)*BP(3,3)+PB(43)*BP(5,3)+PB(111)*BP(7,3)
     $          +PB(98)*BP(4,3)
     $   +(PB(34)*BP(4,3)+PB(148)*BP(6,3))*CD14*SWC(5)
     $   +(PB(37)*BP(4,3))*CD19B*SWC(6)
     $  )*C2STL*F8B
      ENDIF
      F8C=1.+PC(90)*DFC*SWC(1)
      IF(WW(1).NE.9898) THEN
       WC(1,8)=-(PC(6)*BP(4,3)+PC(42)*BP(6,3)+PC(110)*BP(8,3)
     $          +PC(96)*BP(3,3)
     $   +(PC(24)*BP(3,3)+PC(147)*BP(5,3))*CD14*SWC(5)
     $   +(PC(36)*BP(3,3))*CD19B*SWC(6)
     $  )*S2STL*F8C
     $ +(PC(9)*BP(4,3)+PC(43)*BP(6,3)+PC(111)*BP(8,3)
     $          +PC(98)*BP(3,3)
     $   +(PC(34)*BP(3,3)+PC(148)*BP(5,3))*CD14*SWC(5)
     $   +(PC(37)*BP(3,3))*CD19B*SWC(6)
     $  )*C2STL*F8C
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,8)=-(PC(9)*BT(4,3)+PC(43)*BT(6,3)+PC(111)*BT(8,3)
     $          +PC(98)*BT(3,3)
     $   +(PC(34)*BT(3,3)+PC(148)*BT(5,3))*CD14*SWC(5)
     $   +(PC(37)*BT(3,3))*CD19B*SWC(6)
     $  )*S2STL*F8C
     $ - (PC(6)*BT(4,3)+PC(42)*BT(6,3)
     $          +PC(96)*BT(3,3)
     $   +PC(110)*BT(8,3)
     $   +(PC(24)*BT(3,3)+PC(147)*BT(5,3))*CD14*SWC(5)
     $   +(PC(36)*BT(3,3))*CD19B*SWC(6)
     $  )*C2STL*F8C
      ENDIF
  210 CONTINUE
C        TERDIURNAL
C        MAGNETIC ACTIVITY
      IF(SW(9).EQ.0) GOTO 40
      IF(SW9.EQ.-1.) GOTO 30
C           daily AP
      APD=AP(1)-4.
      APDF=(APD+(PB(45)-1.)*(APD+(EXP(-PB(44)*APD)-1.)/PB(44)))
C      APDFC=(APD+(PC(45)-1.)*(APD+(EXP(-PC(44)*APD)-1.)/PC(44)))
      APDFC=APDF
      IF(APD.EQ.0) GOTO 40
      IF(WW(1).NE.9898) THEN
       WB(1,9)=(PB(46)*BT(3,1)+PB(35)*BT(5,1))*APDF
     $    +(PB(122)*BT(2,2)+PB(123)*BT(4,2)+PB(124)*BT(6,2)
     $       )*COS(HR*(STL-PB(125)))*APDF*SWC(7)
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,9)=
     $   (PB(122)*BP(2,2)+PB(123)*BP(4,2)+PB(124)*BP(6,2)
     $     )*COS(HR*(STL-PB(125)+6.))*APDF*SWC(7)
      ENDIF
      IF(WW(1).NE.9898) THEN
       WC(1,9)=
     $   (PC(122)*BP(3,2)+PC(123)*BP(5,2)+PC(124)*BP(7,2)
     $       )*COS(HR*(STL-PC(125)))*APDFC*SWC(7)
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,9)=-(PC(46)*BT(2,1)+PC(35)*BT(4,1))*APDFC
     $  +(PC(122)*BT(3,2)+PC(123)*BT(5,2)+PC(124)*BT(7,2)
     $       )*COS(HR*(STL-PC(125)+6.))*APDFC*SWC(7)
      ENDIF
      GO TO 40
   30 CONTINUE
      IF(PB(25).LT.1.E-4) PB(25)=1.E-4
c Problem not solved because AP is of dimension 1
c     APT=G0(AP(2))
       APT=0
      IF(APT.EQ.0) GOTO 40
      IF(WW(1).NE.9898) THEN
       WB(1,9)=(PB(97)*BT(3,1)+PB(55)*BT(5,1))*APT
     $    +(PB(129)*BT(2,2)+PB(130)*BT(4,2)+PB(131)*BT(6,2)
     $       )*COS(HR*(STL-PB(132)))*APT*SWC(7)
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,9)=
     $   (PB(129)*BP(2,2)+PB(130)*BP(4,2)+PB(131)*BP(6,2)
     $     )*COS(HR*(STL-PB(132)+6.))*APT*SWC(7)
      ENDIF
      IF(WW(1).NE.9898) THEN
       WC(1,9)=
     $   (PC(129)*BP(3,2)+PC(130)*BP(5,2)+PC(131)*BP(7,2)
     $       )*COS(HR*(STL-PC(132)))*APT*SWC(7)
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,9)=-(PC(97)*BT(2,1)+PC(55)*BT(4,1))*APT
     $  +(PC(129)*BT(3,2)+PC(130)*BT(5,2)+PC(131)*BT(7,2)
     $       )*COS(HR*(STL-PC(132)+6.))*APT*SWC(7)
      ENDIF
  40  CONTINUE
      WBT(1)=0
      WBT(2)=0
      WCT(1)=0
      WCT(2)=0
C       SUM WINDS AND CHANGE MERIDIONAL SIGN TO + NORTH
      DO 50 K=1,NSW
        WBT(1)=WBT(1)-ABS(SW(K))*WB(1,K)
        WCT(1)=WCT(1)-ABS(SW(K))*WC(1,K)
        WBT(2)=WBT(2)+ABS(SW(K))*WB(2,K)
        WCT(2)=WCT(2)+ABS(SW(K))*WC(2,K)
   50 CONTINUE
      IF(WW(1).NE.9898) WW(1)=WBT(1)*SW(24)+WCT(1)*SW(25)
      IF(WW(2).NE.9898) WW(2)=WBT(2)*SW(24)+WCT(2)*SW(25)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GLBW5S(IYD,LAT,LONG,STL,PB,PC,WW)
      REAL LAT,LONG
      DIMENSION WB(2,15),WC(2,15),PB(100),PC(100),WW(2)
      COMMON/CSW/SW(25),ISW,SWC(25)
      COMMON/HWMC/WBT(2),WCT(2)
      COMMON/VPOLY2/XVL,LVL,MVL,CLAT,SLAT,BT(20,20),BP(20,20)
      COMMON/LTCOMP/TLL,NSVL,CSTL,SSTL,C2STL,S2STL,C3STL,S3STL
      COMMON/LGCOMP/XLL,NGVL,CLONG,SLONG,C2LONG,S2LONG
      SAVE
      DATA DGTR/.017453/,SR/7.2722E-5/,HR/.2618/,DR/1.72142E-2/
      DATA PB14/-1./,PB18/-1./,PC14/-1./,PC18/-1./,PSET/5./
      DATA NSW/14/,WB/30*0/,WC/30*0/
C       CONFIRM PARAMETER SET
      IF(PB(100).EQ.0) PB(100)=PSET
      IF(PB(100).NE.PSET) THEN
        WRITE(6,900) PSET,PB(100),PC(100)
  900   FORMAT(1X,'WRONG PARAMETER SET FOR GLBW5S',3F10.1)
        STOP
      ENDIF
C
      DO 10 J=1,NSW
        WB(1,J)=0
        WB(2,J)=0
        WC(1,J)=0
        WC(2,J)=0
   10 CONTINUE
      IYR = IYD/1000
      DAY = IYD - IYR*1000
C
      LV=7
      MV=2
      IF(XVL.NE.LAT.OR.LV.GT.LVL.OR.MV.GT.MVL) THEN
        SLAT=SIN(DGTR*LAT)
        CLAT=COS(DGTR*LAT)
        CALL VSPHR1(SLAT,CLAT,LV,MV,BT,BP,20)
        PLG10=SLAT
        PLG30=(5.*SLAT*SLAT-3.)*SLAT/2.
        XVL=LAT
        LVL=LV
        MVL=MV
      ENDIF
C
      NSV=2
      IF(TLL.NE.STL.OR.NSV.GT.NSVL)  THEN
        SSTL = SIN(HR*STL)
        CSTL = COS(HR*STL)
        S2STL = SIN(2.*HR*STL)
        C2STL = COS(2.*HR*STL)
        TLL = STL
        NSVL=NSV
      ENDIF
      IF(DAY.NE.DAYL.OR.PB(14).NE.PB14) CD14B=COS(DR*(DAY-PB(14)))
      IF(DAY.NE.DAYL.OR.PC(14).NE.PC14) CD14C=COS(DR*(DAY-PC(14)))
      IF(DAY.NE.DAYL.OR.PB(18).NE.PB18) CD18B=COS(2.*DR*(DAY-PB(18)))
      IF(DAY.NE.DAYL.OR.PC(18).NE.PC18) CD18C=COS(2.*DR*(DAY-PC(18)))
      IF(DAY.NE.DAYL.OR.PB(19).NE.PB19) CD19B=COS(2.*DR*(DAY-PB(19)))
      IF(DAY.NE.DAYL.OR.PB(25).NE.PB25) CD25B=COS(DR*(DAY-PB(25)))
C      IF(DAY.NE.DAYL.OR.PC(25).NE.PC25) CD25C=COS(DR*(DAY-PC(25)))
      IF(DAY.NE.DAYL.OR.PB(26).NE.PB26) CD26B=COS(DR*(DAY-PB(26)))
C      IF(DAY.NE.DAYL.OR.PC(26).NE.PC26) CD26C=COS(DR*(DAY-PC(26)))
      IF(DAY.NE.DAYL.OR.PC(32).NE.PC32) CD32C=COS(DR*(DAY-PC(32)))
      IF(DAY.NE.DAYL.OR.PC(39).NE.PC39) CD39C=COS(2.*DR*(DAY-PC(39)))
      IF(DAY.NE.DAYL.OR.PC(64).NE.PC64) CD64C=COS(DR*(DAY-PC(64)))
      IF(DAY.NE.DAYL.OR.PC(87).NE.PC87) CD87C=COS(2.*DR*(DAY-PC(87)))
      DAYL=DAY           
      PB14=PB(14)
      PC14=PC(14)
      PB18=PB(18)
      PC18=PC(18)
      PB19=PB(19)
      PB25=PB(25)
      PC25=PC(25)
      PB26=PB(26)
      PC26=PC(26)
      PC32=PC(32)
      PC39=PC(39)
      PC64=PC(64)
      PC87=PC(87)
C
      NGV=1
      IF(XLL.NE.LONG.OR.NGV.GT.NGVL) THEN
        SLONG=SIN(DGTR*LONG)
        CLONG=COS(DGTR*LONG)
        XLL=LONG
        NGVL=NGV
      ENDIF
C       TIME INDEPENDENT
      F1B=1.
      IF(WW(1).NE.9898) THEN
       WB(1,2)=(PB(2)*BT(3,1)+PB(3)*BT(5,1)+PB(23)*BT(7,1))*F1B
      ENDIF
      WB(2,2)=0.
      F1C=1.
      WC(1,2)=0.
      IF(WW(2).NE.9898) THEN
       WC(2,2)=-(PC(2)*BT(2,1)+PC(3)*BT(4,1)+PC(23)*BT(6,1))*F1C
     $ -(PC(27)*BT(3,1)+PC(15)*BT(5,1)+PC(60)*BT(7,1))*F1C
      ENDIF
C       SYMMETRICAL ANNUAL
      IF(WW(2).NE.9898) THEN
       WC(2,3)=-(PC(48)*BT(2,1)+PC(30)*BT(4,1))*CD32C
      ENDIF
C       SYMMETRICAL SEMIANNUAL
      IF(WW(1).NE.9898) THEN
       WB(1,4)=(PB(17)*BT(3,1)+PB(31)*BT(5,1))*CD18B
      ENDIF
      WB(2,4)=0
      WC(1,4)=0
      IF(WW(2).NE.9898) THEN
       WC(2,4)=-(PC(17)*BT(2,1)+PC(31)*BT(4,1)+PC(50)*BT(6,1))*CD18C
      ENDIF
C       ASYMMETRICAL ANNUAL
      F5B=1.
      IF(WW(1).NE.9898) THEN
       WB(1,5)=(PB(10)*BT(2,1)+PB(11)*BT(4,1))*CD14B*F5B
      ENDIF
      WB(2,5)=0
      F5C=1.
      WC(1,5)=0
      IF(WW(2).NE.9898) THEN
       WC(2,5)=-(PC(10)*BT(3,1)+PC(11)*BT(5,1)+PC(21)*BT(7,1))*CD14C*F5C
      ENDIF
C       ASYMMETRICAL SEMIANNUAL
      IF(WW(2).NE.9898) THEN
       WC(2,6)=-(PC(38)*BT(3,1)+PC(99)*BT(5,1))*CD39C
      ENDIF
C       DIURNAL      
      IF(SW(7).EQ.0) GOTO 200
      F7B=1.
      F75B=1.
      IF(WW(1).NE.9898) THEN
       WB(1,7)=(PB(7)*BT(2,2)+PB(8)*BT(4,2)+PB(29)*BT(6,2)
     $         +PB(89)*BT(3,2)
     $  )*SSTL*F7B
     $ +(PB(13)*BT(3,2))
     $    *CD25B*SSTL*F75B*SWC(5)
     $ + (PB(4)*BT(2,2)+PB(5)*BT(4,2)+PB(28)*BT(6,2)
     $         +PB(88)*BT(3,2)
     $  )*CSTL*F7B
     $ +(PB(12)*BT(3,2))
     $      *CD25B*CSTL*F75B*SWC(5)
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,7)=-(PB(4)*BP(2,2)+PB(5)*BP(4,2)+PB(28)*BP(6,2)
     $         +PB(88)*BP(3,2)
     $  )*SSTL*F7B
     $ -(PB(12)*BP(3,2))
     $    *CD25B*SSTL*F75B*SWC(5)
     $ + (PB(7)*BP(2,2)+PB(8)*BP(4,2)+PB(29)*BP(6,2)
     $         +PB(89)*BP(3,2)
     $  )*CSTL*F7B
     $ +(PB(13)*BP(3,2))
     $    *CD25B*CSTL*F75B*SWC(5)
      ENDIF
      F7C=1.
      F75C=1.
      IF(WW(1).NE.9898) THEN
       WC(1,7)=-(PC(4)*BP(3,2)+PC(5)*BP(5,2)+PC(28)*BP(7,2)
     $         +PC(88)*BP(2,2)
     $  )*SSTL*F7C
     $ -(PC(12)*BP(2,2))
     $    *CD25B*SSTL
     $   *F75C*SWC(5)
     $ +(PC(7)*BP(3,2)+PC(8)*BP(5,2)+PC(29)*BP(7,2)
     $         +PC(89)*BP(2,2)
     $  )*CSTL*F7C
     $ +(PC(13)*BP(2,2))
     $     *CD25B*CSTL
     $   *F75C*SWC(5)
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,7)=-(PC(7)*BT(3,2)+PC(8)*BT(5,2)+PC(29)*BT(7,2)
     $         +PC(89)*BT(2,2)
     $  )*SSTL*F7C
     $ -(PC(13)*BT(2,2))
     $    *CD25B*SSTL
     $   *F75C*SWC(5)
     $ -(PC(4)*BT(3,2)+PC(5)*BT(5,2)+PC(28)*BT(7,2)
     $         +PC(88)*BT(2,2)
     $  )*CSTL*F7C
     $ -(PC(12)*BT(2,2))
     $    *CD25B*CSTL
     $   *F75C*SWC(5)
      ENDIF
  200 CONTINUE
C       SEMIDIURNAL
      IF(SW(8).EQ.0) GOTO 210
      F8B=1.
      IF(WW(1).NE.9898) THEN
       WB(1,8)=(PB(9)*BT(3,3)+PB(43)*BT(5,3)+PB(35)*BT(7,3)
     $         +PB(98)*BT(4,3)
     $   +(PB(34)*BT(4,3))*CD26B*SWC(5)
     $   +(PB(37)*BT(4,3))*CD19B*SWC(6)
     $  )*S2STL*F8B
     $ +(PB(6)*BT(3,3)+PB(42)*BT(5,3)+PB(33)*BT(7,3)
     $         +PB(96)*BT(4,3)
     $   +(PB(24)*BT(4,3))*CD26B*SWC(5)
     $   +(PB(36)*BT(4,3))*CD19B*SWC(6)
     $  )*C2STL*F8B
      ENDIF
      IF(WW(2).NE.9898) THEN
       WB(2,8)=-(PB(6)*BP(3,3)+PB(42)*BP(5,3)+PB(33)*BP(7,3)
     $          +PB(96)*BP(4,3)
     $   +(PB(24)*BP(4,3))*CD26B*SWC(5)
     $   +(PB(36)*BP(4,3))*CD19B*SWC(6)
     $  )*S2STL*F8B
     $   + (PB(9)*BP(3,3)+PB(43)*BP(5,3)+PB(35)*BP(7,3)
     $          +PB(98)*BP(4,3)
     $   +(PB(34)*BP(4,3))*CD26B*SWC(5)
     $   +(PB(37)*BP(4,3))*CD19B*SWC(6)
     $  )*C2STL*F8B
      ENDIF
      F8C=1.
      IF(WW(1).NE.9898) THEN
       WC(1,8)=-(PC(6)*BP(4,3)+PC(42)*BP(6,3)+PC(33)*BP(8,3)
     $          +PC(96)*BP(3,3)
     $   +(PC(24)*BP(3,3))*CD26B*SWC(5)
     $   +(PC(36)*BP(3,3))*CD19B*SWC(6)
     $  )*S2STL*F8C
     $ +(PC(9)*BP(4,3)+PC(43)*BP(6,3)+PC(35)*BP(8,3)
     $          +PC(98)*BP(3,3)
     $   +(PC(34)*BP(3,3))*CD26B*SWC(5)
     $   +(PC(37)*BP(3,3))*CD19B*SWC(6)
     $  )*C2STL*F8C
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,8)=-(PC(9)*BT(4,3)+PC(43)*BT(6,3)+PC(35)*BT(8,3)
     $          +PC(98)*BT(3,3)
     $   +(PC(34)*BT(3,3))*CD26B*SWC(5)
     $   +(PC(37)*BT(3,3))*CD19B*SWC(6)
     $  )*S2STL*F8C
     $ - (PC(6)*BT(4,3)+PC(42)*BT(6,3)+PC(33)*BT(8,3)
     $          +PC(96)*BT(3,3)
     $   +(PC(24)*BT(3,3))*CD26B*SWC(5)
     $   +(PC(36)*BT(3,3))*CD19B*SWC(6)
     $  )*C2STL*F8C
      ENDIF
  210 CONTINUE
C        LONGITUDINAL
      IF(SW(10).EQ.0.OR.SW(11).EQ.0) GOTO 230
      IF(WW(1).NE.9898) THEN
       WC(1,11)=
     $ - (PC(65)*BP(2,2)+PC(66)*BP(4,2)+PC(67)*BP(6,2)
     $   +PC(75)*BP(3,2)+PC(76)*BP(5,2)+ PC(77)*BP(7,2)
     $   +(PC(57)*BP(2,2)+PC(59)*BP(4,2)+PC(62)*BP(6,2)
     $    +PC(51)*BP(3,2)+PC(53)*BP(5,2)+PC(55)*BP(7,2))
     $     *CD64C*SWC(3)
     $   +(PC(74)*BP(2,2)+PC(82)*BP(4,2)+PC(85)*BP(6,2)
     $    +PC(68)*BP(3,2)+PC(70)*BP(5,2)+PC(72)*BP(7,2))
     $     *CD87C*SWC(4)
     $  )*SLONG
     $ + (PC(91)*BP(2,2)+PC(92)*BP(4,2)+PC(93)*BP(6,2)
     $   +PC(78)*BP(3,2)+PC(79)*BP(5,2)+PC(80)*BP(7,2)
     $   +(PC(58)*BP(2,2)+PC(61)*BP(4,2)+PC(63)*BP(6,2)
     $    +PC(52)*BP(3,2)+PC(54)*BP(5,2)+PC(56)*BP(7,2))
     $     *CD64C*SWC(3)
     $   +(PC(81)*BP(2,2)+PC(84)*BP(4,2)+PC(86)*BP(6,2)
     $    +PC(69)*BP(3,2)+PC(71)*BP(5,2)+PC(73)*BP(7,2))
     $     *CD87C*SWC(4)
     $  )*CLONG
      ENDIF
      IF(WW(2).NE.9898) THEN
       WC(2,11)=
     $ - (PC(91)*BT(2,2)+PC(92)*BT(4,2)+PC(93)*BT(6,2)
     $   +PC(78)*BT(3,2)+PC(79)*BT(5,2)+PC(80)*BT(7,2)
     $   +(PC(58)*BT(2,2)+PC(61)*BT(4,2)+PC(63)*BT(6,2)
     $    +PC(52)*BT(3,2)+PC(54)*BT(5,2)+PC(56)*BT(7,2))
     $     *CD64C*SWC(3)
     $   +(PC(81)*BT(2,2)+PC(84)*BT(4,2)+PC(86)*BT(6,2)
     $    +PC(69)*BT(3,2)+PC(71)*BT(5,2)+PC(73)*BT(7,2))
     $     *CD87C*SWC(4)
     $  )*SLONG
     $ - (PC(65)*BT(2,2)+PC(66)*BT(4,2)+PC(67)*BT(6,2)
     $   +PC(75)*BT(3,2)+PC(76)*BT(5,2)+PC(77)*BT(7,2)
     $   +(PC(57)*BT(2,2)+PC(59)*BT(4,2)+PC(62)*BT(6,2)
     $    +PC(51)*BT(3,2)+PC(53)*BT(5,2)+PC(55)*BT(7,2))
     $     *CD64C*SWC(3)
     $   +(PC(74)*BT(2,2)+PC(82)*BT(4,2)+PC(85)*BT(6,2)
     $    +PC(68)*BT(3,2)+PC(70)*BT(5,2)+PC(72)*BT(7,2))
     $     *CD87C*SWC(4)
     $  )*CLONG
      ENDIF
  230 CONTINUE
      WBT(1)=0
      WBT(2)=0
      WCT(1)=0
      WCT(2)=0
C       SUM WINDS AND CHANGE MERIDIONAL SIGN TO + NORTH
      DO 50 K=1,NSW
        WBT(1)=WBT(1)-ABS(SW(K))*WB(1,K)
        WCT(1)=WCT(1)-ABS(SW(K))*WC(1,K)
        WBT(2)=WBT(2)+ABS(SW(K))*WB(2,K)
        WCT(2)=WCT(2)+ABS(SW(K))*WC(2,K)
   50 CONTINUE
      IF(WW(1).NE.9898) WW(1)=WBT(1)*SW(24)+WCT(1)*SW(25)
      IF(WW(2).NE.9898) WW(2)=WBT(2)*SW(24)+WCT(2)*SW(25)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE TSELEC(SV)
C        SET SWITCHES
C        SW FOR MAIN TERMS, SWC FOR CROSS TERMS
      DIMENSION SV(1),SAV(25),SVV(1)
      COMMON/CSW/SW(25),ISW,SWC(25)
      DO 100 I = 1,25
        SAV(I)=SV(I)
        SW(I)=AMOD(SV(I),2.)
        IF(ABS(SV(I)).EQ.1.OR.ABS(SV(I)).EQ.2.) THEN
          SWC(I)=1.
        ELSE
          SWC(I)=0.
        ENDIF
  100 CONTINUE
      ISW=64999
      RETURN
      ENTRY TRETRV(SVV)
      DO 200 I=1,25
        SVV(I)=SAV(I)
  200 CONTINUE
      END
C-----------------------------------------------------------------------
      SUBROUTINE VSPHR1(C,S,L,M,BT,BP,LMAX)
C      CALCULATE VECTOR SPHERICAL HARMONIC B FIELD THETA AND PHI
C      FUNCTIONS BT,BP THROUGH ORDER L,M FOR COLATITUDE (THETA)
C      WITH COSINE C AND SINE S OF COLATITUDE
C      BT(L+1,M+1)= [(L-M+1) P(L+1,M) - (L+1) P(L,M) COS(THETA)] /
C                [SQRT(L(L+1)) SIN(THETA)]
C      BP(L+1,M+1)= M P(L,M) /[SQRT(L(L+1)) SIN(THETA)]
C       RESULT FOR GIVEN L,M SAVED IN BT AND BP AT ONE HIGHER INDEX NUM
      DIMENSION BT(LMAX,1),BP(LMAX,1),PLG(20,20)
      SAVE
      DATA DGTR/1.74533E-2/
      IF(M.GT.L.OR.L.GT.LMAX-1) THEN
        WRITE(6,100) L,M,LMAX
  100   FORMAT('ILLEGAL INDICIES TO VSPHER',3I6)
        RETURN
      ENDIF
      BT(1,1)=0
      BP(1,1)=0
      IF(L.EQ.0.AND.M.EQ.0) RETURN
      CALL LEGPL1(C,S,L+1,M,PLG,20)
      IF(ABS(S).LT.1.E-5) THEN
        IC=SIGN(1.,S)
        S=0
      ENDIF
      DO 20 LL=1,L
        SQT=SQRT(FLOAT(LL)*(FLOAT(LL)+1))
        LMX=MIN(LL,M)
        DO 15 MM=0,LMX
          IF(S.EQ.0) THEN
            IF(MM.NE.1) THEN
              BT(LL+1,MM+1)=0
              BP(LL+1,MM+1)=0
            ELSE
              BT(LL+1,MM+1)=(LL*(LL+1)*(LL+2)*.5*(IC)**(LL+2)
     $           -(LL+1)*C*LL*(LL+1)*.5*(IC)**(LL+1))/SQT
              BP(LL+1,MM+1)=MM*LL*(LL+1)*.5*(IC)**(LL+1)/SQT
            ENDIF
          ELSE
            BT(LL+1,MM+1)=((LL-MM+1)*PLG(LL+2,MM+1)
     $      -(LL+1)*C*PLG(LL+1,MM+1))/(S*SQT)
            BP(LL+1,MM+1)=MM*PLG(LL+1,MM+1)/(S*SQT)
          ENDIF
   15   CONTINUE
   20 CONTINUE
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGPL1(C,S,L,M,PLG,LMAX)
C      CALCULATE LEGENDRE POLYNOMIALS PLG(L+1,M+1) THROUGH ORDER L,M 
C      FOR COSINE C AND SINE S OF COLATITUDE
      DIMENSION PLG(LMAX,1)
      SAVE
      DATA DGTR/1.74533E-2/
      IF(M.GT.L.OR.L.GT.LMAX-1) THEN
        WRITE(6,99) L,M,LMAX
   99 FORMAT(1X,'ILLEGAL INDICIES TO LEGPOL',3I5)
        RETURN
      ENDIF
      PLG(1,1)=1.
      IF(L.EQ.0.AND.M.EQ.0) RETURN
C      CALCULATE L=M CASE AND L=M+1
      DO 10 MM=0,M
        IF(MM.GT.0) PLG(MM+1,MM+1)=PLG(MM,MM)*(2.*MM-1.)*S
        IF(L.GT.MM) PLG(MM+2,MM+1)=PLG(MM+1,MM+1)*(2.*MM+1)*C
   10 CONTINUE
      IF(L.EQ.1) RETURN
      MMX=MIN(M,L-2)
      DO 30 MM=0,MMX
        DO 20 LL=MM+2,L
          PLG(LL+1,MM+1)=((2.*LL-1.)*C*PLG(LL,MM+1)-
     $     (LL+MM-1.)*PLG(LL-1,MM+1))/(LL-MM)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
C        CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
C        X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        N: SIZE OF ARRAYS X,Y
C        YP1,YPN: SPECIFIED DERIVATIVES AT X(1) AND X(N); VALUES
C                 >= 1E30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
C        Y2: OUTPUT ARRAY OF SECOND DERIVATIVES
      PARAMETER (NMAX=100)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      SAVE
      IF(YP1.GT..99E30) THEN
        Y2(1)=0
        U(1)=0
      ELSE
        Y2(1)=-.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     $    /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
   11 CONTINUE
      IF(YPN.GT..99E30) THEN
        QN=0
        UN=0
      ELSE
        QN=.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
   12 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
C        CALCULATE CUBIC SPLINE INTERP VALUE
C        XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        Y2A: ARRAY OF SECOND DERIVATIVES
C        N: SIZE OF ARRAYS XA,YA,Y2A
C        X: ABSCISSA FOR INTERPOLATION
C        Y: OUTPUT VALUE
      DIMENSION XA(N),YA(N),Y2A(N)
      SAVE
      KLO=1
      KHI=N
    1 CONTINUE
      IF(KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X) THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
        GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF(H.EQ.0) WRITE(6,*) 'BAD XA INPUT TO SPLINT'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     $  ((A*A*A-A)*Y2A(KLO)+(B*B*B-B)*Y2A(KHI))*H*H/6.
      RETURN
      END
C-----------------------------------------------------------------------
      BLOCK DATA INITW5
C       For wind model GWS
      COMMON/CSW/SW(25),ISW,SWC(25)
      COMMON/VPOLY2/XVL,LVL,MVL,CLAT,SLAT,BT(20,20),BP(20,20)
      COMMON/LTCOMP/TLL,NSVL,CSTL,SSTL,C2STL,S2STL,C3STL,S3STL
      COMMON/LGCOMP/XLL,NGVL,CLONG,SLONG,C2LONG,S2LONG
      DATA ISW/0/
      DATA XVL/-999./,LVL/-1/,MVL/-1/
      DATA TLL/-999./,NSVL/-1/
      DATA XLL/-999./,NGVL/-1/
      END
C-----------------------------------------------------------------------
      BLOCK DATA GWSBK5
C          HWM93    28-JAN-93   
      COMMON/PARMW5/PBA1(50),PBA2(50),PBA3(50),PBA4(50),
     $PCA1(50),PCA2(50),PCA3(50),PCA4(50),
     $PBB1(50),PBB2(50),PBB3(50),PCB1(50),PCB2(50),PCB3(50),
     $PBC1(50),PBC2(50),PBC3(50),PCC1(50),PCC2(50),PCC3(50),
     $PBD1(50),PBD2(50),PBD3(50),PCD1(50),PCD2(50),PCD3(50),
     $PBE1(50),PBE2(50),PBE3(50),PCE1(50),PCE2(50),PCE3(50),
     $PBF1(50),PBF2(50),PBF3(50),PCF1(50),PCF2(50),PCF3(50),
     $PBG1(50),PBG2(50),PBG3(50),PCG1(50),PCG2(50),PCG3(50),
     $PBH1(50),PBH2(50),PBH3(50),PCH1(50),PCH2(50),PCH3(50),
     $PBI1(50),PBI2(50),PCI1(50),PCI2(50),PBJ1(50),PBJ2(50),
     $PCJ1(50),PCJ2(50),PBK1(50),PBK2(50),PCK1(50),PCK2(50),
     $PBL1(50),PBL2(50),PCL1(50),PCL2(50),PBM1(50),PBM2(50),
     $PCM1(50),PCM2(50),PBN1(50),PBN2(50),PCN1(50),PCN2(50),
     $PBO1(50),PBO2(50),PCO1(50),PCO2(50),PBP1(50),PBP2(50),
     $PCP1(50),PCP2(50),PBQ1(50),PBQ2(50),PCQ1(50),PCQ2(50),
     $PBR1(50),PBR2(50),PCR1(50),PCR2(50),PBS1(50),PBS2(50),
     $PCS1(50),PCS2(50),PBT1(50),PBT2(50),PCT1(50),PCT2(50),
     $PBU1(50),PBU2(50),PCU1(50),PCU2(50)
      COMMON/DATW/ISDATE(3),ISTIME(2),NAME(2)
      character*4 ISDATE,ISTIME,NAME
      DATA ISDATE/'28-J','AN-9','3   '/,ISTIME/'20:3','5:39'/
      DATA NAME/'HWM9','3   '/
C         WINF
      DATA PBA1/
     *  0.00000E+00,-1.31640E+01,-1.52352E+01, 1.00718E+02, 3.94962E+00,
     *  2.19452E-01, 8.03296E+01,-1.02032E+00,-2.02149E-01, 5.67263E+01,
     *  0.00000E+00,-6.05459E+00, 6.68106E+00,-8.49486E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 8.39399E+01, 0.00000E+00, 9.96285E-02,
     *  0.00000E+00,-2.66243E-02, 0.00000E+00,-1.32373E+00, 1.39396E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 3.36523E+01,-7.42795E-01,-3.89352E+00,-7.81354E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 3.76631E+00,-1.22024E+00,
     * -5.47580E-01, 1.09146E+00, 9.06245E-01, 2.21119E-02, 0.00000E+00,
     *  7.73919E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBA2/
     * -3.82415E-01, 0.00000E+00, 1.76202E-01, 0.00000E+00,-6.77651E-01,
     *  1.10357E+00, 2.25732E+00, 0.00000E+00, 1.54237E+04, 0.00000E+00,
     *  1.27411E-01,-2.84314E-03, 4.62562E-01,-5.34596E+01,-7.23808E+00,
     *  0.00000E+00, 0.00000E+00, 4.52770E-01,-8.50922E+00,-2.85389E-01,
     *  2.12000E+01, 6.80171E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -2.72552E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.64109E+03,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-1.47320E+00,-2.98179E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.05412E-02,
     *  4.93452E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 7.98332E-02,-5.30954E+01, 2.10211E-02, 3.00000E+00/
      DATA PBA3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.79843E-01,
     *  1.81152E-01, 0.00000E+00, 0.00000E+00,-6.24673E-02,-5.37589E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-8.94418E-02, 3.70413E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-4.84645E+00,
     *  4.24178E-01, 0.00000E+00, 0.00000E+00, 1.86494E-01,-9.56931E-02/
      DATA PBA4/
     *  2.08426E+00, 1.53714E+00,-2.87496E-01, 4.06380E-01,-3.59788E-01,
     * -1.87814E-01, 0.00000E+00, 0.00000E+00, 2.01362E-01,-1.21604E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 7.86304E+00,
     *  2.51878E+00, 2.91455E+00, 4.32308E+00, 6.77054E-02,-2.39125E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.57976E+00,-5.44598E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -5.30593E-01,-5.02237E-01,-2.05258E-01, 2.62263E-01,-2.50195E-01,
     *  4.28151E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         WINF
      DATA PCA1/
     *  0.00000E+00, 1.31026E+01,-4.93171E+01, 2.51045E+01,-1.30531E+01,
     *  6.56421E-01, 2.75633E+01, 4.36433E+00, 1.04638E+00, 5.77365E+01,
     *  0.00000E+00,-6.27766E+00, 2.33010E+00,-1.41351E+01, 2.49653E-01,
     *  0.00000E+00, 0.00000E+00, 8.00000E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 1.03817E-02,-1.70950E+01,-1.92295E+00, 0.00000E+00,
     *  0.00000E+00,-1.17490E+01,-7.14788E-01, 6.72649E+00, 0.00000E+00,
     *  0.00000E+00,-1.57793E+02,-1.70815E+00,-7.92416E+00,-1.67372E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.87973E-01,
     * -1.61602E-01,-1.13832E-01,-7.22447E-01, 2.21119E-02, 0.00000E+00,
     * -3.01967E+00,-1.72798E-01,-5.15055E-03,-1.23477E-02, 3.60805E-03/
      DATA PCA2/
     * -1.36730E+00, 0.00000E+00, 1.24390E-02, 0.00000E+00,-1.36577E+00,
     *  3.18101E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.39334E+01,
     *  1.42088E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,-4.72219E+00,
     * -7.47970E+00,-4.96528E+00, 0.00000E+00, 1.24712E+00,-2.56833E+01,
     * -4.26630E+01, 3.92431E+04,-2.57155E+00,-4.35589E-02, 0.00000E+00,
     *  0.00000E+00, 2.02425E+00,-1.48131E+00,-7.72242E-01, 2.99008E+04,
     *  4.50148E-03, 5.29718E-03,-1.26697E-02, 3.20909E-02, 0.00000E+00,
     *  0.00000E+00, 7.01739E+00, 3.11204E+00, 0.00000E+00, 0.00000E+00,
     * -2.13088E+00, 1.32789E+01, 5.07958E+00, 7.26537E-02, 2.87495E-01,
     *  9.97311E-03,-2.56440E+00, 0.00000E+00, 0.00000E+00, 3.00000E+00/
      DATA PCA3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-9.90073E-03,-3.27333E-02,
     * -4.30379E+01,-2.87643E+01,-5.91793E+00,-1.50460E+02, 0.00000E+00,
     *  0.00000E+00, 6.55038E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 6.18051E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.40484E+00, 5.54554E+00, 0.00000E+00, 0.00000E+00, 7.93810E+00,
     *  1.57192E+00, 1.03971E+00, 9.88279E-01,-4.37662E-02,-2.15763E-02/
      DATA PCA4/
     * -2.31583E+00, 4.32633E+00,-1.12716E+00, 3.38459E-01, 4.66956E-01,
     *  7.18403E-01, 5.80836E-02, 4.12653E-01, 1.04111E-01,-8.30672E-02,
     * -5.55541E+00,-4.97473E+00,-2.03007E+01, 0.00000E+00,-6.06235E-01,
     * -1.73121E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 9.29850E-02,-6.38131E-02,
     *  3.93037E-02, 5.21942E-02, 2.26578E-02, 4.13157E-02, 0.00000E+00,
     *  6.28524E+00, 4.43721E+00,-4.31270E+00, 2.32787E+00, 2.55591E-01,
     *  1.60098E+00,-1.20649E+00, 3.05042E+00,-1.88944E+00, 5.35561E+00,
     *  2.02391E-01, 4.62950E-02, 3.39155E-01, 7.94007E-02, 6.30345E-01,
     *  1.93554E-01, 3.93238E-01, 1.76254E-01,-2.51359E-01,-7.06879E-01/
C       UGN1(1)
      DATA PBB1/
     *  6.22831E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  5.90566E+00, 0.00000E+00, 0.00000E+00,-3.20571E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-8.30368E-01, 1.39396E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.40657E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-4.80790E+00,-1.62744E+00, 2.21119E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBB2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.00000E+00/
      DATA PBB3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.10531E-01,
     * -8.94829E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         UGN1(1)
      DATA PCB1/
     *  5.45009E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -3.60304E+00, 0.00000E+00, 0.00000E+00,-5.04071E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 5.62113E-01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.14657E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 4.65483E-01, 1.73636E+00, 2.21119E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PCB2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.00000E+00/
      DATA PCB3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-8.30769E-01,
     *  7.73649E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         UN1(1)
      DATA PBC1/
     *  6.09940E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.39396E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.21119E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBC2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.00000E+00/
      DATA PBC3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C        UN1(1)
      DATA PCC1/
     *  5.46739E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.21119E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PCC2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.00000E+00/
      DATA PCC3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C        UN1(2)
      DATA PBD1/
     *  4.99007E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  2.59994E+00, 0.00000E+00, 0.00000E+00,-1.78418E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-5.24986E+00, 1.39396E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.77918E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.21119E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBD2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.00000E+00/
      DATA PBD3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.68996E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         UN1(2)
      DATA PCD1/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -7.26156E+00, 0.00000E+00, 0.00000E+00,-4.12416E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-2.88934E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 3.65720E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.21119E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PCD2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.00000E+00/
      DATA PCD3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.01835E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C        UN1(3)
      DATA PBE1/
     *  0.00000E+00,-1.37217E+01, 0.00000E+00, 2.38712E-01,-3.92230E+00,
     *  6.11035E+00,-1.57794E+00,-5.87709E-01, 1.21178E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 5.23202E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-2.22836E+03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-3.94006E+00, 1.39396E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 3.99844E-01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.38936E+00, 2.22534E+00, 2.21119E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBE2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.00000E+00/
      DATA PBE3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 4.35518E-01, 8.40051E-01, 0.00000E+00,-8.88181E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 6.81729E-01, 9.67536E-01,
     *  0.00000E+00,-9.67836E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          UN1(3)
      DATA PCE1/
     *  0.00000E+00,-2.75655E+01,-6.61134E+00, 4.85118E+00, 8.15375E-01,
     * -2.62856E+00, 2.99508E-02,-2.00532E-01,-9.35618E+00, 1.17196E+01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-2.43848E+00, 1.90065E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-3.37525E-01, 1.76471E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PCE2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.00000E+00/
      DATA PCE3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-9.23682E-01,-8.84150E-02, 0.00000E+00,-9.88578E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-1.00747E+00,-1.07468E-02,
     *  0.00000E+00,-3.66376E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C           UN1(4)
      DATA PBF1/
     *  0.00000E+00, 1.02709E+01, 0.00000E+00,-1.42016E+00,-4.90438E+00,
     * -9.11544E+00,-3.80570E+00,-2.09013E+00, 1.32939E+01,-1.28062E+01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.23024E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 3.92126E+02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.39396E-02,
     *  0.00000E+00, 0.00000E+00,-5.56532E+00,-1.27046E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-3.03553E+00,-9.09832E-01, 2.21119E-02, 0.00000E+00,
     *  8.89965E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBF2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 9.19210E-01, 0.00000E+00, 0.00000E+00, 4.00000E+00/
      DATA PBF3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-2.46693E-01, 7.44650E-02, 3.84661E-01, 9.44052E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-2.25083E-01, 1.54206E-01,
     *  4.41303E-01, 8.74742E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C           UN1(4)
      DATA PCF1/
     *  0.00000E+00, 3.61143E+00,-8.24679E+00, 1.70751E+00, 1.16676E+00,
     *  6.24821E+00,-5.68968E-01, 8.53046E-01,-6.94168E+00, 1.04152E+01,
     * -3.70861E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.23336E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 5.33958E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-6.43682E-01,-1.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.00000E+00, 0.00000E+00,-5.47300E-01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-8.58764E-01, 4.72310E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PCF2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.00000E+00/
      DATA PCF3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 3.37325E-01,-3.57698E-02,-6.97393E-01, 1.35387E+01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.78162E-01,-2.33383E-01,
     * -7.12994E-01, 1.29234E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         UN1(5)
      DATA PBG1/
     *  0.00000E+00,-1.71856E+00, 5.32877E+00, 5.33548E-01,-2.66034E+00,
     *  6.76192E-01, 2.25618E+00,-5.78954E-01,-2.69685E+00, 1.21933E+00,
     * -6.13650E+00, 7.79531E-01, 1.63652E+00, 3.63835E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 7.51539E+00,-5.27337E-01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.06625E-01, 1.39396E-02,
     *  0.00000E+00, 0.00000E+00,-1.07240E+00,-8.31257E-01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 7.04016E-01, 0.00000E+00,
     *  7.56158E-01,-4.21268E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 1.02843E+00, 5.21034E-01, 2.21119E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBG2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 4.12504E+00, 1.08459E-01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -3.16261E-01, 0.00000E+00,-1.44288E-01, 0.00000E+00, 4.00000E+00/
      DATA PBG3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.36181E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         UN1(5)
      DATA PCG1/
     *  0.00000E+00, 3.47155E+00, 1.76102E+01, 2.80371E+00,-2.08556E+00,
     *  1.10473E+00, 6.74582E+00,-5.75988E-01, 1.02708E+00,-2.23381E+01,
     *  8.60171E+00, 5.12046E-01,-8.12861E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 9.11036E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 3.89742E+00, 2.01725E-01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 5.06308E-01, 2.04596E-01, 0.00000E+00,
     *  4.40377E+00, 0.00000E+00, 0.00000E+00, 2.20760E+00, 0.00000E+00,
     * -1.36478E+00, 2.38097E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-7.08949E-02,-1.61277E-01, 2.21119E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PCG2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-2.16898E+00,-5.31596E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  2.53060E+00, 0.00000E+00,-7.17287E-01, 0.00000E+00, 4.00000E+00/
      DATA PCG3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.91762E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          UGN1(2 
      DATA PBH1/
     *  0.00000E+00,-7.70936E-01, 1.58158E+00, 3.61790E+00,-1.51748E+00,
     * -5.66098E-01, 1.69393E+00,-4.60489E-01,-8.31527E-01,-4.66437E-01,
     * -1.21750E+00, 0.00000E+00, 0.00000E+00, 1.56505E+02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-5.19321E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.39396E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 3.09223E-01, 1.33715E-01, 2.21119E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBH2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.00000E+00/
      DATA PBH3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          UGN1(2)
      DATA PCH1/
     *  0.00000E+00, 1.72324E-01, 3.08033E-01, 4.55771E-01, 1.46516E-01,
     *  1.97176E-01,-1.53329E-01, 6.91877E-02,-3.07184E-01, 2.65686E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-2.24369E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 4.04079E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  4.99627E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-7.83317E-03,-6.88967E-02, 2.21119E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PCH2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.00000E+00/
      DATA PCH3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          UN2(2)
      DATA PBI1/
     *  0.00000E+00,-7.99767E-01,-3.24774E-01, 7.70975E-01, 6.71796E-01,
     *  5.65483E-01,-2.99727E+00, 3.32448E+00,-9.15018E-01, 5.97656E+00,
     *  0.00000E+00,-1.19515E+00,-8.30457E-01, 3.26074E+00, 0.00000E+00,
     *  0.00000E+00,-1.58365E+00, 7.44825E-02, 5.91372E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-1.41511E-01,-3.01048E+00,
     *  2.35960E+01, 0.00000E+00,-1.70352E+00,-2.39746E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.30488E+00, 0.00000E+00,
     *  5.95132E-01, 5.64301E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 5.30317E-01, 5.66569E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBI2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 5.72367E+00, 1.58411E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.04557E-01, 0.00000E+00,-2.04710E-01, 0.00000E+00, 5.00000E+00/
C          UN2(2)
      DATA PCI1/
     *  0.00000E+00, 6.34487E+00, 9.84162E+00, 3.42136E+00,-5.10607E+00,
     * -8.58745E-02, 3.11501E+00, 5.34570E-01, 1.18027E+00, 4.28027E+00,
     *  4.75123E+00, 6.40947E-01,-4.15165E+00,-1.38154E+01, 0.00000E+00,
     *  0.00000E+00, 1.13145E+01,-5.15954E+00, 0.00000E+00, 0.00000E+00,
     *  1.35576E+01, 0.00000E+00,-5.78982E+00,-2.22043E+00, 3.36776E+00,
     *  3.04791E+01, 0.00000E+00, 2.94709E+00,-4.17536E-01,-1.59855E+00,
     * -2.18320E+00, 1.68269E+01, 0.00000E+00, 1.00829E+00, 0.00000E+00,
     * -6.85096E-01, 2.07822E-01, 3.50168E-01,-3.03662E+01, 0.00000E+00,
     *  0.00000E+00,-1.65726E-01,-8.97831E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-5.24159E+00, 0.00000E+00,-3.52218E+00/
      DATA PCI2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 5.69093E-01,-7.44918E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  2.10865E+00, 0.00000E+00, 1.76776E-01, 1.54755E+00, 5.00000E+00/
C          UN2(3)
      DATA PBJ1/
     *  0.00000E+00, 2.28657E+00, 4.96548E-01, 6.99915E+00,-2.31540E+00,
     * -1.82163E-01,-5.00779E-01, 3.18199E-01,-6.14645E-01, 6.34816E+00,
     *  0.00000E+00, 7.94635E-01,-5.55565E-01, 3.85494E+00, 0.00000E+00,
     *  0.00000E+00,-3.96109E+00, 1.90775E-01, 4.51396E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-5.04618E-01,-4.14385E+00,
     *  2.30244E+01, 0.00000E+00, 1.00689E+00, 5.75680E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 8.56741E-01, 0.00000E+00,
     *  9.54921E-02, 5.56659E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 1.38503E-01, 4.50415E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBJ2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 2.22813E-01,-8.63549E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.37970E-01, 0.00000E+00,-3.25612E-01, 0.00000E+00, 5.00000E+00/
C          UN2(3)
      DATA PCJ1/
     *  0.00000E+00, 5.07608E+00, 3.31479E+00, 3.01548E-01,-1.12100E+00,
     * -7.63711E-02, 2.29748E+00,-1.36699E+00, 7.53433E-01, 3.60702E+01,
     * -1.55266E+00, 1.47382E+00,-2.53895E+00,-1.47720E+01, 0.00000E+00,
     *  0.00000E+00, 1.11787E+01,-1.06256E+01, 0.00000E+00, 0.00000E+00,
     *  7.86391E+00, 0.00000E+00,-8.61020E+00,-1.59313E+00,-5.17013E+00,
     *  1.20468E+00, 0.00000E+00, 5.76509E-01, 9.96195E-01,-1.45539E+00,
     * -1.79950E+01, 8.76957E+00, 0.00000E+00,-1.22863E+00, 0.00000E+00,
     * -6.19019E-01,-1.09571E-01,-4.31325E-02,-4.21981E+01, 0.00000E+00,
     *  0.00000E+00,-1.51519E-01,-1.24067E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-6.39248E+00, 0.00000E+00, 6.64508E-01/
      DATA PCJ2/
     * -7.33184E-01,-9.72031E-03, 1.36789E+00,-8.62311E-01,-3.06395E-03,
     *  2.53354E-01,-2.40918E-01,-4.06932E-02,-5.82223E-01, 0.00000E+00,
     * -8.70285E-01, 7.72318E-01,-6.54213E-01,-2.19231E+01,-1.56509E-01,
     *  2.71745E-01, 5.93538E-01, 2.27757E-01,-5.98215E-01, 3.96457E-01,
     *  2.98705E-01, 1.78618E-01,-5.24538E-01, 1.16439E-01, 7.56829E-02,
     * -4.26809E-01, 5.77187E-01, 8.65450E-01,-7.53614E-01, 1.38381E-01,
     * -1.82265E-01, 2.85263E-01, 4.51322E-01, 1.02775E-01, 3.55731E-01,
     * -4.60896E-01,-3.13037E+01,-2.70818E+00,-7.84847E-01, 0.00000E+00,
     * -1.03473E-01,-3.87649E-01,-1.22613E-01, 0.00000E+00, 0.00000E+00,
     *  8.91325E-01, 0.00000E+00, 1.06189E-01, 9.13030E-02, 5.00000E+00/
C          UN2(4)
      DATA PBK1/
     *  0.00000E+00, 2.94921E+00, 2.79238E+00, 2.58949E+00, 3.56459E-01,
     *  3.12952E-01, 3.34337E+00,-2.83209E+00,-1.05979E+00, 3.92313E+00,
     *  0.00000E+00, 1.73703E-01,-3.23441E-01, 4.15836E+00, 0.00000E+00,
     *  0.00000E+00,-1.77156E+00, 6.44113E-01, 1.88743E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-4.64778E-01,-4.23560E+00,
     *  2.27271E+01, 0.00000E+00,-4.89468E-01, 1.82689E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 4.38217E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 8.62449E-02, 4.46041E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBK2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.40393E-01, 1.01821E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.00000E+00/
C          UN2(4)
      DATA PCK1/
     *  0.00000E+00, 6.04465E+00, 4.50924E+00, 3.84425E-02,-8.70772E-01,
     * -9.55408E-02, 2.28287E+00,-4.37834E-01, 3.57839E-01, 7.20721E+01,
     * -4.41757E+00,-9.13648E-01,-8.71866E-01,-6.26173E+00, 0.00000E+00,
     *  0.00000E+00, 5.92817E+00, 6.15853E+00, 0.00000E+00, 0.00000E+00,
     * -4.89060E+00, 0.00000E+00,-8.30378E+00, 1.07462E-01, 1.08471E+02,
     *  3.39150E+01,-4.57863E+00,-7.18349E-02,-2.71703E-01,-8.96297E+00,
     * -2.37986E+01, 4.11880E+00, 0.00000E+00,-9.95820E-01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-8.91622E+00,-6.85950E+01, 0.00000E+00,
     *  0.00000E+00,-3.62769E-02,-1.65893E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-2.94563E+00, 0.00000E+00, 1.23581E+00/
      DATA PCK2/
     * -6.06026E-01,-6.50229E-01, 1.91330E+00,-1.00314E+00, 1.13346E-01,
     *  4.21885E-01,-3.97688E-01,-2.77437E-01,-6.65893E-01, 0.00000E+00,
     * -1.37646E+00, 1.35171E+00,-9.55595E-01,-1.96450E+01,-2.50039E-01,
     *  5.93389E-01, 9.87131E-01, 5.43559E-01,-1.04322E+00, 6.32546E-01,
     *  3.73259E-01, 5.22657E-01,-5.81400E-01,-1.26425E-01,-1.29843E-01,
     * -5.36598E-01, 8.02402E-01, 9.04347E-01,-1.10799E+00, 1.24800E-01,
     *  1.62487E-02, 2.84237E-01,-1.68866E+00, 5.07723E-01, 5.14161E-01,
     * -4.71292E-01,-3.03487E+01, 4.17455E-01,-1.12591E+00, 0.00000E+00,
     * -3.03544E-01,-6.60313E-01,-1.48606E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.00607E+01, 5.00000E+00/
C          UN2(5)
      DATA PBL1/
     *  0.00000E+00, 2.52207E+00, 3.84550E+00, 1.68023E+00, 7.93489E-01,
     *  3.93725E-02,-2.79707E+00,-4.76621E-01,-1.19972E-01, 3.20454E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 4.17146E+00, 0.00000E+00,
     *  0.00000E+00,-5.30699E-01, 9.14373E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-4.84434E-02, 1.85902E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBL2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.00000E+00/
C          UN2(5)
      DATA PCL1/
     *  0.00000E+00, 1.55386E+01, 4.21418E+00,-9.70151E-01,-8.77326E-01,
     *  2.65813E-02, 1.40164E+00,-9.03874E-01, 3.17281E-03, 9.26891E+01,
     * -4.96004E+00, 0.00000E+00, 0.00000E+00,-4.17851E+00, 0.00000E+00,
     *  0.00000E+00,-1.14760E+01, 2.67744E+00, 0.00000E+00, 0.00000E+00,
     * -1.60056E+01, 0.00000E+00,-7.14647E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-2.89639E+00, 0.00000E+00, 0.00000E+00,-3.88601E+00,
     * -1.65784E+01, 8.44796E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-3.75324E+00,-6.24047E+01, 0.00000E+00,
     *  0.00000E+00,-2.86808E-02,-1.95891E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-3.10534E-01, 0.00000E+00,-3.37448E+00/
      DATA PCL2/
     *  1.63964E-02,-1.45191E+00, 1.85618E+00,-9.61979E-01, 3.93783E-01,
     *  4.21681E-01,-5.30254E-01,-2.96232E-01,-7.55211E-01, 0.00000E+00,
     * -1.85443E+00, 1.88047E+00,-1.07818E+00,-1.35373E+01,-3.05785E-01,
     *  7.82159E-01, 1.32586E+00, 2.34413E-01,-7.47152E-01, 9.92893E-01,
     * -2.80110E-02, 3.61747E-01,-4.16280E-01,-3.46427E-01,-5.76431E-01,
     * -2.13906E-01, 9.51184E-01, 3.69403E-01,-1.35563E+00, 6.59534E-02,
     *  1.39764E-01, 4.50687E-01,-1.22025E+00, 5.73280E-02, 7.49303E-01,
     * -8.37947E-01,-3.01332E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -4.36697E-01,-7.76068E-01,-1.41680E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.21958E+01, 5.00000E+00/
C          UN2(6)
      DATA PBM1/
     *  0.00000E+00, 3.13842E+00,-8.20417E-01, 3.72282E+00,-5.20477E-01,
     * -3.61867E-01,-2.92604E+00, 3.13013E-01,-1.38865E-01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.30060E+01, 0.00000E+00,
     *  0.00000E+00, 1.67696E+00, 9.85990E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-8.46922E-02, 5.59429E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBM2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.00000E+00/
C          UN2(6)
      DATA PCM1/
     *  0.00000E+00, 1.78539E+01, 1.07314E+01,-1.13212E+00, 1.59867E-02,
     *  1.53736E-01, 2.25710E+00,-9.39080E-01,-9.72620E-02, 9.89789E+01,
     * -5.17469E+00, 0.00000E+00, 0.00000E+00,-2.98597E+00, 0.00000E+00,
     *  0.00000E+00,-2.04707E+01, 4.92899E+00, 0.00000E+00, 0.00000E+00,
     * -1.44316E+01, 0.00000E+00,-3.31557E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-6.22743E+00, 0.00000E+00, 0.00000E+00,-4.34344E+00,
     * -8.29640E+00,-3.03800E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 2.79387E+00,-5.23752E+01, 0.00000E+00,
     *  0.00000E+00,-2.59963E-02,-1.73426E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-5.37220E+00, 0.00000E+00,-6.53478E-01/
      DATA PCM2/
     *  3.48181E-01,-1.88980E+00, 1.47787E+00,-7.92670E-01, 6.49224E-01,
     *  5.96079E-01,-1.04901E+00,-5.24003E-01,-6.77311E-01, 0.00000E+00,
     * -2.26873E+00, 2.80910E+00,-9.84994E-01,-6.79661E+00,-3.71975E-01,
     *  1.13310E+00, 1.57164E+00, 2.15176E-01,-5.58583E-01, 1.16045E+00,
     *  2.05395E-02, 2.27714E-01, 1.41203E-01,-3.92231E-01,-8.82859E-01,
     *  4.90400E-01, 1.14013E+00,-2.25250E-01,-1.64930E+00, 5.73434E-02,
     *  1.89857E-01, 4.31221E-01,-1.35345E+00,-2.94189E-01, 6.87530E-01,
     * -7.78284E-01,-2.88975E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -3.98115E-01,-7.40699E-01,-8.28264E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.02069E+00, 5.00000E+00/
C          UN2(7)
      DATA PBN1/
     *  0.00000E+00, 2.08818E+00,-1.96235E+00, 4.55317E+00,-1.76012E+00,
     * -4.75258E-01,-1.44220E+00,-3.28566E-01,-1.41177E-01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.49146E+01, 0.00000E+00,
     *  0.00000E+00, 1.73222E+00, 9.91286E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.35468E-01, 1.91833E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBN2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.00000E+00/
C          UN2(7)
      DATA PCN1/
     *  0.00000E+00, 1.25645E+01, 2.43937E+01,-4.89691E-01,-5.46437E-01,
     *  1.22200E-01, 2.89309E+00,-2.85509E-01,-2.27122E-01, 9.54192E+01,
     * -4.07394E+00, 0.00000E+00, 0.00000E+00,-3.04354E+00, 0.00000E+00,
     *  0.00000E+00,-2.36547E+01, 1.04903E+01, 0.00000E+00, 0.00000E+00,
     * -8.32274E+00, 0.00000E+00,-3.34712E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-7.95953E+00, 0.00000E+00, 0.00000E+00,-5.83474E+00,
     * -1.48074E+00, 1.02268E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 6.19470E+00,-3.90767E+01, 0.00000E+00,
     *  0.00000E+00,-3.58136E-03, 1.22289E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-8.49787E+00, 0.00000E+00,-3.97498E+00/
      DATA PCN2/
     *  3.79580E-01,-1.93595E+00, 2.89114E+00,-4.73457E-01, 7.67548E-01,
     *  5.66859E-01,-1.28683E+00,-8.37174E-01,-3.48022E-01, 0.00000E+00,
     * -2.62865E+00, 3.50575E+00,-7.93257E-01,-8.10692E-01,-4.99450E-01,
     *  1.56654E+00, 1.63039E+00, 7.58900E-02,-4.30952E-01, 1.23068E+00,
     *  1.06404E-01, 4.73870E-02, 5.50559E-01,-4.11375E-01,-9.94162E-01,
     *  1.35025E+00, 1.26053E+00,-7.34502E-01,-2.01952E+00, 2.05398E-01,
     * -4.77248E-02, 2.41549E-01,-9.32522E-01,-5.63663E-01, 5.34833E-01,
     * -5.77563E-01,-2.65033E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -2.42317E-01,-7.33679E-01,-7.85537E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.56842E-01, 5.00000E+00/
C          UN2(8)
      DATA PBO1/
     *  0.00000E+00, 7.00409E-01,-4.17017E-01, 3.24757E+00,-1.28352E+00,
     * -4.23875E-01, 1.64346E+00,-1.20855E+00,-7.65316E-01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-3.39417E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 2.68534E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.56444E-01,-4.60043E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBO2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.00000E+00/
C          UN2(8)
      DATA PCO1/
     *  0.00000E+00, 7.30129E+00, 3.14811E+01,-7.06834E-02,-2.96193E-01,
     *  1.73817E-01, 1.62127E+00,-2.71556E-01,-2.05844E-01, 8.02088E+01,
     * -1.86956E-01, 0.00000E+00, 0.00000E+00,-9.43641E-01,-3.24716E+00,
     *  0.00000E+00,-2.32748E+01, 1.96724E+01, 0.00000E+00, 0.00000E+00,
     * -3.95949E+00, 0.00000E+00, 5.44787E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.00161E+01, 0.00000E+00, 0.00000E+00,-4.57422E+00,
     *  4.31304E+00, 1.49868E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 5.99489E+00,-2.82120E+01, 0.00000E+00,
     *  0.00000E+00, 4.03624E-02, 1.19463E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.39050E+01, 0.00000E+00,-2.65634E+00/
      DATA PCO2/
     *  6.37036E-01,-1.77461E+00, 3.03103E+00,-1.49839E-01, 7.02027E-01,
     *  6.08841E-01,-9.27289E-01,-8.52362E-01, 5.61723E-01, 0.00000E+00,
     * -2.72061E+00, 3.66183E+00,-2.54943E-01, 2.94668E+00,-3.57898E-01,
     *  1.71858E+00, 1.58782E+00,-2.42995E-01,-3.57783E-01, 1.20157E+00,
     *  2.58895E-01,-1.05773E-01, 5.79397E-01,-3.30395E-01,-4.03569E-01,
     *  1.99175E+00, 1.21688E+00,-8.64350E-01,-1.95569E+00, 4.61136E-01,
     * -8.61382E-02, 3.38859E-01, 0.00000E+00,-5.78864E-01, 4.46659E-01,
     * -4.57428E-01,-1.99920E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -1.19841E-01,-4.56968E-01, 2.00180E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-1.07368E+00, 5.00000E+00/
C          UN2(9)
      DATA PBP1/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.75863E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 3.18522E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBP2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.00000E+00/
C          UN2(9)
      DATA PCP1/
     *  0.00000E+00, 4.61019E-02, 3.50615E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 6.15349E+01,
     *  4.28634E+00, 0.00000E+00, 0.00000E+00, 6.03982E+00,-4.72305E+00,
     *  0.00000E+00,-1.43678E+01, 3.62580E+01, 0.00000E+00, 0.00000E+00,
     *  1.26574E+00, 0.00000E+00,-2.77285E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.14802E+01, 0.00000E+00, 0.00000E+00,-1.11940E+01,
     * -1.39535E+00, 2.63070E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.53024E+00,-2.14609E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.26956E+01, 0.00000E+00, 5.49926E+00/
      DATA PCP2/
     *  9.80142E-01,-1.19016E+00, 2.75110E+00, 4.23423E-01, 5.89893E-01,
     *  4.94288E-01,-5.25954E-01,-8.51760E-01, 1.62676E+00, 0.00000E+00,
     * -1.90027E+00, 3.19950E+00, 4.72739E-01, 7.04179E+00,-1.43685E-03,
     *  1.43219E+00, 1.32136E+00,-2.92744E-03,-3.43680E-01, 7.75735E-01,
     *  6.92202E-01,-1.45519E-01, 6.97813E-02,-3.11588E-01, 6.65750E-01,
     *  2.33809E+00, 1.06694E+00,-5.77590E-01,-1.33717E+00, 8.13367E-01,
     * -5.05737E-01, 5.99169E-01,-8.83386E-01,-4.38123E-01, 2.63649E-01,
     * -3.03448E-01,-1.28190E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.45478E-02, 1.45491E-01, 2.40080E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-3.86910E+00, 5.00000E+00/
C          UN2(10)
      DATA PBQ1/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.10647E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 3.13252E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBQ2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.00000E+00/
C          UN2(10)
      DATA PCQ1/
     *  0.00000E+00,-3.03260E+00, 3.15488E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.42798E+01,
     *  7.08849E+00, 0.00000E+00, 0.00000E+00, 1.64773E+01,-6.86505E+00,
     *  0.00000E+00,-6.27112E+00, 3.78373E+01, 0.00000E+00, 0.00000E+00,
     *  2.97763E+00, 0.00000E+00,-3.44134E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.19424E+01, 0.00000E+00, 0.00000E+00,-1.64645E+01,
     * -2.27053E+00, 3.82330E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 1.33140E-01,-2.08131E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-7.04687E+00, 0.00000E+00, 6.52184E+00/
      DATA PCQ2/
     *  7.31799E-01,-2.75395E-01, 1.92467E+00, 8.71269E-01, 3.72836E-01,
     *  3.04967E-01, 7.72480E-02,-5.08596E-01, 1.99828E+00, 0.00000E+00,
     * -5.51169E-01, 2.12420E+00, 8.96069E-01, 1.12092E+01,-4.30438E-02,
     *  7.38391E-01, 6.12050E-01, 3.62981E-02,-1.02054E-01, 1.82404E-01,
     *  3.70643E-01,-1.68899E-01,-1.79628E-01,-1.21117E-01, 1.45823E+00,
     *  2.04352E+00, 7.83711E-01,-3.42979E-02,-2.31363E-01, 7.11253E-01,
     * -3.16353E-01, 6.21069E-01,-1.05676E+00,-4.03488E-01, 4.11595E-01,
     * -2.12535E-01,-6.51453E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.48238E-01, 6.38716E-01, 2.99311E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-1.01846E+00, 5.00000E+00/
C          UN2(11)
      DATA PBR1/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.21764E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 6.77475E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBR2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.00000E+00/
C          UN2(11)
      DATA PCR1/
     *  0.00000E+00,-1.74115E+00, 2.66621E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.13017E+01,
     *  6.86985E+00, 0.00000E+00, 0.00000E+00, 2.08835E+01,-7.86030E+00,
     *  0.00000E+00,-3.77141E+00, 3.87788E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 1.31580E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-9.98927E+00, 0.00000E+00, 0.00000E+00,-1.71002E+01,
     * -9.88358E-01, 4.47756E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 5.95029E-01,-2.11313E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-3.84164E+00, 0.00000E+00, 0.00000E+00/
      DATA PCR2/
     *  3.07191E-01, 4.79094E-02, 6.72159E-01, 5.54185E-01, 1.82847E-01,
     * -1.23768E-02, 1.91637E-01,-2.89429E-02, 1.18297E+00, 0.00000E+00,
     *  2.37450E-01, 9.23551E-01, 6.05670E-01, 1.35990E+01,-1.64210E-01,
     *  5.38355E-03,-4.91246E-02,-1.06966E-01,-2.09635E-01,-3.23023E-02,
     * -3.41663E-02,-3.48871E-02,-2.62450E-01, 2.21492E-01, 1.43749E+00,
     *  1.08677E+00, 3.97778E-01, 3.61526E-01, 5.55950E-01, 3.53058E-01,
     * -5.93339E-02, 4.14203E-01,-6.05024E-01,-1.38714E-01, 2.78897E-01,
     * -8.92889E-02,-3.59033E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  9.90623E-02, 4.36170E-01, 7.95418E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-1.11426E+00, 5.00000E+00/
C          UN2(12)
      DATA PBS1/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 3.07320E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 1.60738E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBS2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.00000E+00/
C          UN2(12)
      DATA PCS1/
     *  0.00000E+00, 1.26217E+01, 2.30787E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.00029E+01,
     * -2.88682E+00, 0.00000E+00, 0.00000E+00, 2.09439E+01,-4.56923E+00,
     *  0.00000E+00,-2.15929E+00, 3.87149E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-7.98039E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-6.63423E+00, 0.00000E+00, 0.00000E+00,-5.84850E+00,
     *  3.72111E+00, 4.52300E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 3.21872E-01, 0.00000E+00, 0.00000E+00/
      DATA PCS2/
     *  1.09405E-02,-4.35341E-02, 8.00586E-02, 1.48577E-01, 1.01602E-01,
     * -1.01104E-01,-1.98993E-02, 3.51174E-02, 2.41112E-01, 0.00000E+00,
     *  2.76479E-01, 1.97043E-01, 2.68708E-01, 1.39832E+01,-1.56638E-01,
     * -2.39101E-01,-1.50605E-01,-2.17139E-01,-2.59057E-01,-4.36362E-01,
     * -1.43496E-01, 7.51305E-02,-2.40850E-01, 1.34858E-01, 7.59193E-01,
     *  3.52708E-01, 1.29922E-01, 3.27957E-01, 5.35491E-01, 1.19120E-01,
     * -2.94029E-02, 1.76113E-01,-6.51597E-01, 3.61575E-02, 4.26836E-02,
     * -2.29297E-02,-4.27373E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -2.78548E-02, 5.77322E-02,-1.02411E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.00000E+00/
C          UN2(13)
      DATA PBT1/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 3.69447E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 2.34073E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBT2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.00000E+00/
C          UN2(13)
      DATA PCT1/
     *  0.00000E+00, 1.22096E+01, 1.92342E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 8.13667E+00,
     * -6.19078E+00, 0.00000E+00, 0.00000E+00, 2.37009E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-7.87365E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.12371E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-2.76047E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 1.85864E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PCT2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.00000E+00/
C          UN2(14)
      DATA PBU1/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.01008E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 2.21469E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PBU2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.00000E+00/
C          UN2(14)
      DATA PCU1/
     *  0.00000E+00,-1.40697E+00, 6.88709E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 3.67624E+02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 1.58312E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-2.46486E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.90327E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 1.13248E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PCU2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.00000E+00/
      END
      
C ****************************************************************
C *****************           MSISE00            *****************
C ****************************************************************
      SUBROUTINE GTD7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
C
C     NRLMSISE-00
C     -----------
C        Neutral Atmosphere Empirical Model from the surface to lower
C        exosphere
C
C        NEW FEATURES:
C          *Extensive satellite drag database used in model generation
C          *Revised O2 (and O) in lower thermosphere
C          *Additional nonlinear solar activity term
C          *"ANOMALOUS OXYGEN" NUMBER DENSITY, OUTPUT D(9)
C           At high altitudes (> 500 km), hot atomic oxygen or ionized
C           oxygen can become appreciable for some ranges of subroutine
C           inputs, thereby affecting drag on satellites and debris. We
C           group these species under the term "anomalous oxygen," since
C           their individual variations are not presently separable with
C           the drag data used to define this model component.
C
C        SUBROUTINES FOR SPECIAL OUTPUTS:
C        
C        HIGH ALTITUDE DRAG: EFFECTIVE TOTAL MASS DENSITY 
C        (SUBROUTINE GTD7D, OUTPUT D(6))
C           For atmospheric drag calculations at altitudes above 500 km,
C           call SUBROUTINE GTD7D to compute the "effective total mass
C           density" by including contributions from "anomalous oxygen."
C           See "NOTES ON OUTPUT VARIABLES" below on D(6).
C
C        PRESSURE GRID (SUBROUTINE GHP7)
C          See subroutine GHP7 to specify outputs at a pressure level
C          rather than at an altitude.
C
C        OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)
C 
C     INPUT VARIABLES:
C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
C              (Year ignored in current model)
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C                 MASS 17 IS Anomalous O ONLY.)
C
C     NOTES ON INPUT VARIABLES: 
C        UT, Local Time, and Longitude are used independently in the
C        model and are not of equal importance for every situation.  
C        For the most physically realistic calculation these three
C        variables should be consistent (STL=SEC/3600+GLONG/15).
C        The Equation of Time departures from the above formula
C        for apparent local time can be included if available but
C        are of minor importance.
c
C        F107 and F107A values used to generate the model correspond
C        to the 10.7 cm radio flux at the actual distance of the Earth
C        from the Sun rather than the radio flux at 1 AU. The following
C        site provides both classes of values:
C        ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
C
C        F107, F107A, and AP effects are neither large nor well
C        established below 80 km and these parameters should be set to
C        150., 150., and 4. respectively.
C
C     OUTPUT VARIABLES:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)                       
C        D(6) - TOTAL MASS DENSITY(GM/CM3)
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
C     NOTES ON OUTPUT VARIABLES:
C        TO GET OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.) 
C
C        O, H, and N are set to zero below 72.5 km
C
C        T(1), Exospheric temperature, is set to global average for
C        altitudes below 120 km. The 120 km gradient is left at global
C        average value for altitudes below 72 km.
C
C        D(6), TOTAL MASS DENSITY, is NOT the same for subroutines GTD7 
C        and GTD7D
C
C          SUBROUTINE GTD7 -- D(6) is the sum of the mass densities of the
C          species labeled by indices 1-5 and 7-8 in output variable D.
C          This includes He, O, N2, O2, Ar, H, and N but does NOT include
C          anomalous oxygen (species index 9).
C
C          SUBROUTINE GTD7D -- D(6) is the "effective total mass density
C          for drag" and is the sum of the mass densities of all species
C          in this model, INCLUDING anomalous oxygen.
C        
C     SWITCHES: The following is for test and special purposes:
C          
C        TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SW),
C        WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1. 
C        FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
C        FOR THE FOLLOWING VARIATIONS
C               1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
C               3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
C               5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
C               7 - DIURNAL               8 - SEMIDIURNAL
C               9 - DAILY AP             10 - ALL UT/LONG EFFECTS
C              11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
C              13 - MIXED AP/UT/LONG     14 - TERDIURNAL
C              15 - DEPARTURES FROM DIFFUSIVE EQUILIBRIUM
C              16 - ALL TINF VAR         17 - ALL TLB VAR
C              18 - ALL TN1 VAR           19 - ALL S VAR
C              20 - ALL TN2 VAR           21 - ALL NLB VAR
C              22 - ALL TN3 VAR           23 - TURBO SCALE HEIGHT VAR
C
C        To get current values of SW: CALL TRETRV(SW)
C
      DIMENSION D(9),T(2),AP(7),DS(9),TS(2)
      DIMENSION ZN3(5),ZN2(4),SV(25)
      COMMON/GTS3C/TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
     & ,G0,RL,DD,DB14,TR12
      COMMON/MESO7/TN1(5),TN2(4),TN3(5),TGN1(2),TGN2(2),TGN3(2)
      COMMON/LOWER7/PTM(10),PDM(10,8)
      COMMON/PARM7/PT(150),PD(150,9),PS(150),PDL(25,2),PTL(100,4),
     $ PMA(100,10),SAM(100)
      COMMON/DATIM7/ISD(3),IST(2),NAM(2)
      COMMON/DATIME/ISDATE(3),ISTIME(2),NAME(2)
      COMMON/CSW/SW(25),ISW,SWC(25)
      COMMON/MAVG7/PAVGM(10)
      COMMON/DMIX/DM04,DM16,DM28,DM32,DM40,DM01,DM14
      COMMON/PARMB/GSURF,RE
      COMMON/METSEL/IMR
      SAVE
      EXTERNAL GTD7BK
      DATA MN3/5/,ZN3/32.5,20.,15.,10.,0./
      DATA MN2/4/,ZN2/72.5,55.,45.,32.5/
      DATA ZMIX/62.5/,ALAST/99999./,MSSL/-999/
      DATA SV/25*1./
      IF(ISW.NE.64999) CALL TSELEC(SV)
C      Put identification data into common/datime/
      DO 1 I=1,3
        ISDATE(I)=ISD(I)
    1 CONTINUE
      DO 2 I=1,2
        ISTIME(I)=IST(I)
        NAME(I)=NAM(I)
    2 CONTINUE
C
C        Test for changed input
      V1=VTST7(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,1)
C       Latitude variation of gravity (none for SW(2)=0)
      XLAT=GLAT
      IF(SW(2).EQ.0) XLAT=45.
      CALL GLATF(XLAT,GSURF,RE)
C
      XMM=PDM(5,3)
C
C       THERMOSPHERE/MESOSPHERE (above ZN2(1))
      ALTT=AMAX1(ALT,ZN2(1))
      MSS=MASS
C       Only calculate N2 in thermosphere if alt in mixed region
      IF(ALT.LT.ZMIX.AND.MASS.GT.0) MSS=28
C       Only calculate thermosphere if input parameters changed
C         or altitude above ZN2(1) in mesosphere
      IF(V1.EQ.1..OR.ALT.GT.ZN2(1).OR.ALAST.GT.ZN2(1).OR.MSS.NE.MSSL)
     $ THEN
        CALL GTS7(IYD,SEC,ALTT,GLAT,GLONG,STL,F107A,F107,AP,MSS,DS,TS)
        DM28M=DM28
C         metric adjustment
        IF(IMR.EQ.1) DM28M=DM28*1.E6
        MSSL=MSS
      ENDIF
      T(1)=TS(1)
      T(2)=TS(2)
      IF(ALT.GE.ZN2(1)) THEN
        DO 5 J=1,9
          D(J)=DS(J)
    5   CONTINUE
        GOTO 10
      ENDIF
C
C       LOWER MESOSPHERE/UPPER STRATOSPHERE [between ZN3(1) and ZN2(1)]
C         Temperature at nodes and gradients at end nodes
C         Inverse temperature a linear function of spherical harmonics
C         Only calculate nodes if input changed
       IF(V1.EQ.1..OR.ALAST.GE.ZN2(1)) THEN
        TGN2(1)=TGN1(2)
        TN2(1)=TN1(5)
        TN2(2)=PMA(1,1)*PAVGM(1)/(1.-SW(20)*GLOB7S(PMA(1,1)))
        TN2(3)=PMA(1,2)*PAVGM(2)/(1.-SW(20)*GLOB7S(PMA(1,2)))
        TN2(4)=PMA(1,3)*PAVGM(3)/(1.-SW(20)*SW(22)*GLOB7S(PMA(1,3)))
        TGN2(2)=PAVGM(9)*PMA(1,10)*(1.+SW(20)*SW(22)*GLOB7S(PMA(1,10)))
     $  *TN2(4)*TN2(4)/(PMA(1,3)*PAVGM(3))**2
        TN3(1)=TN2(4)
       ENDIF
       IF(ALT.GE.ZN3(1)) GOTO 6
C
C       LOWER STRATOSPHERE AND TROPOSPHERE [below ZN3(1)]
C         Temperature at nodes and gradients at end nodes
C         Inverse temperature a linear function of spherical harmonics
C         Only calculate nodes if input changed
        IF(V1.EQ.1..OR.ALAST.GE.ZN3(1)) THEN
         TGN3(1)=TGN2(2)
         TN3(2)=PMA(1,4)*PAVGM(4)/(1.-SW(22)*GLOB7S(PMA(1,4)))
         TN3(3)=PMA(1,5)*PAVGM(5)/(1.-SW(22)*GLOB7S(PMA(1,5)))
         TN3(4)=PMA(1,6)*PAVGM(6)/(1.-SW(22)*GLOB7S(PMA(1,6)))
         TN3(5)=PMA(1,7)*PAVGM(7)/(1.-SW(22)*GLOB7S(PMA(1,7)))
         TGN3(2)=PMA(1,8)*PAVGM(8)*(1.+SW(22)*GLOB7S(PMA(1,8)))
     $   *TN3(5)*TN3(5)/(PMA(1,7)*PAVGM(7))**2
        ENDIF
    6   CONTINUE
        IF(MASS.EQ.0) GOTO 50
C          LINEAR TRANSITION TO FULL MIXING BELOW ZN2(1)
        DMC=0
        IF(ALT.GT.ZMIX) DMC=1.-(ZN2(1)-ALT)/(ZN2(1)-ZMIX)
        DZ28=DS(3)
C      ***** N2 DENSITY ****
        DMR=DS(3)/DM28M-1.
        D(3)=DENSM(ALT,DM28M,XMM,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)
        D(3)=D(3)*(1.+DMR*DMC)
C      ***** HE DENSITY ****
        D(1)=0
        IF(MASS.NE.4.AND.MASS.NE.48) GOTO 204
          DMR=DS(1)/(DZ28*PDM(2,1))-1.
          D(1)=D(3)*PDM(2,1)*(1.+DMR*DMC)
  204   CONTINUE
C      **** O DENSITY ****
        D(2)=0
        D(9)=0
  216   CONTINUE
C      ***** O2 DENSITY ****
        D(4)=0
        IF(MASS.NE.32.AND.MASS.NE.48) GOTO 232
          DMR=DS(4)/(DZ28*PDM(2,4))-1.
          D(4)=D(3)*PDM(2,4)*(1.+DMR*DMC)
  232   CONTINUE
C      ***** AR DENSITY ****
        D(5)=0
        IF(MASS.NE.40.AND.MASS.NE.48) GOTO 240
          DMR=DS(5)/(DZ28*PDM(2,5))-1.
          D(5)=D(3)*PDM(2,5)*(1.+DMR*DMC)
  240   CONTINUE
C      ***** HYDROGEN DENSITY ****
        D(7)=0
C      ***** ATOMIC NITROGEN DENSITY ****
        D(8)=0
C
C       TOTAL MASS DENSITY
C
        IF(MASS.EQ.48) THEN
         D(6) = 1.66E-24*(4.*D(1)+16.*D(2)+28.*D(3)+32.*D(4)+40.*D(5)+
     &       D(7)+14.*D(8))  
         IF(IMR.EQ.1) D(6)=D(6)/1000.
         ENDIF
         T(2)=TZ
   10 CONTINUE
      GOTO 90
   50 CONTINUE
      DD=DENSM(ALT,1.,0,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)                
      T(2)=TZ
   90 CONTINUE
      ALAST=ALT
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GTD7D(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,
     $ D,T)
C
C     NRLMSISE-00
C     -----------
C        This subroutine provides Effective Total Mass Density for
C        output D(6) which includes contributions from "anomalous
C        oxygen" which can affect satellite drag above 500 km.  This
C        subroutine is part of the distribution package for the 
C        Neutral Atmosphere Empirical Model from the surface to lower
C        exosphere.  See subroutine GTD7 for more extensive comments.
C
C     INPUT VARIABLES:
C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
C              (Year ignored in current model)
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C                 MASS 17 IS Anomalous O ONLY.)
C
C     NOTES ON INPUT VARIABLES: 
C        UT, Local Time, and Longitude are used independently in the
C        model and are not of equal importance for every situation.  
C        For the most physically realistic calculation these three
C        variables should be consistent (STL=SEC/3600+GLONG/15).
C        The Equation of Time departures from the above formula
C        for apparent local time can be included if available but
C        are of minor importance.
c
C        F107 and F107A values used to generate the model correspond
C        to the 10.7 cm radio flux at the actual distance of the Earth
C        from the Sun rather than the radio flux at 1 AU.
C
C     OUTPUT VARIABLES:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)                       
C        D(6) - TOTAL MASS DENSITY(GM/CM3) [includes anomalous oxygen]
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
      DIMENSION D(9),T(2),AP(7),DS(9),TS(2)
      COMMON/METSEL/IMR
      CALL GTD7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
C       TOTAL MASS DENSITY
C
        IF(MASS.EQ.48) THEN
         D(6) = 1.66E-24*(4.*D(1)+16.*D(2)+28.*D(3)+32.*D(4)+40.*D(5)+
     &       D(7)+14.*D(8)+16.*D(9))  
         IF(IMR.EQ.1) D(6)=D(6)/1000.
         ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GHP7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,
     $  D,T,PRESS)
C       FIND ALTITUDE OF PRESSURE SURFACE (PRESS) FROM GTD7
C     INPUT:
C        IYD - YEAR AND DAY AS YYDDD
C        SEC - UT(SEC)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS)
C        F107A - 3 MONTH AVERAGE OF F10.7 FLUX
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 59 HRS PRIOR
C                    TO CURRENT TIME
C        PRESS - PRESSURE LEVEL(MB)
C     OUTPUT:
C        ALT - ALTITUDE(KM) 
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)
C        D(6) - TOTAL MASS DENSITY(GM/CM3)
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        D(9) - HOT O NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
      COMMON/PARMB/GSURF,RE
      COMMON/METSEL/IMR
      DIMENSION D(9),T(2),AP(7)
      SAVE
      DATA BM/1.3806E-19/,RGAS/831.4/
      DATA TEST/.00043/,LTEST/12/
      PL=ALOG10(PRESS)
C      Initial altitude estimate
      IF(PL.GE.-5.) THEN
         IF(PL.GT.2.5) ZI=18.06*(3.00-PL)
         IF(PL.GT..75.AND.PL.LE.2.5) ZI=14.98*(3.08-PL)
         IF(PL.GT.-1..AND.PL.LE..75) ZI=17.8*(2.72-PL)
         IF(PL.GT.-2..AND.PL.LE.-1.) ZI=14.28*(3.64-PL)
         IF(PL.GT.-4..AND.PL.LE.-2.) ZI=12.72*(4.32-PL)
         IF(PL.LE.-4.) ZI=25.3*(.11-PL)
         IDAY=MOD(IYD,1000)
         CL=GLAT/90.
         CL2=CL*CL
         IF(IDAY.LT.182) CD=1.-IDAY/91.25
         IF(IDAY.GE.182) CD=IDAY/91.25-3.
         CA=0
         IF(PL.GT.-1.11.AND.PL.LE.-.23) CA=1.0
         IF(PL.GT.-.23) CA=(2.79-PL)/(2.79+.23)
         IF(PL.LE.-1.11.AND.PL.GT.-3.) CA=(-2.93-PL)/(-2.93+1.11)
         Z=ZI-4.87*CL*CD*CA-1.64*CL2*CA+.31*CA*CL
      ENDIF
      IF(PL.LT.-5.) Z=22.*(PL+4.)**2+110
C      ITERATION LOOP
      L=0
   10 CONTINUE
        L=L+1
        CALL GTD7(IYD,SEC,Z,GLAT,GLONG,STL,F107A,F107,AP,48,D,T)
        XN=D(1)+D(2)+D(3)+D(4)+D(5)+D(7)+D(8)
        P=BM*XN*T(2)
        IF(IMR.EQ.1) P=P*1.E-6
        DIFF=PL-ALOG10(P)
        IF(ABS(DIFF).LT.TEST .OR. L.EQ.LTEST) GOTO 20
        XM=D(6)/XN/1.66E-24
        IF(IMR.EQ.1) XM = XM*1.E3
        G=GSURF/(1.+Z/RE)**2
        SH=RGAS*T(2)/(XM*G)
C         New altitude estimate using scale height
        IF(L.LT.6) THEN
          Z=Z-SH*DIFF*2.302
        ELSE
          Z=Z-SH*DIFF
        ENDIF
        GOTO 10
   20 CONTINUE
      IF(L.EQ.LTEST) WRITE(6,100) PRESS,DIFF
  100 FORMAT(1X,29HGHP7 NOT CONVERGING FOR PRESS, 1PE12.2,E12.2)
      ALT=Z
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GLATF(LAT,GV,REFF)
C      CALCULATE LATITUDE VARIABLE GRAVITY (GV) AND EFFECTIVE
C      RADIUS (REFF)
      REAL LAT
      SAVE
      DATA DGTR/1.74533E-2/
      C2 = COS(2.*DGTR*LAT)
      GV = 980.616*(1.-.0026373*C2)
      REFF = 2.*GV/(3.085462E-6 + 2.27E-9*C2)*1.E-5
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION VTST7(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,IC)
C       Test if geophysical variables or switches changed and save
C       Return 0 if unchanged and 1 if changed
      DIMENSION AP(7),IYDL(2),SECL(2),GLATL(2),GLL(2),STLL(2)
      DIMENSION FAL(2),FL(2),APL(7,2),SWL(25,2),SWCL(25,2)
      COMMON/CSW/SW(25),ISW,SWC(25)
      SAVE
      DATA IYDL/2*-999/,SECL/2*-999./,GLATL/2*-999./,GLL/2*-999./
      DATA STLL/2*-999./,FAL/2*-999./,FL/2*-999./,APL/14*-999./
      DATA SWL/50*-999./,SWCL/50*-999./
      VTST7=0
      IF(IYD.NE.IYDL(IC)) GOTO 10
      IF(SEC.NE.SECL(IC)) GOTO 10
      IF(GLAT.NE.GLATL(IC)) GOTO 10
      IF(GLONG.NE.GLL(IC)) GOTO 10
      IF(STL.NE.STLL(IC)) GOTO 10
      IF(F107A.NE.FAL(IC)) GOTO 10
      IF(F107.NE.FL(IC)) GOTO 10
      DO 5 I=1,7
        IF(AP(I).NE.APL(I,IC)) GOTO 10
    5 CONTINUE
      DO 7 I=1,25
        IF(SW(I).NE.SWL(I,IC)) GOTO 10
        IF(SWC(I).NE.SWCL(I,IC)) GOTO 10
    7 CONTINUE
      GOTO 20
   10 CONTINUE
      VTST7=1
      IYDL(IC)=IYD
      SECL(IC)=SEC
      GLATL(IC)=GLAT
      GLL(IC)=GLONG
      STLL(IC)=STL
      FAL(IC)=F107A
      FL(IC)=F107
      DO 15 I=1,7
        APL(I,IC)=AP(I)
   15 CONTINUE
      DO 16 I=1,25
        SWL(I,IC)=SW(I)
        SWCL(I,IC)=SWC(I)
   16 CONTINUE
   20 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GTS7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
C
C     Thermospheric portion of NRLMSISE-00
C     See GTD7 for more extensive comments
C
C        OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)
C 
C     INPUT VARIABLES:
C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
C              (Year ignored in current model)
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM) (>72.5 km)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C                 MASS 17 IS Anomalous O ONLY.)
C
C     NOTES ON INPUT VARIABLES: 
C        UT, Local Time, and Longitude are used independently in the
C        model and are not of equal importance for every situation.  
C        For the most physically realistic calculation these three
C        variables should be consistent (STL=SEC/3600+GLONG/15).
C        The Equation of Time departures from the above formula
C        for apparent local time can be included if available but
C        are of minor importance.
c
C        F107 and F107A values used to generate the model correspond
C        to the 10.7 cm radio flux at the actual distance of the Earth
C        from the Sun rather than the radio flux at 1 AU. The following
C        site provides both classes of values:
C        ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
C
C        F107, F107A, and AP effects are neither large nor well
C        established below 80 km and these parameters should be set to
C        150., 150., and 4. respectively.
C
C     OUTPUT VARIABLES:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)                       
C        D(6) - TOTAL MASS DENSITY(GM/CM3) [Anomalous O NOT included]
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
      DIMENSION ZN1(5),ALPHA(9)
      COMMON/GTS3C/TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
     & ,G0,RL,DD,DB14,TR12
      COMMON/MESO7/TN1(5),TN2(4),TN3(5),TGN1(2),TGN2(2),TGN3(2)
      DIMENSION D(9),T(2),MT(11),AP(1),ALTL(8)
      COMMON/LOWER7/PTM(10),PDM(10,8)
      COMMON/PARM7/PT(150),PD(150,9),PS(150),PDL(25,2),PTL(100,4),
     $ PMA(100,10),SAM(100)
      COMMON/CSW/SW(25),ISW,SWC(25)
      COMMON/TTEST/TINFG,GB,ROUT,TT(15)
      COMMON/DMIX/DM04,DM16,DM28,DM32,DM40,DM01,DM14
      COMMON/METSEL/IMR
      SAVE
      DATA MT/48,0,4,16,28,32,40,1,49,14,17/
      DATA ALTL/200.,300.,160.,250.,240.,450.,320.,450./
      DATA MN1/5/,ZN1/120.,110.,100.,90.,72.5/
      DATA DGTR/1.74533E-2/,DR/1.72142E-2/,ALAST/-999./
      DATA ALPHA/-0.38,0.,0.,0.,0.17,0.,-0.38,0.,0./
C        Test for changed input
      V2=VTST7(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,2)
C
      YRD=IYD
      ZA=PDL(16,2)
      ZN1(1)=ZA
      DO 2 J=1,9
        D(J)=0.
    2 CONTINUE
C        TINF VARIATIONS NOT IMPORTANT BELOW ZA OR ZN1(1)
      IF(ALT.GT.ZN1(1)) THEN
        IF(V2.EQ.1..OR.ALAST.LE.ZN1(1)) TINF=PTM(1)*PT(1)
     $  *(1.+SW(16)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PT))
      ELSE
        TINF=PTM(1)*PT(1)
      ENDIF
      T(1)=TINF
C          GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5)
      IF(ALT.GT.ZN1(5)) THEN
        IF(V2.EQ.1.OR.ALAST.LE.ZN1(5)) G0=PTM(4)*PS(1)
     $   *(1.+SW(19)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PS))
      ELSE
        G0=PTM(4)*PS(1)
      ENDIF
C      Calculate these temperatures only if input changed
      IF(V2.EQ.1. .OR. ALT.LT.300.)
     $  TLB=PTM(2)*(1.+SW(17)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,
     $  F107A,F107,AP,PD(1,4)))*PD(1,4)
       S=G0/(TINF-TLB)
C       Lower thermosphere temp variations not significant for
C        density above 300 km
       IF(ALT.LT.300.) THEN
        IF(V2.EQ.1..OR.ALAST.GE.300.) THEN
         TN1(2)=PTM(7)*PTL(1,1)/(1.-SW(18)*GLOB7S(PTL(1,1)))
         TN1(3)=PTM(3)*PTL(1,2)/(1.-SW(18)*GLOB7S(PTL(1,2)))
         TN1(4)=PTM(8)*PTL(1,3)/(1.-SW(18)*GLOB7S(PTL(1,3)))
         TN1(5)=PTM(5)*PTL(1,4)/(1.-SW(18)*SW(20)*GLOB7S(PTL(1,4)))
         TGN1(2)=PTM(9)*PMA(1,9)*(1.+SW(18)*SW(20)*GLOB7S(PMA(1,9)))
     $   *TN1(5)*TN1(5)/(PTM(5)*PTL(1,4))**2
        ENDIF
       ELSE
        TN1(2)=PTM(7)*PTL(1,1)
        TN1(3)=PTM(3)*PTL(1,2)
        TN1(4)=PTM(8)*PTL(1,3)
        TN1(5)=PTM(5)*PTL(1,4)
        TGN1(2)=PTM(9)*PMA(1,9)
     $  *TN1(5)*TN1(5)/(PTM(5)*PTL(1,4))**2
       ENDIF
C
      Z0=ZN1(4)
      T0=TN1(4)
      TR12=1.
C
      IF(MASS.EQ.0) GO TO 50
C       N2 variation factor at Zlb
      G28=SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107, 
     & AP,PD(1,3))
      DAY=AMOD(YRD,1000.)
C        VARIATION OF TURBOPAUSE HEIGHT
      ZHF=PDL(25,2)
     $    *(1.+SW(5)*PDL(25,1)*SIN(DGTR*GLAT)*COS(DR*(DAY-PT(14))))
      YRD=IYD
      T(1)=TINF
      XMM=PDM(5,3)
      Z=ALT
C
      DO 10 J = 1,11
      IF(MASS.EQ.MT(J))   GO TO 15
   10 CONTINUE
      WRITE(6,100) MASS
      GO TO 90
   15 IF(Z.GT.ALTL(6).AND.MASS.NE.28.AND.MASS.NE.48) GO TO 17
C
C       **** N2 DENSITY ****
C
C      Diffusive density at Zlb
      DB28 = PDM(1,3)*EXP(G28)*PD(1,3)
C      Diffusive density at Alt
      D(3)=DENSU(Z,DB28,TINF,TLB, 28.,ALPHA(3),T(2),PTM(6),S,MN1,ZN1,
     & TN1,TGN1)
      DD=D(3)
C      Turbopause
      ZH28=PDM(3,3)*ZHF
      ZHM28=PDM(4,3)*PDL(6,2) 
      XMD=28.-XMM
C      Mixed density at Zlb
      B28=DENSU(ZH28,DB28,TINF,TLB,XMD,ALPHA(3)-1.,TZ,PTM(6),S,MN1,
     & ZN1,TN1,TGN1)
      IF(Z.GT.ALTL(3).OR.SW(15).EQ.0.) GO TO 17
C      Mixed density at Alt
      DM28=DENSU(Z,B28,TINF,TLB,XMM,ALPHA(3),TZ,PTM(6),S,MN1,
     & ZN1,TN1,TGN1)
C      Net density at Alt
      D(3)=DNET(D(3),DM28,ZHM28,XMM,28.)
   17 CONTINUE
      GO TO (20,50,20,25,90,35,40,45,25,48,46),  J
   20 CONTINUE
C
C       **** HE DENSITY ****
C
C       Density variation factor at Zlb
      G4 = SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,1))
C      Diffusive density at Zlb
      DB04 = PDM(1,1)*EXP(G4)*PD(1,1)
C      Diffusive density at Alt
      D(1)=DENSU(Z,DB04,TINF,TLB, 4.,ALPHA(1),T(2),PTM(6),S,MN1,ZN1,
     & TN1,TGN1)
      DD=D(1)
      IF(Z.GT.ALTL(1).OR.SW(15).EQ.0.) GO TO 24
C      Turbopause
      ZH04=PDM(3,1)
C      Mixed density at Zlb
      B04=DENSU(ZH04,DB04,TINF,TLB,4.-XMM,ALPHA(1)-1.,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM04=DENSU(Z,B04,TINF,TLB,XMM,0.,T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
      ZHM04=ZHM28
C      Net density at Alt
      D(1)=DNET(D(1),DM04,ZHM04,XMM,4.)
C      Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,1)/B04)
      ZC04=PDM(5,1)*PDL(1,2)
      HC04=PDM(6,1)*PDL(2,2)
C      Net density corrected at Alt
      D(1)=D(1)*CCOR(Z,RL,HC04,ZC04)
   24 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   25 CONTINUE
C
C      **** O DENSITY ****
C
C       Density variation factor at Zlb
      G16= SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,2))
C      Diffusive density at Zlb
      DB16 =  PDM(1,2)*EXP(G16)*PD(1,2)
C       Diffusive density at Alt
      D(2)=DENSU(Z,DB16,TINF,TLB, 16.,ALPHA(2),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      DD=D(2)
      IF(Z.GT.ALTL(2).OR.SW(15).EQ.0.) GO TO 34
C  Corrected from PDM(3,1) to PDM(3,2)  12/2/85
C       Turbopause
      ZH16=PDM(3,2)
C      Mixed density at Zlb
      B16=DENSU(ZH16,DB16,TINF,TLB,16-XMM,ALPHA(2)-1.,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM16=DENSU(Z,B16,TINF,TLB,XMM,0.,T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
      ZHM16=ZHM28
C      Net density at Alt
      D(2)=DNET(D(2),DM16,ZHM16,XMM,16.)
C   3/16/99 Change form to match O2 departure from diff equil near 150
C   km and add dependence on F10.7
C      RL=ALOG(B28*PDM(2,2)*ABS(PDL(17,2))/B16)
      RL=PDM(2,2)*PDL(17,2)*(1.+SW(1)*PDL(24,1)*(F107A-150.))
      HC16=PDM(6,2)*PDL(4,2)
      ZC16=PDM(5,2)*PDL(3,2)
      HC216=PDM(6,2)*PDL(5,2)
      D(2)=D(2)*CCOR2(Z,RL,HC16,ZC16,HC216)
C       Chemistry correction
      HCC16=PDM(8,2)*PDL(14,2)
      ZCC16=PDM(7,2)*PDL(13,2)
      RC16=PDM(4,2)*PDL(15,2)
C      Net density corrected at Alt
      D(2)=D(2)*CCOR(Z,RC16,HCC16,ZCC16)
   34 CONTINUE
      IF(MASS.NE.48.AND.MASS.NE.49) GO TO 90
   35 CONTINUE
C
C       **** O2 DENSITY ****
C
C       Density variation factor at Zlb
      G32= SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,5))
C      Diffusive density at Zlb
      DB32 = PDM(1,4)*EXP(G32)*PD(1,5)
C       Diffusive density at Alt
      D(4)=DENSU(Z,DB32,TINF,TLB, 32.,ALPHA(4),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      IF(MASS.EQ.49) THEN
         DD=DD+2.*D(4)
      ELSE
         DD=D(4)
      ENDIF
      IF(SW(15).EQ.0.) GO TO 39
      IF(Z.GT.ALTL(4)) GO TO 38
C       Turbopause
      ZH32=PDM(3,4)
C      Mixed density at Zlb
      B32=DENSU(ZH32,DB32,TINF,TLB,32.-XMM,ALPHA(4)-1.,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM32=DENSU(Z,B32,TINF,TLB,XMM,0.,T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
      ZHM32=ZHM28
C      Net density at Alt
      D(4)=DNET(D(4),DM32,ZHM32,XMM,32.)
C       Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,4)/B32)
      HC32=PDM(6,4)*PDL(8,2)
      ZC32=PDM(5,4)*PDL(7,2)
      D(4)=D(4)*CCOR(Z,RL,HC32,ZC32)
   38 CONTINUE
C      Correction for general departure from diffusive equilibrium above Zlb
      HCC32=PDM(8,4)*PDL(23,2)
      HCC232=PDM(8,4)*PDL(23,1)
      ZCC32=PDM(7,4)*PDL(22,2)
      RC32=PDM(4,4)*PDL(24,2)*(1.+SW(1)*PDL(24,1)*(F107A-150.))
C      Net density corrected at Alt
      D(4)=D(4)*CCOR2(Z,RC32,HCC32,ZCC32,HCC232)
   39 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   40 CONTINUE
C
C       **** AR DENSITY ****
C
C       Density variation factor at Zlb
      G40= SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,6))
C      Diffusive density at Zlb
      DB40 = PDM(1,5)*EXP(G40)*PD(1,6)
C       Diffusive density at Alt
      D(5)=DENSU(Z,DB40,TINF,TLB, 40.,ALPHA(5),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      DD=D(5)
      IF(Z.GT.ALTL(5).OR.SW(15).EQ.0.) GO TO 44
C       Turbopause
      ZH40=PDM(3,5)
C      Mixed density at Zlb
      B40=DENSU(ZH40,DB40,TINF,TLB,40.-XMM,ALPHA(5)-1.,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM40=DENSU(Z,B40,TINF,TLB,XMM,0.,T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
      ZHM40=ZHM28
C      Net density at Alt
      D(5)=DNET(D(5),DM40,ZHM40,XMM,40.)
C       Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,5)/B40)
      HC40=PDM(6,5)*PDL(10,2)
      ZC40=PDM(5,5)*PDL(9,2)
C      Net density corrected at Alt
      D(5)=D(5)*CCOR(Z,RL,HC40,ZC40)
   44 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   45 CONTINUE
C
C        **** HYDROGEN DENSITY ****
C
C       Density variation factor at Zlb
      G1 = SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,7))
C      Diffusive density at Zlb
      DB01 = PDM(1,6)*EXP(G1)*PD(1,7)
C       Diffusive density at Alt
      D(7)=DENSU(Z,DB01,TINF,TLB,1.,ALPHA(7),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      DD=D(7)
      IF(Z.GT.ALTL(7).OR.SW(15).EQ.0.) GO TO 47
C       Turbopause
      ZH01=PDM(3,6)
C      Mixed density at Zlb
      B01=DENSU(ZH01,DB01,TINF,TLB,1.-XMM,ALPHA(7)-1.,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM01=DENSU(Z,B01,TINF,TLB,XMM,0.,T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
      ZHM01=ZHM28
C      Net density at Alt
      D(7)=DNET(D(7),DM01,ZHM01,XMM,1.)
C       Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,6)*ABS(PDL(18,2))/B01)
      HC01=PDM(6,6)*PDL(12,2)
      ZC01=PDM(5,6)*PDL(11,2)
      D(7)=D(7)*CCOR(Z,RL,HC01,ZC01)
C       Chemistry correction
      HCC01=PDM(8,6)*PDL(20,2)
      ZCC01=PDM(7,6)*PDL(19,2)
      RC01=PDM(4,6)*PDL(21,2)
C      Net density corrected at Alt
      D(7)=D(7)*CCOR(Z,RC01,HCC01,ZCC01)
   47 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   48 CONTINUE
C
C        **** ATOMIC NITROGEN DENSITY ****
C
C       Density variation factor at Zlb
      G14 = SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,8))
C      Diffusive density at Zlb
      DB14 = PDM(1,7)*EXP(G14)*PD(1,8)
C       Diffusive density at Alt
      D(8)=DENSU(Z,DB14,TINF,TLB,14.,ALPHA(8),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      DD=D(8)
      IF(Z.GT.ALTL(8).OR.SW(15).EQ.0.) GO TO 49
C       Turbopause
      ZH14=PDM(3,7)
C      Mixed density at Zlb
      B14=DENSU(ZH14,DB14,TINF,TLB,14.-XMM,ALPHA(8)-1.,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM14=DENSU(Z,B14,TINF,TLB,XMM,0.,T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
      ZHM14=ZHM28
C      Net density at Alt
      D(8)=DNET(D(8),DM14,ZHM14,XMM,14.)
C       Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,7)*ABS(PDL(3,1))/B14)
      HC14=PDM(6,7)*PDL(2,1)
      ZC14=PDM(5,7)*PDL(1,1)
      D(8)=D(8)*CCOR(Z,RL,HC14,ZC14)
C       Chemistry correction
      HCC14=PDM(8,7)*PDL(5,1)
      ZCC14=PDM(7,7)*PDL(4,1)
      RC14=PDM(4,7)*PDL(6,1)
C      Net density corrected at Alt
      D(8)=D(8)*CCOR(Z,RC14,HCC14,ZCC14)
   49 CONTINUE
      IF(MASS.NE.48) GO TO 90
   46 CONTINUE
C
C        **** Anomalous OXYGEN DENSITY ****
C
      G16H = SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,9))
      DB16H = PDM(1,8)*EXP(G16H)*PD(1,9)
      THO=PDM(10,8)*PDL(7,1)
      DD=DENSU(Z,DB16H,THO,THO,16.,ALPHA(9),T2,PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      ZSHT=PDM(6,8)
      ZMHO=PDM(5,8)
      ZSHO=SCALH(ZMHO,16.,THO)
      D(9)=DD*EXP(-ZSHT/ZSHO*(EXP(-(Z-ZMHO)/ZSHT)-1.))
      IF(MASS.NE.48) GO TO 90
C
C       TOTAL MASS DENSITY
C
      D(6) = 1.66E-24*(4.*D(1)+16.*D(2)+28.*D(3)+32.*D(4)+40.*D(5)+
     &       D(7)+14.*D(8))
      DB48=1.66E-24*(4.*DB04+16.*DB16+28.*DB28+32.*DB32+40.*DB40+DB01+
     &        14.*DB14)
      GO TO 90
C       TEMPERATURE AT ALTITUDE
   50 CONTINUE
      Z=ABS(ALT)
      DDUM  = DENSU(Z,1., TINF,TLB,0.,0.,T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
   90 CONTINUE
C       ADJUST DENSITIES FROM CGS TO KGM
      IF(IMR.EQ.1) THEN
        DO 95 I=1,9
          D(I)=D(I)*1.E6
   95   CONTINUE
        D(6)=D(6)/1000.
      ENDIF
      ALAST=ALT
      RETURN
  100 FORMAT(1X,'MASS', I5, '  NOT VALID')
      END
C-----------------------------------------------------------------------
      SUBROUTINE METERS(METER)
C      Convert outputs to Kg & Meters if METER true
      LOGICAL METER
      COMMON/METSEL/IMR
      SAVE
      IMR=0
      IF(METER) IMR=1
      END
C-----------------------------------------------------------------------
      FUNCTION SCALH(ALT,XM,TEMP)
C      Calculate scale height (km)
      COMMON/PARMB/GSURF,RE
      SAVE
      DATA RGAS/831.4/
      G=GSURF/(1.+ALT/RE)**2
      SCALH=RGAS*TEMP/(G*XM)
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION GLOBE7(YRD,SEC,LAT,LONG,TLOC,F107A,F107,AP,P)
C       CALCULATE G(L) FUNCTION 
C       Upper Thermosphere Parameters
      REAL LAT, LONG
      DIMENSION P(*),SV(25),AP(*)
      COMMON/TTEST/TINF,GB,ROUT,T(15)
      COMMON/CSW/SW(25),ISW,SWC(25)
      COMMON/LPOLY/PLG(9,4),CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC,
     $ IYR,DAY,DF,DFA,APD,APDF,APT(4),XLONG
      SAVE
      DATA DGTR/1.74533E-2/,DR/1.72142E-2/, XL/1000./,TLL/1000./
      DATA SW9/1./,DAYL/-1./,P14/-1000./,P18/-1000./,P32/-1000./
      DATA HR/.2618/,SR/7.2722E-5/,SV/25*1./,NSW/14/,P39/-1000./
C       3hr Magnetic activity functions
C      Eq. A24d
      G0(A)=(A-4.+(P(26)-1.)*(A-4.+(EXP(-ABS(P(25))*(A-4.))-1.)/ABS(P(25
     *))))
C       Eq. A24c
       SUMEX(EX)=1.+(1.-EX**19)/(1.-EX)*EX**(.5)
C       Eq. A24a
      SG0(EX)=(G0(AP(2))+(G0(AP(3))*EX+G0(AP(4))*EX*EX+G0(AP(5))*EX**3
     $ +(G0(AP(6))*EX**4+G0(AP(7))*EX**12)*(1.-EX**8)/(1.-EX))
     $ )/SUMEX(EX)
      IF(ISW.NE.64999) CALL TSELEC(SV)
      DO 10 J=1,14
       T(J)=0
   10 CONTINUE
      IF(SW(9).GT.0) SW9=1.
      IF(SW(9).LT.0) SW9=-1.
      IYR = YRD/1000.
      DAY = YRD - IYR*1000.
      XLONG=LONG
C      Eq. A22 (remainder of code)
      IF(XL.EQ.LAT)   GO TO 15
C          CALCULATE LEGENDRE POLYNOMIALS
      C = SIN(LAT*DGTR)
      S = COS(LAT*DGTR)
      C2 = C*C
      C4 = C2*C2
      S2 = S*S
      PLG(2,1) = C
      PLG(3,1) = 0.5*(3.*C2 -1.)
      PLG(4,1) = 0.5*(5.*C*C2-3.*C)
      PLG(5,1) = (35.*C4 - 30.*C2 + 3.)/8.
      PLG(6,1) = (63.*C2*C2*C - 70.*C2*C + 15.*C)/8.
      PLG(7,1) = (11.*C*PLG(6,1) - 5.*PLG(5,1))/6.
C     PLG(8,1) = (13.*C*PLG(7,1) - 6.*PLG(6,1))/7.
      PLG(2,2) = S
      PLG(3,2) = 3.*C*S
      PLG(4,2) = 1.5*(5.*C2-1.)*S
      PLG(5,2) = 2.5*(7.*C2*C-3.*C)*S
      PLG(6,2) = 1.875*(21.*C4 - 14.*C2 +1.)*S
      PLG(7,2) = (11.*C*PLG(6,2)-6.*PLG(5,2))/5.
C     PLG(8,2) = (13.*C*PLG(7,2)-7.*PLG(6,2))/6.
C     PLG(9,2) = (15.*C*PLG(8,2)-8.*PLG(7,2))/7.
      PLG(3,3) = 3.*S2
      PLG(4,3) = 15.*S2*C
      PLG(5,3) = 7.5*(7.*C2 -1.)*S2
      PLG(6,3) = 3.*C*PLG(5,3)-2.*PLG(4,3)
      PLG(7,3)=(11.*C*PLG(6,3)-7.*PLG(5,3))/4.
      PLG(8,3)=(13.*C*PLG(7,3)-8.*PLG(6,3))/5.
      PLG(4,4) = 15.*S2*S
      PLG(5,4) = 105.*S2*S*C 
      PLG(6,4)=(9.*C*PLG(5,4)-7.*PLG(4,4))/2.
      PLG(7,4)=(11.*C*PLG(6,4)-8.*PLG(5,4))/3.
      XL=LAT
   15 CONTINUE
      IF(TLL.EQ.TLOC)   GO TO 16
      IF(SW(7).EQ.0.AND.SW(8).EQ.0.AND.SW(14).EQ.0) GOTO 16
      STLOC = SIN(HR*TLOC)
      CTLOC = COS(HR*TLOC)
      S2TLOC = SIN(2.*HR*TLOC)
      C2TLOC = COS(2.*HR*TLOC)
      S3TLOC = SIN(3.*HR*TLOC)
      C3TLOC = COS(3.*HR*TLOC)
      TLL = TLOC
   16 CONTINUE
      IF(DAY.NE.DAYL.OR.P(14).NE.P14) CD14=COS(DR*(DAY-P(14)))
      IF(DAY.NE.DAYL.OR.P(18).NE.P18) CD18=COS(2.*DR*(DAY-P(18)))
      IF(DAY.NE.DAYL.OR.P(32).NE.P32) CD32=COS(DR*(DAY-P(32)))
      IF(DAY.NE.DAYL.OR.P(39).NE.P39) CD39=COS(2.*DR*(DAY-P(39)))
      DAYL = DAY
      P14 = P(14)
      P18 = P(18)
      P32 = P(32)
      P39 = P(39)
C         F10.7 EFFECT
      DF = F107 - F107A
      DFA=F107A-150.
      T(1) =  P(20)*DF*(1.+P(60)*DFA) + P(21)*DF*DF + P(22)*DFA
     & + P(30)*DFA**2
      F1 = 1. + (P(48)*DFA +P(20)*DF+P(21)*DF*DF)*SWC(1)
      F2 = 1. + (P(50)*DFA+P(20)*DF+P(21)*DF*DF)*SWC(1)
C        TIME INDEPENDENT
      T(2) =
     1  (P(2)*PLG(3,1) + P(3)*PLG(5,1)+P(23)*PLG(7,1))
     & +(P(15)*PLG(3,1))*DFA*SWC(1)
     2 +P(27)*PLG(2,1)
C        SYMMETRICAL ANNUAL
      T(3) =
     1 (P(19) )*CD32
C        SYMMETRICAL SEMIANNUAL
      T(4) =
     1 (P(16)+P(17)*PLG(3,1))*CD18
C        ASYMMETRICAL ANNUAL
      T(5) =  F1*
     1  (P(10)*PLG(2,1)+P(11)*PLG(4,1))*CD14
C         ASYMMETRICAL SEMIANNUAL
      T(6) =    P(38)*PLG(2,1)*CD39
C        DIURNAL
      IF(SW(7).EQ.0) GOTO 200
      T71 = (P(12)*PLG(3,2))*CD14*SWC(5)
      T72 = (P(13)*PLG(3,2))*CD14*SWC(5)
      T(7) = F2*
     1 ((P(4)*PLG(2,2) + P(5)*PLG(4,2) + P(28)*PLG(6,2)
     2 + T71)*CTLOC
     4 + (P(7)*PLG(2,2) + P(8)*PLG(4,2) +P(29)*PLG(6,2)
     5 + T72)*STLOC)
  200 CONTINUE
C        SEMIDIURNAL
      IF(SW(8).EQ.0) GOTO 210
      T81 = (P(24)*PLG(4,3)+P(36)*PLG(6,3))*CD14*SWC(5) 
      T82 = (P(34)*PLG(4,3)+P(37)*PLG(6,3))*CD14*SWC(5)
      T(8) = F2*
     1 ((P(6)*PLG(3,3) + P(42)*PLG(5,3) + T81)*C2TLOC
     3 +(P(9)*PLG(3,3) + P(43)*PLG(5,3) + T82)*S2TLOC)
  210 CONTINUE
C        TERDIURNAL
      IF(SW(14).EQ.0) GOTO 220
      T(14) = F2*
     1 ((P(40)*PLG(4,4)+(P(94)*PLG(5,4)+P(47)*PLG(7,4))*CD14*SWC(5))*
     $ S3TLOC
     2 +(P(41)*PLG(4,4)+(P(95)*PLG(5,4)+P(49)*PLG(7,4))*CD14*SWC(5))*
     $ C3TLOC)
  220 CONTINUE
C          MAGNETIC ACTIVITY BASED ON DAILY AP

      IF(SW9.EQ.-1.) GO TO 30
      APD=(AP(1)-4.)
      P44=P(44)
      P45=P(45)
      IF(P44.LT.0) P44=1.E-5
      APDF = APD+(P45-1.)*(APD+(EXP(-P44  *APD)-1.)/P44)
      IF(SW(9).EQ.0) GOTO 40
      T(9)=APDF*(P(33)+P(46)*PLG(3,1)+P(35)*PLG(5,1)+
     $ (P(101)*PLG(2,1)+P(102)*PLG(4,1)+P(103)*PLG(6,1))*CD14*SWC(5)+
     $ (P(122)*PLG(2,2)+P(123)*PLG(4,2)+P(124)*PLG(6,2))*SWC(7)*
     $ COS(HR*(TLOC-P(125))))
      GO TO 40
   30 CONTINUE
      IF(P(52).EQ.0) GO TO 40
      EXP1 = EXP(-10800.*ABS(P(52))/(1.+P(139)*(45.-ABS(LAT))))
      IF(EXP1.GT..99999) EXP1=.99999
      IF(P(25).LT.1.E-4) P(25)=1.E-4
      APT(1)=SG0(EXP1)
C      APT(2)=SG2(EXP1)
c      APT(3)=SG0(EXP2)
C      APT(4)=SG2(EXP2)
      IF(SW(9).EQ.0) GOTO 40
      T(9) = APT(1)*(P(51)+P(97)*PLG(3,1)+P(55)*PLG(5,1)+
     $ (P(126)*PLG(2,1)+P(127)*PLG(4,1)+P(128)*PLG(6,1))*CD14*SWC(5)+
     $ (P(129)*PLG(2,2)+P(130)*PLG(4,2)+P(131)*PLG(6,2))*SWC(7)*
     $ COS(HR*(TLOC-P(132))))
  40  CONTINUE
      IF(SW(10).EQ.0.OR.LONG.LE.-1000.) GO TO 49
C        LONGITUDINAL
      IF(SW(11).EQ.0) GOTO 230
      T(11)= (1.+P(81)*DFA*SWC(1))*
     $((P(65)*PLG(3,2)+P(66)*PLG(5,2)+P(67)*PLG(7,2)
     $ +P(104)*PLG(2,2)+P(105)*PLG(4,2)+P(106)*PLG(6,2)
     $ +SWC(5)*(P(110)*PLG(2,2)+P(111)*PLG(4,2)+P(112)*PLG(6,2))*CD14)*
     $     COS(DGTR*LONG)
     $ +(P(91)*PLG(3,2)+P(92)*PLG(5,2)+P(93)*PLG(7,2)
     $ +P(107)*PLG(2,2)+P(108)*PLG(4,2)+P(109)*PLG(6,2)
     $ +SWC(5)*(P(113)*PLG(2,2)+P(114)*PLG(4,2)+P(115)*PLG(6,2))*CD14)*
     $  SIN(DGTR*LONG))
  230 CONTINUE
C        UT AND MIXED UT,LONGITUDE
      IF(SW(12).EQ.0) GOTO 240
      T(12)=(1.+P(96)*PLG(2,1))*(1.+P(82)*DFA*SWC(1))*
     $(1.+P(120)*PLG(2,1)*SWC(5)*CD14)*
     $((P(69)*PLG(2,1)+P(70)*PLG(4,1)+P(71)*PLG(6,1))*
     $     COS(SR*(SEC-P(72))))
      T(12)=T(12)+SWC(11)*
     $ (P(77)*PLG(4,3)+P(78)*PLG(6,3)+P(79)*PLG(8,3))*
     $     COS(SR*(SEC-P(80))+2.*DGTR*LONG)*(1.+P(138)*DFA*SWC(1))
  240 CONTINUE
C        UT,LONGITUDE MAGNETIC ACTIVITY
      IF(SW(13).EQ.0) GOTO 48
      IF(SW9.EQ.-1.) GO TO 45
      T(13)= APDF*SWC(11)*(1.+P(121)*PLG(2,1))*
     $((P( 61)*PLG(3,2)+P( 62)*PLG(5,2)+P( 63)*PLG(7,2))*
     $     COS(DGTR*(LONG-P( 64))))
     $ +APDF*SWC(11)*SWC(5)*
     $ (P(116)*PLG(2,2)+P(117)*PLG(4,2)+P(118)*PLG(6,2))*
     $     CD14*COS(DGTR*(LONG-P(119)))
     $ + APDF*SWC(12)*
     $ (P( 84)*PLG(2,1)+P( 85)*PLG(4,1)+P( 86)*PLG(6,1))*
     $     COS(SR*(SEC-P( 76)))
      GOTO 48
   45 CONTINUE
      IF(P(52).EQ.0) GOTO 48
      T(13)=APT(1)*SWC(11)*(1.+P(133)*PLG(2,1))*
     $((P(53)*PLG(3,2)+P(99)*PLG(5,2)+P(68)*PLG(7,2))*
     $     COS(DGTR*(LONG-P(98))))
     $ +APT(1)*SWC(11)*SWC(5)*
     $ (P(134)*PLG(2,2)+P(135)*PLG(4,2)+P(136)*PLG(6,2))*
     $     CD14*COS(DGTR*(LONG-P(137)))
     $ +APT(1)*SWC(12)*
     $ (P(56)*PLG(2,1)+P(57)*PLG(4,1)+P(58)*PLG(6,1))*
     $     COS(SR*(SEC-P(59)))
   48 CONTINUE
C  PARMS NOT USED: 83, 90,100,140-150
   49 CONTINUE
      TINF=P(31)
      DO 50 I = 1,NSW
   50 TINF = TINF + ABS(SW(I))*T(I)
      GLOBE7 = TINF
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION GLOB7S(P)
C      VERSION OF GLOBE FOR LOWER ATMOSPHERE 10/26/99
      REAL LONG
      COMMON/LPOLY/PLG(9,4),CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC,
     $ IYR,DAY,DF,DFA,APD,APDF,APT(4),LONG
      COMMON/CSW/SW(25),ISW,SWC(25)
      DIMENSION P(*),T(14)
      SAVE
      DATA DR/1.72142E-2/,DGTR/1.74533E-2/,PSET/2./
      DATA DAYL/-1./,P32,P18,P14,P39/4*-1000./
C       CONFIRM PARAMETER SET
      IF(P(100).EQ.0) P(100)=PSET
      IF(P(100).NE.PSET) THEN
        WRITE(6,900) PSET,P(100)
  900   FORMAT(1X,'WRONG PARAMETER SET FOR GLOB7S',3F10.1)
        STOP
      ENDIF
      DO 10 J=1,14
        T(J)=0.
   10 CONTINUE
      IF(DAY.NE.DAYL.OR.P32.NE.P(32)) CD32=COS(DR*(DAY-P(32)))
      IF(DAY.NE.DAYL.OR.P18.NE.P(18)) CD18=COS(2.*DR*(DAY-P(18)))       
      IF(DAY.NE.DAYL.OR.P14.NE.P(14)) CD14=COS(DR*(DAY-P(14)))
      IF(DAY.NE.DAYL.OR.P39.NE.P(39)) CD39=COS(2.*DR*(DAY-P(39)))
      DAYL=DAY
      P32=P(32)
      P18=P(18)
      P14=P(14)
      P39=P(39)
C
C       F10.7
      T(1)=P(22)*DFA
C       TIME INDEPENDENT
      T(2)=P(2)*PLG(3,1)+P(3)*PLG(5,1)+P(23)*PLG(7,1)
     $     +P(27)*PLG(2,1)+P(15)*PLG(4,1)+P(60)*PLG(6,1)
C       SYMMETRICAL ANNUAL
      T(3)=(P(19)+P(48)*PLG(3,1)+P(30)*PLG(5,1))*CD32
C       SYMMETRICAL SEMIANNUAL
      T(4)=(P(16)+P(17)*PLG(3,1)+P(31)*PLG(5,1))*CD18
C       ASYMMETRICAL ANNUAL
      T(5)=(P(10)*PLG(2,1)+P(11)*PLG(4,1)+P(21)*PLG(6,1))*CD14
C       ASYMMETRICAL SEMIANNUAL
      T(6)=(P(38)*PLG(2,1))*CD39
C        DIURNAL
      IF(SW(7).EQ.0) GOTO 200
      T71 = P(12)*PLG(3,2)*CD14*SWC(5)
      T72 = P(13)*PLG(3,2)*CD14*SWC(5)
      T(7) = 
     1 ((P(4)*PLG(2,2) + P(5)*PLG(4,2)
     2 + T71)*CTLOC
     4 + (P(7)*PLG(2,2) + P(8)*PLG(4,2)
     5 + T72)*STLOC)
  200 CONTINUE
C        SEMIDIURNAL
      IF(SW(8).EQ.0) GOTO 210
      T81 = (P(24)*PLG(4,3)+P(36)*PLG(6,3))*CD14*SWC(5) 
      T82 = (P(34)*PLG(4,3)+P(37)*PLG(6,3))*CD14*SWC(5)
      T(8) = 
     1 ((P(6)*PLG(3,3) + P(42)*PLG(5,3) + T81)*C2TLOC
     3 +(P(9)*PLG(3,3) + P(43)*PLG(5,3) + T82)*S2TLOC)
  210 CONTINUE
C        TERDIURNAL
      IF(SW(14).EQ.0) GOTO 220
      T(14) = P(40)*PLG(4,4)*S3TLOC
     $ +P(41)*PLG(4,4)*C3TLOC
  220 CONTINUE
C       MAGNETIC ACTIVITY
      IF(SW(9).EQ.0) GOTO 40
      IF(SW(9).EQ.1)
     $ T(9)=APDF*(P(33)+P(46)*PLG(3,1)*SWC(2))
      IF(SW(9).EQ.-1)
     $ T(9)=(P(51)*APT(1)+P(97)*PLG(3,1)*APT(1)*SWC(2))
   40 CONTINUE
      IF(SW(10).EQ.0.OR.SW(11).EQ.0.OR.LONG.LE.-1000.) GO TO 49
C        LONGITUDINAL
      T(11)= (1.+PLG(2,1)*(P(81)*SWC(5)*COS(DR*(DAY-P(82)))
     $           +P(86)*SWC(6)*COS(2.*DR*(DAY-P(87))))
     $        +P(84)*SWC(3)*COS(DR*(DAY-P(85)))
     $           +P(88)*SWC(4)*COS(2.*DR*(DAY-P(89))))
     $ *((P(65)*PLG(3,2)+P(66)*PLG(5,2)+P(67)*PLG(7,2)
     $   +P(75)*PLG(2,2)+P(76)*PLG(4,2)+P(77)*PLG(6,2)
     $    )*COS(DGTR*LONG)
     $  +(P(91)*PLG(3,2)+P(92)*PLG(5,2)+P(93)*PLG(7,2)
     $   +P(78)*PLG(2,2)+P(79)*PLG(4,2)+P(80)*PLG(6,2)
     $    )*SIN(DGTR*LONG))
   49 CONTINUE
      TT=0.
      DO 50 I=1,14
   50 TT=TT+ABS(SW(I))*T(I)
      GLOB7S=TT
      RETURN
      END
C--------------------------------------------------------------------
      FUNCTION DENSU(ALT,DLB,TINF,TLB,XM,ALPHA,TZ,ZLB,S2,
     $  MN1,ZN1,TN1,TGN1)
C       Calculate Temperature and Density Profiles for MSIS models
C       New lower thermo polynomial 10/30/89
      DIMENSION ZN1(MN1),TN1(MN1),TGN1(2),XS(5),YS(5),Y2OUT(5)
      COMMON/PARMB/GSURF,RE
      COMMON/LSQV/MP,II,JG,LT,QPB(50),IERR,IFUN,N,J,DV(60)
      SAVE
      DATA RGAS/831.4/
      ZETA(ZZ,ZL)=(ZZ-ZL)*(RE+ZL)/(RE+ZZ)
CCCCCCWRITE(6,*) 'DB',ALT,DLB,TINF,TLB,XM,ALPHA,ZLB,S2,MN1,ZN1,TN1
      DENSU=1.
C        Joining altitude of Bates and spline
      ZA=ZN1(1)
      Z=AMAX1(ALT,ZA)
C      Geopotential altitude difference from ZLB
      ZG2=ZETA(Z,ZLB)
C      Bates temperature
      TT=TINF-(TINF-TLB)*EXP(-S2*ZG2)
      TA=TT
      TZ=TT
      DENSU=TZ
      IF(ALT.GE.ZA) GO TO 10
C
C       CALCULATE TEMPERATURE BELOW ZA
C      Temperature gradient at ZA from Bates profile
      DTA=(TINF-TA)*S2*((RE+ZLB)/(RE+ZA))**2
      TGN1(1)=DTA 
      TN1(1)=TA
      Z=AMAX1(ALT,ZN1(MN1))
      MN=MN1
      Z1=ZN1(1)
      Z2=ZN1(MN)
      T1=TN1(1)
      T2=TN1(MN)
C      Geopotental difference from Z1
      ZG=ZETA(Z,Z1)
      ZGDIF=ZETA(Z2,Z1)
C       Set up spline nodes
      DO 20 K=1,MN
        XS(K)=ZETA(ZN1(K),Z1)/ZGDIF
        YS(K)=1./TN1(K)
   20 CONTINUE
C        End node derivatives
      YD1=-TGN1(1)/(T1*T1)*ZGDIF
      YD2=-TGN1(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2
C       Calculate spline coefficients
      CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT(XS,YS,Y2OUT,MN,X,Y)
C       temperature at altitude
      TZ=1./Y
      DENSU=TZ
   10 IF(XM.EQ.0.) GO TO 50
C
C      CALCULATE DENSITY ABOVE ZA
      GLB=GSURF/(1.+ZLB/RE)**2
      GAMMA=XM*GLB/(S2*RGAS*TINF)
      EXPL=EXP(-S2*GAMMA*ZG2)
      IF(EXPL.GT.50.OR.TT.LE.0.) THEN
        EXPL=50.
      ENDIF
C       Density at altitude
      DENSA=DLB*(TLB/TT)**(1.+ALPHA+GAMMA)*EXPL
      DENSU=DENSA
      IF(ALT.GE.ZA) GO TO 50
C
C      CALCULATE DENSITY BELOW ZA
      GLB=GSURF/(1.+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C       integrate spline temperatures
      CALL SPLINI(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50..OR.TZ.LE.0.) THEN
        EXPL=50.
      ENDIF
C       Density at altitude
      DENSU=DENSU*(T1/TZ)**(1.+ALPHA)*EXP(-EXPL)
   50 CONTINUE
      RETURN
      END
C--------------------------------------------------------------------
      FUNCTION DENSM(ALT,D0,XM,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)
C       Calculate Temperature and Density Profiles for lower atmos.
      DIMENSION ZN3(MN3),TN3(MN3),TGN3(2),XS(10),YS(10),Y2OUT(10)
      DIMENSION ZN2(MN2),TN2(MN2),TGN2(2)
      COMMON/PARMB/GSURF,RE
      COMMON/FIT/TAF
      COMMON/LSQV/MP,II,JG,LT,QPB(50),IERR,IFUN,N,J,DV(60)
      SAVE
      DATA RGAS/831.4/
      ZETA(ZZ,ZL)=(ZZ-ZL)*(RE+ZL)/(RE+ZZ)
      DENSM=D0
      IF(ALT.GT.ZN2(1)) GOTO 50
C      STRATOSPHERE/MESOSPHERE TEMPERATURE
      Z=AMAX1(ALT,ZN2(MN2))
      MN=MN2
      Z1=ZN2(1)
      Z2=ZN2(MN)
      T1=TN2(1)
      T2=TN2(MN)
      ZG=ZETA(Z,Z1)
      ZGDIF=ZETA(Z2,Z1)
C       Set up spline nodes
      DO 210 K=1,MN
        XS(K)=ZETA(ZN2(K),Z1)/ZGDIF
        YS(K)=1./TN2(K)
  210 CONTINUE
      YD1=-TGN2(1)/(T1*T1)*ZGDIF
      YD2=-TGN2(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2
C       Calculate spline coefficients
      CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT(XS,YS,Y2OUT,MN,X,Y)
C       Temperature at altitude
      TZ=1./Y
      IF(XM.EQ.0.) GO TO 20
C
C      CALCULATE STRATOSPHERE/MESOSPHERE DENSITY 
      GLB=GSURF/(1.+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C       Integrate temperature profile
      CALL SPLINI(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50.) EXPL=50.
C       Density at altitude
      DENSM=DENSM*(T1/TZ)*EXP(-EXPL)
   20 CONTINUE
      IF(ALT.GT.ZN3(1)) GOTO 50
C
C      TROPOSPHERE/STRATOSPHERE TEMPERATURE
      Z=ALT
      MN=MN3
      Z1=ZN3(1)
      Z2=ZN3(MN)
      T1=TN3(1)
      T2=TN3(MN)
      ZG=ZETA(Z,Z1)
      ZGDIF=ZETA(Z2,Z1)
C       Set up spline nodes
      DO 220 K=1,MN
        XS(K)=ZETA(ZN3(K),Z1)/ZGDIF
        YS(K)=1./TN3(K)
  220 CONTINUE
      YD1=-TGN3(1)/(T1*T1)*ZGDIF
      YD2=-TGN3(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2
C       Calculate spline coefficients
      CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT(XS,YS,Y2OUT,MN,X,Y)
C       temperature at altitude
      TZ=1./Y
      IF(XM.EQ.0.) GO TO 30
C
C      CALCULATE TROPOSPHERIC/STRATOSPHERE DENSITY 
C     
      GLB=GSURF/(1.+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C        Integrate temperature profile
      CALL SPLINI(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50.) EXPL=50.
C        Density at altitude
      DENSM=DENSM*(T1/TZ)*EXP(-EXPL)
   30 CONTINUE
   50 CONTINUE
      IF(XM.EQ.0) DENSM=TZ
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SPLINI(XA,YA,Y2A,N,X,YI)
C       INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
C        XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        Y2A: ARRAY OF SECOND DERIVATIVES
C        N: SIZE OF ARRAYS XA,YA,Y2A
C        X: ABSCISSA ENDPOINT FOR INTEGRATION
C        Y: OUTPUT VALUE
      DIMENSION XA(N),YA(N),Y2A(N)
      SAVE
      YI=0
      KLO=1
      KHI=2
    1 CONTINUE
      IF(X.GT.XA(KLO).AND.KHI.LE.N) THEN
        XX=X
        IF(KHI.LT.N) XX=AMIN1(X,XA(KHI))
        H=XA(KHI)-XA(KLO)
        A=(XA(KHI)-XX)/H
        B=(XX-XA(KLO))/H
        A2=A*A
        B2=B*B
        YI=YI+((1.-A2)*YA(KLO)/2.+B2*YA(KHI)/2.+
     $     ((-(1.+A2*A2)/4.+A2/2.)*Y2A(KLO)+
     $     (B2*B2/4.-B2/2.)*Y2A(KHI))*H*H/6.)*H
        KLO=KLO+1
        KHI=KHI+1
        GOTO 1
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION DNET(DD,DM,ZHM,XMM,XM)
C       TURBOPAUSE CORRECTION FOR MSIS MODELS
C         Root mean density
C       8/20/80
C          DD - diffusive density
C          DM - full mixed density
C          ZHM - transition scale length
C          XMM - full mixed molecular weight
C          XM  - species molecular weight
C          DNET - combined density
      SAVE
      A=ZHM/(XMM-XM)
      IF(DM.GT.0.AND.DD.GT.0) GOTO 5
        WRITE(6,*) 'DNET LOG ERROR',DM,DD,XM
        IF(DD.EQ.0.AND.DM.EQ.0) DD=1.
        IF(DM.EQ.0) GOTO 10
        IF(DD.EQ.0) GOTO 20
    5 CONTINUE
      YLOG=A*ALOG(DM/DD)
      IF(YLOG.LT.-10.) GO TO 10
      IF(YLOG.GT.10.)  GO TO 20
        DNET=DD*(1.+EXP(YLOG))**(1/A)
        GO TO 50
   10 CONTINUE
        DNET=DD
        GO TO 50
   20 CONTINUE
        DNET=DM
        GO TO 50
   50 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION  CCOR(ALT, R,H1,ZH)
C        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
C        ALT - altitude
C        R - target ratio
C        H1 - transition scale length
C        ZH - altitude of 1/2 R
      SAVE
      E=(ALT-ZH)/H1
      IF(E.GT.70.) GO TO 20
      IF(E.LT.-70.) GO TO 10
        EX=EXP(E)
        CCOR=R/(1.+EX)
        GO TO 50
   10   CCOR=R
        GO TO 50
   20   CCOR=0.
        GO TO 50
   50 CONTINUE
      CCOR=EXP(CCOR)
       RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION  CCOR2(ALT, R,H1,ZH,H2)
C       O&O2 CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
      E1=(ALT-ZH)/H1
      E2=(ALT-ZH)/H2
      IF(E1.GT.70. .OR. E2.GT.70.) GO TO 20
      IF(E1.LT.-70. .AND. E2.LT.-70) GO TO 10
        EX1=EXP(E1)
        EX2=EXP(E2)
        CCOR2=R/(1.+.5*(EX1+EX2))
        GO TO 50
   10   CCOR2=R
        GO TO 50
   20   CCOR2=0.
        GO TO 50
   50 CONTINUE
      CCOR2=EXP(CCOR2)
       RETURN
      END
C-----------------------------------------------------------------------
      BLOCK DATA GTD7BK
C          MSISE-00 01-FEB-02   
      COMMON/PARM7/PT1(50),PT2(50),PT3(50),PA1(50),PA2(50),PA3(50),
     $ PB1(50),PB2(50),PB3(50),PC1(50),PC2(50),PC3(50),
     $ PD1(50),PD2(50),PD3(50),PE1(50),PE2(50),PE3(50),
     $ PF1(50),PF2(50),PF3(50),PG1(50),PG2(50),PG3(50),
     $ PH1(50),PH2(50),PH3(50),PI1(50),PI2(50),PI3(50),
     $ PJ1(50),PJ2(50),PJ3(50),PK1(50),PL1(50),PL2(50),
     $ PM1(50),PM2(50),PN1(50),PN2(50),PO1(50),PO2(50),
     $ PP1(50),PP2(50),PQ1(50),PQ2(50),PR1(50),PR2(50),
     $ PS1(50),PS2(50),PU1(50),PU2(50),PV1(50),PV2(50),
     $ PW1(50),PW2(50),PX1(50),PX2(50),PY1(50),PY2(50),
     $ PZ1(50),PZ2(50),PAA1(50),PAA2(50)
      COMMON/LOWER7/PTM(10),PDM(10,8)
      COMMON/MAVG7/PAVGM(10)
      COMMON/DATIM7/ISDATE(3),ISTIME(2),NAME(2)
      COMMON/METSEL/IMR
      character*4 ISDATE,ISTIME,NAME
      DATA IMR/0/
      DATA ISDATE/'01-F','EB-0','2   '/,ISTIME/'15:4','9:27'/
      DATA NAME/'MSIS','E-00'/
C         TEMPERATURE
      DATA PT1/
     *  9.86573E-01, 1.62228E-02, 1.55270E-02,-1.04323E-01,-3.75801E-03,
     * -1.18538E-03,-1.24043E-01, 4.56820E-03, 8.76018E-03,-1.36235E-01,
     * -3.52427E-02, 8.84181E-03,-5.92127E-03,-8.61650E+00, 0.00000E+00,
     *  1.28492E-02, 0.00000E+00, 1.30096E+02, 1.04567E-02, 1.65686E-03,
     * -5.53887E-06, 2.97810E-03, 0.00000E+00, 5.13122E-03, 8.66784E-02,
     *  1.58727E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,-7.27026E-06,
     *  0.00000E+00, 6.74494E+00, 4.93933E-03, 2.21656E-03, 2.50802E-03,
     *  0.00000E+00, 0.00000E+00,-2.08841E-02,-1.79873E+00, 1.45103E-03,
     *  2.81769E-04,-1.44703E-03,-5.16394E-05, 8.47001E-02, 1.70147E-01,
     *  5.72562E-03, 5.07493E-05, 4.36148E-03, 1.17863E-04, 4.74364E-03/
      DATA PT2/
     *  6.61278E-03, 4.34292E-05, 1.44373E-03, 2.41470E-05, 2.84426E-03,
     *  8.56560E-04, 2.04028E-03, 0.00000E+00,-3.15994E+03,-2.46423E-03,
     *  1.13843E-03, 4.20512E-04, 0.00000E+00,-9.77214E+01, 6.77794E-03,
     *  5.27499E-03, 1.14936E-03, 0.00000E+00,-6.61311E-03,-1.84255E-02,
     * -1.96259E-02, 2.98618E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  6.44574E+02, 8.84668E-04, 5.05066E-04, 0.00000E+00, 4.02881E+03,
     * -1.89503E-03, 0.00000E+00, 0.00000E+00, 8.21407E-04, 2.06780E-03,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -1.20410E-02,-3.63963E-03, 9.92070E-05,-1.15284E-04,-6.33059E-05,
     * -6.05545E-01, 8.34218E-03,-9.13036E+01, 3.71042E-04, 0.00000E+00/
      DATA PT3/
     *  4.19000E-04, 2.70928E-03, 3.31507E-03,-4.44508E-03,-4.96334E-03,
     * -1.60449E-03, 3.95119E-03, 2.48924E-03, 5.09815E-04, 4.05302E-03,
     *  2.24076E-03, 0.00000E+00, 6.84256E-03, 4.66354E-04, 0.00000E+00,
     * -3.68328E-04, 0.00000E+00, 0.00000E+00,-1.46870E+02, 0.00000E+00,
     *  0.00000E+00, 1.09501E-03, 4.65156E-04, 5.62583E-04, 3.21596E+00,
     *  6.43168E-04, 3.14860E-03, 3.40738E-03, 1.78481E-03, 9.62532E-04,
     *  5.58171E-04, 3.43731E+00,-2.33195E-01, 5.10289E-04, 0.00000E+00,
     *  0.00000E+00,-9.25347E+04, 0.00000E+00,-1.99639E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         HE DENSITY
      DATA PA1/
     *  1.09979E+00,-4.88060E-02,-1.97501E-01,-9.10280E-02,-6.96558E-03,
     *  2.42136E-02, 3.91333E-01,-7.20068E-03,-3.22718E-02, 1.41508E+00,
     *  1.68194E-01, 1.85282E-02, 1.09384E-01,-7.24282E+00, 0.00000E+00,
     *  2.96377E-01,-4.97210E-02, 1.04114E+02,-8.61108E-02,-7.29177E-04,
     *  1.48998E-06, 1.08629E-03, 0.00000E+00, 0.00000E+00, 8.31090E-02,
     *  1.12818E-01,-5.75005E-02,-1.29919E-02,-1.78849E-02,-2.86343E-06,
     *  0.00000E+00,-1.51187E+02,-6.65902E-03, 0.00000E+00,-2.02069E-03,
     *  0.00000E+00, 0.00000E+00, 4.32264E-02,-2.80444E+01,-3.26789E-03,
     *  2.47461E-03, 0.00000E+00, 0.00000E+00, 9.82100E-02, 1.22714E-01,
     * -3.96450E-02, 0.00000E+00,-2.76489E-03, 0.00000E+00, 1.87723E-03/
      DATA PA2/
     * -8.09813E-03, 4.34428E-05,-7.70932E-03, 0.00000E+00,-2.28894E-03,
     * -5.69070E-03,-5.22193E-03, 6.00692E-03,-7.80434E+03,-3.48336E-03,
     * -6.38362E-03,-1.82190E-03, 0.00000E+00,-7.58976E+01,-2.17875E-02,
     * -1.72524E-02,-9.06287E-03, 0.00000E+00, 2.44725E-02, 8.66040E-02,
     *  1.05712E-01, 3.02543E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -6.01364E+03,-5.64668E-03,-2.54157E-03, 0.00000E+00, 3.15611E+02,
     * -5.69158E-03, 0.00000E+00, 0.00000E+00,-4.47216E-03,-4.49523E-03,
     *  4.64428E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  4.51236E-02, 2.46520E-02, 6.17794E-03, 0.00000E+00, 0.00000E+00,
     * -3.62944E-01,-4.80022E-02,-7.57230E+01,-1.99656E-03, 0.00000E+00/
      DATA PA3/
     * -5.18780E-03,-1.73990E-02,-9.03485E-03, 7.48465E-03, 1.53267E-02,
     *  1.06296E-02, 1.18655E-02, 2.55569E-03, 1.69020E-03, 3.51936E-02,
     * -1.81242E-02, 0.00000E+00,-1.00529E-01,-5.10574E-03, 0.00000E+00,
     *  2.10228E-03, 0.00000E+00, 0.00000E+00,-1.73255E+02, 5.07833E-01,
     * -2.41408E-01, 8.75414E-03, 2.77527E-03,-8.90353E-05,-5.25148E+00,
     * -5.83899E-03,-2.09122E-02,-9.63530E-03, 9.77164E-03, 4.07051E-03,
     *  2.53555E-04,-5.52875E+00,-3.55993E-01,-2.49231E-03, 0.00000E+00,
     *  0.00000E+00, 2.86026E+01, 0.00000E+00, 3.42722E-04, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         O DENSITY
      DATA PB1/
     *  1.02315E+00,-1.59710E-01,-1.06630E-01,-1.77074E-02,-4.42726E-03,
     *  3.44803E-02, 4.45613E-02,-3.33751E-02,-5.73598E-02, 3.50360E-01,
     *  6.33053E-02, 2.16221E-02, 5.42577E-02,-5.74193E+00, 0.00000E+00,
     *  1.90891E-01,-1.39194E-02, 1.01102E+02, 8.16363E-02, 1.33717E-04,
     *  6.54403E-06, 3.10295E-03, 0.00000E+00, 0.00000E+00, 5.38205E-02,
     *  1.23910E-01,-1.39831E-02, 0.00000E+00, 0.00000E+00,-3.95915E-06,
     *  0.00000E+00,-7.14651E-01,-5.01027E-03, 0.00000E+00,-3.24756E-03,
     *  0.00000E+00, 0.00000E+00, 4.42173E-02,-1.31598E+01,-3.15626E-03,
     *  1.24574E-03,-1.47626E-03,-1.55461E-03, 6.40682E-02, 1.34898E-01,
     * -2.42415E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 6.13666E-04/
      DATA PB2/
     * -5.40373E-03, 2.61635E-05,-3.33012E-03, 0.00000E+00,-3.08101E-03,
     * -2.42679E-03,-3.36086E-03, 0.00000E+00,-1.18979E+03,-5.04738E-02,
     * -2.61547E-03,-1.03132E-03, 1.91583E-04,-8.38132E+01,-1.40517E-02,
     * -1.14167E-02,-4.08012E-03, 1.73522E-04,-1.39644E-02,-6.64128E-02,
     * -6.85152E-02,-1.34414E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  6.07916E+02,-4.12220E-03,-2.20996E-03, 0.00000E+00, 1.70277E+03,
     * -4.63015E-03, 0.00000E+00, 0.00000E+00,-2.25360E-03,-2.96204E-03,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  3.92786E-02, 1.31186E-02,-1.78086E-03, 0.00000E+00, 0.00000E+00,
     * -3.90083E-01,-2.84741E-02,-7.78400E+01,-1.02601E-03, 0.00000E+00/
      DATA PB3/
     * -7.26485E-04,-5.42181E-03,-5.59305E-03, 1.22825E-02, 1.23868E-02,
     *  6.68835E-03,-1.03303E-02,-9.51903E-03, 2.70021E-04,-2.57084E-02,
     * -1.32430E-02, 0.00000E+00,-3.81000E-02,-3.16810E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-9.05762E-04,-2.14590E-03,-1.17824E-03, 3.66732E+00,
     * -3.79729E-04,-6.13966E-03,-5.09082E-03,-1.96332E-03,-3.08280E-03,
     * -9.75222E-04, 4.03315E+00,-2.52710E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         N2 DENSITY
      DATA PC1/
     *  1.16112E+00, 0.00000E+00, 0.00000E+00, 3.33725E-02, 0.00000E+00,
     *  3.48637E-02,-5.44368E-03, 0.00000E+00,-6.73940E-02, 1.74754E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.74712E+02, 0.00000E+00,
     *  1.26733E-01, 0.00000E+00, 1.03154E+02, 5.52075E-02, 0.00000E+00,
     *  0.00000E+00, 8.13525E-04, 0.00000E+00, 0.00000E+00, 8.66784E-02,
     *  1.58727E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-2.50482E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.48894E-03,
     *  6.16053E-04,-5.79716E-04, 2.95482E-03, 8.47001E-02, 1.70147E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PC2/
     *  0.00000E+00, 2.47425E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PC3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         TLB
      DATA PD1/
     *  9.44846E-01, 0.00000E+00, 0.00000E+00,-3.08617E-02, 0.00000E+00,
     * -2.44019E-02, 6.48607E-03, 0.00000E+00, 3.08181E-02, 4.59392E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.74712E+02, 0.00000E+00,
     *  2.13260E-02, 0.00000E+00,-3.56958E+02, 0.00000E+00, 1.82278E-04,
     *  0.00000E+00, 3.07472E-04, 0.00000E+00, 0.00000E+00, 8.66784E-02,
     *  1.58727E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 3.83054E-03, 0.00000E+00, 0.00000E+00,
     * -1.93065E-03,-1.45090E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.23493E-03, 1.36736E-03, 8.47001E-02, 1.70147E-01,
     *  3.71469E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PD2/
     *  5.10250E-03, 2.47425E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 3.68756E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PD3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         O2 DENSITY
      DATA PE1/
     *  1.35580E+00, 1.44816E-01, 0.00000E+00, 6.07767E-02, 0.00000E+00,
     *  2.94777E-02, 7.46900E-02, 0.00000E+00,-9.23822E-02, 8.57342E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.38636E+01, 0.00000E+00,
     *  7.71653E-02, 0.00000E+00, 8.18751E+01, 1.87736E-02, 0.00000E+00,
     *  0.00000E+00, 1.49667E-02, 0.00000E+00, 0.00000E+00, 8.66784E-02,
     *  1.58727E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-3.67874E+02, 5.48158E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 8.47001E-02, 1.70147E-01,
     *  1.22631E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PE2/
     *  8.17187E-03, 3.71617E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.10826E-03,
     * -3.13640E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -7.35742E-02,-5.00266E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 1.94965E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PE3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         AR DENSITY
      DATA PF1/
     *  1.04761E+00, 2.00165E-01, 2.37697E-01, 3.68552E-02, 0.00000E+00,
     *  3.57202E-02,-2.14075E-01, 0.00000E+00,-1.08018E-01,-3.73981E-01,
     *  0.00000E+00, 3.10022E-02,-1.16305E-03,-2.07596E+01, 0.00000E+00,
     *  8.64502E-02, 0.00000E+00, 9.74908E+01, 5.16707E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 8.66784E-02,
     *  1.58727E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 3.46193E+02, 1.34297E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-3.48509E-03,
     * -1.54689E-04, 0.00000E+00, 0.00000E+00, 8.47001E-02, 1.70147E-01,
     *  1.47753E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PF2/
     *  1.89320E-02, 3.68181E-05, 1.32570E-02, 0.00000E+00, 0.00000E+00,
     *  3.59719E-03, 7.44328E-03,-1.00023E-03,-6.50528E+03, 0.00000E+00,
     *  1.03485E-02,-1.00983E-03,-4.06916E-03,-6.60864E+01,-1.71533E-02,
     *  1.10605E-02, 1.20300E-02,-5.20034E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -2.62769E+03, 7.13755E-03, 4.17999E-03, 0.00000E+00, 1.25910E+04,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-2.23595E-03, 4.60217E-03,
     *  5.71794E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -3.18353E-02,-2.35526E-02,-1.36189E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.03522E-02,-6.67837E+01,-1.09724E-03, 0.00000E+00/
      DATA PF3/
     * -1.38821E-02, 1.60468E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.51574E-02,
     * -5.44470E-04, 0.00000E+00, 7.28224E-02, 6.59413E-02, 0.00000E+00,
     * -5.15692E-03, 0.00000E+00, 0.00000E+00,-3.70367E+03, 0.00000E+00,
     *  0.00000E+00, 1.36131E-02, 5.38153E-03, 0.00000E+00, 4.76285E+00,
     * -1.75677E-02, 2.26301E-02, 0.00000E+00, 1.76631E-02, 4.77162E-03,
     *  0.00000E+00, 5.39354E+00, 0.00000E+00,-7.51710E-03, 0.00000E+00,
     *  0.00000E+00,-8.82736E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          H DENSITY
      DATA PG1/
     *  1.26376E+00,-2.14304E-01,-1.49984E-01, 2.30404E-01, 2.98237E-02,
     *  2.68673E-02, 2.96228E-01, 2.21900E-02,-2.07655E-02, 4.52506E-01,
     *  1.20105E-01, 3.24420E-02, 4.24816E-02,-9.14313E+00, 0.00000E+00,
     *  2.47178E-02,-2.88229E-02, 8.12805E+01, 5.10380E-02,-5.80611E-03,
     *  2.51236E-05,-1.24083E-02, 0.00000E+00, 0.00000E+00, 8.66784E-02,
     *  1.58727E-01,-3.48190E-02, 0.00000E+00, 0.00000E+00, 2.89885E-05,
     *  0.00000E+00, 1.53595E+02,-1.68604E-02, 0.00000E+00, 1.01015E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.84552E-04,
     * -1.22181E-03, 0.00000E+00, 0.00000E+00, 8.47001E-02, 1.70147E-01,
     * -1.04927E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,-5.91313E-03/
      DATA PG2/
     * -2.30501E-02, 3.14758E-05, 0.00000E+00, 0.00000E+00, 1.26956E-02,
     *  8.35489E-03, 3.10513E-04, 0.00000E+00, 3.42119E+03,-2.45017E-03,
     * -4.27154E-04, 5.45152E-04, 1.89896E-03, 2.89121E+01,-6.49973E-03,
     * -1.93855E-02,-1.48492E-02, 0.00000E+00,-5.10576E-02, 7.87306E-02,
     *  9.51981E-02,-1.49422E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  2.65503E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 6.37110E-03, 3.24789E-04,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  6.14274E-02, 1.00376E-02,-8.41083E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.27099E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PG3/
     * -3.94077E-03,-1.28601E-02,-7.97616E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-6.71465E-03,-1.69799E-03, 1.93772E-03, 3.81140E+00,
     * -7.79290E-03,-1.82589E-02,-1.25860E-02,-1.04311E-02,-3.02465E-03,
     *  2.43063E-03, 3.63237E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          N DENSITY
      DATA PH1/
     *  7.09557E+01,-3.26740E-01, 0.00000E+00,-5.16829E-01,-1.71664E-03,
     *  9.09310E-02,-6.71500E-01,-1.47771E-01,-9.27471E-02,-2.30862E-01,
     * -1.56410E-01, 1.34455E-02,-1.19717E-01, 2.52151E+00, 0.00000E+00,
     * -2.41582E-01, 5.92939E-02, 4.39756E+00, 9.15280E-02, 4.41292E-03,
     *  0.00000E+00, 8.66807E-03, 0.00000E+00, 0.00000E+00, 8.66784E-02,
     *  1.58727E-01, 9.74701E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 6.70217E+01,-1.31660E-03, 0.00000E+00,-1.65317E-02,
     *  0.00000E+00, 0.00000E+00, 8.50247E-02, 2.77428E+01, 4.98658E-03,
     *  6.15115E-03, 9.50156E-03,-2.12723E-02, 8.47001E-02, 1.70147E-01,
     * -2.38645E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.37380E-03/
      DATA PH2/
     * -8.41918E-03, 2.80145E-05, 7.12383E-03, 0.00000E+00,-1.66209E-02,
     *  1.03533E-04,-1.68898E-02, 0.00000E+00, 3.64526E+03, 0.00000E+00,
     *  6.54077E-03, 3.69130E-04, 9.94419E-04, 8.42803E+01,-1.16124E-02,
     * -7.74414E-03,-1.68844E-03, 1.42809E-03,-1.92955E-03, 1.17225E-01,
     * -2.41512E-02, 1.50521E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.60261E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-3.54403E-04,-1.87270E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  2.76439E-02, 6.43207E-03,-3.54300E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-2.80221E-02, 8.11228E+01,-6.75255E-04, 0.00000E+00/
      DATA PH3/
     * -1.05162E-02,-3.48292E-03,-6.97321E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.45546E-03,-1.31970E-02,-3.57751E-03,-1.09021E+00,
     * -1.50181E-02,-7.12841E-03,-6.64590E-03,-3.52610E-03,-1.87773E-02,
     * -2.22432E-03,-3.93895E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C        HOT O DENSITY
      DATA PI1/
     *  6.04050E-02, 1.57034E+00, 2.99387E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.51018E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-8.61650E+00, 1.26454E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 5.50878E-03, 0.00000E+00, 0.00000E+00, 8.66784E-02,
     *  1.58727E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 6.23881E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 8.47001E-02, 1.70147E-01,
     * -9.45934E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PI2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PI3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          S PARAM  
      DATA PJ1/
     *  9.56827E-01, 6.20637E-02, 3.18433E-02, 0.00000E+00, 0.00000E+00,
     *  3.94900E-02, 0.00000E+00, 0.00000E+00,-9.24882E-03,-7.94023E-03,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.74712E+02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.74677E-03, 0.00000E+00, 1.54951E-02, 8.66784E-02,
     *  1.58727E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-6.99007E-04, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 1.24362E-02,-5.28756E-03, 8.47001E-02, 1.70147E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PJ2/
     *  0.00000E+00, 2.47425E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PJ3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TURBO
      DATA PK1/
     *  1.09930E+00, 3.90631E+00, 3.07165E+00, 9.86161E-01, 1.63536E+01,
     *  4.63830E+00, 1.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 1.28840E+00, 3.10302E-02, 1.18339E-01,
     *  1.00000E+00, 7.00000E-01, 1.15020E+00, 3.44689E+00, 1.28840E+00,
     *  1.00000E+00, 1.08738E+00, 1.22947E+00, 1.10016E+00, 7.34129E-01,
     *  1.15241E+00, 2.22784E+00, 7.95046E-01, 4.01612E+00, 4.47749E+00,
     *  1.23435E+02,-7.60535E-02, 1.68986E-06, 7.44294E-01, 1.03604E+00,
     *  1.72783E+02, 1.15020E+00, 3.44689E+00,-7.46230E-01, 9.49154E-01/
C         LOWER BOUNDARY
      DATA PTM/
     L  1.04130E+03, 3.86000E+02, 1.95000E+02, 1.66728E+01, 2.13000E+02,
     L  1.20000E+02, 2.40000E+02, 1.87000E+02,-2.00000E+00, 0.00000E+00/
      DATA PDM/
     L  2.45600E+07, 6.71072E-06, 1.00000E+02, 0.00000E+00, 1.10000E+02,
     L  1.00000E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  8.59400E+10, 1.00000E+00, 1.05000E+02,-8.00000E+00, 1.10000E+02,
     L  1.00000E+01, 9.00000E+01, 2.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  2.81000E+11, 0.00000E+00, 1.05000E+02, 2.80000E+01, 2.89500E+01,
     L  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  3.30000E+10, 2.68270E-01, 1.05000E+02, 1.00000E+00, 1.10000E+02,
     L  1.00000E+01, 1.10000E+02,-1.00000E+01, 0.00000E+00, 0.00000E+00,
C
     L  1.33000E+09, 1.19615E-02, 1.05000E+02, 0.00000E+00, 1.10000E+02,
     L  1.00000E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  1.76100E+05, 1.00000E+00, 9.50000E+01,-8.00000E+00, 1.10000E+02,
     L  1.00000E+01, 9.00000E+01, 2.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  1.00000E+07, 1.00000E+00, 1.05000E+02,-8.00000E+00, 1.10000E+02,
     L  1.00000E+01, 9.00000E+01, 2.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  1.00000E+06, 1.00000E+00, 1.05000E+02,-8.00000E+00, 5.50000E+02,
     L  7.60000E+01, 9.00000E+01, 2.00000E+00, 0.00000E+00, 4.00000E+03/
C         TN1(2)
      DATA PL1/
     *  1.00858E+00, 4.56011E-02,-2.22972E-02,-5.44388E-02, 5.23136E-04,
     * -1.88849E-02, 5.23707E-02,-9.43646E-03, 6.31707E-03,-7.80460E-02,
     * -4.88430E-02, 0.00000E+00, 0.00000E+00,-7.60250E+00, 0.00000E+00,
     * -1.44635E-02,-1.76843E-02,-1.21517E+02, 2.85647E-02, 0.00000E+00,
     *  0.00000E+00, 6.31792E-04, 0.00000E+00, 5.77197E-03, 8.66784E-02,
     *  1.58727E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-8.90272E+03, 3.30611E-03, 3.02172E-03, 0.00000E+00,
     * -2.13673E-03,-3.20910E-04, 0.00000E+00, 0.00000E+00, 2.76034E-03,
     *  2.82487E-03,-2.97592E-04,-4.21534E-03, 8.47001E-02, 1.70147E-01,
     *  8.96456E-03, 0.00000E+00,-1.08596E-02, 0.00000E+00, 0.00000E+00/
      DATA PL2/
     *  5.57917E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 9.65405E-03, 0.00000E+00, 0.00000E+00, 2.00000E+00/
C         TN1(3)
      DATA PM1/
     *  9.39664E-01, 8.56514E-02,-6.79989E-03, 2.65929E-02,-4.74283E-03,
     *  1.21855E-02,-2.14905E-02, 6.49651E-03,-2.05477E-02,-4.24952E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.19148E+01, 0.00000E+00,
     *  1.18777E-02,-7.28230E-02,-8.15965E+01, 1.73887E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.44691E-02, 2.80259E-04, 8.66784E-02,
     *  1.58727E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.16584E+02, 3.18713E-03, 7.37479E-03, 0.00000E+00,
     * -2.55018E-03,-3.92806E-03, 0.00000E+00, 0.00000E+00,-2.89757E-03,
     * -1.33549E-03, 1.02661E-03, 3.53775E-04, 8.47001E-02, 1.70147E-01,
     * -9.17497E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PM2/
     *  3.56082E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.00902E-02, 0.00000E+00, 0.00000E+00, 2.00000E+00/
C         TN1(4)
      DATA PN1/
     *  9.85982E-01,-4.55435E-02, 1.21106E-02, 2.04127E-02,-2.40836E-03,
     *  1.11383E-02,-4.51926E-02, 1.35074E-02,-6.54139E-03, 1.15275E-01,
     *  1.28247E-01, 0.00000E+00, 0.00000E+00,-5.30705E+00, 0.00000E+00,
     * -3.79332E-02,-6.24741E-02, 7.71062E-01, 2.96315E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 6.81051E-03,-4.34767E-03, 8.66784E-02,
     *  1.58727E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 1.07003E+01,-2.76907E-03, 4.32474E-04, 0.00000E+00,
     *  1.31497E-03,-6.47517E-04, 0.00000E+00,-2.20621E+01,-1.10804E-03,
     * -8.09338E-04, 4.18184E-04, 4.29650E-03, 8.47001E-02, 1.70147E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PN2/
     * -4.04337E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-9.52550E-04,
     *  8.56253E-04, 4.33114E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.21223E-03,
     *  2.38694E-04, 9.15245E-04, 1.28385E-03, 8.67668E-04,-5.61425E-06,
     *  1.04445E+00, 3.41112E+01, 0.00000E+00,-8.40704E-01,-2.39639E+02,
     *  7.06668E-01,-2.05873E+01,-3.63696E-01, 2.39245E+01, 0.00000E+00,
     * -1.06657E-03,-7.67292E-04, 1.54534E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.00000E+00/
C         TN1(5) TN2(1)
      DATA PO1/
     *  1.00320E+00, 3.83501E-02,-2.38983E-03, 2.83950E-03, 4.20956E-03,
     *  5.86619E-04, 2.19054E-02,-1.00946E-02,-3.50259E-03, 4.17392E-02,
     * -8.44404E-03, 0.00000E+00, 0.00000E+00, 4.96949E+00, 0.00000E+00,
     * -7.06478E-03,-1.46494E-02, 3.13258E+01,-1.86493E-03, 0.00000E+00,
     * -1.67499E-02, 0.00000E+00, 0.00000E+00, 5.12686E-04, 8.66784E-02,
     *  1.58727E-01,-4.64167E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  4.37353E-03,-1.99069E+02, 0.00000E+00,-5.34884E-03, 0.00000E+00,
     *  1.62458E-03, 2.93016E-03, 2.67926E-03, 5.90449E+02, 0.00000E+00,
     *  0.00000E+00,-1.17266E-03,-3.58890E-04, 8.47001E-02, 1.70147E-01,
     *  0.00000E+00, 0.00000E+00, 1.38673E-02, 0.00000E+00, 0.00000E+00/
      DATA PO2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.60571E-03,
     *  6.28078E-04, 5.05469E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.57829E-03,
     * -4.00855E-04, 5.04077E-05,-1.39001E-03,-2.33406E-03,-4.81197E-04,
     *  1.46758E+00, 6.20332E+00, 0.00000E+00, 3.66476E-01,-6.19760E+01,
     *  3.09198E-01,-1.98999E+01, 0.00000E+00,-3.29933E+02, 0.00000E+00,
     * -1.10080E-03,-9.39310E-05, 1.39638E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.00000E+00/
C          TN2(2)
      DATA PP1/
     *  9.81637E-01,-1.41317E-03, 3.87323E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-3.58707E-02,
     * -8.63658E-03, 0.00000E+00, 0.00000E+00,-2.02226E+00, 0.00000E+00,
     * -8.69424E-03,-1.91397E-02, 8.76779E+01, 4.52188E-03, 0.00000E+00,
     *  2.23760E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-7.07572E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -4.11210E-03, 3.50060E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-8.36657E-03, 1.61347E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.45130E-02, 0.00000E+00, 0.00000E+00/
      DATA PP2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.24152E-03,
     *  6.43365E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.33255E-03,
     *  2.42657E-03, 1.60666E-03,-1.85728E-03,-1.46874E-03,-4.79163E-06,
     *  1.22464E+00, 3.53510E+01, 0.00000E+00, 4.49223E-01,-4.77466E+01,
     *  4.70681E-01, 8.41861E+00,-2.88198E-01, 1.67854E+02, 0.00000E+00,
     *  7.11493E-04, 6.05601E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.00000E+00/
C          TN2(3)
      DATA PQ1/
     *  1.00422E+00,-7.11212E-03, 5.24480E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-5.28914E-02,
     * -2.41301E-02, 0.00000E+00, 0.00000E+00,-2.12219E+01,-1.03830E-02,
     * -3.28077E-03, 1.65727E-02, 1.68564E+00,-6.68154E-03, 0.00000E+00,
     *  1.45155E-02, 0.00000E+00, 8.42365E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-4.34645E-03, 0.00000E+00, 0.00000E+00, 2.16780E-02,
     *  0.00000E+00,-1.38459E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 7.04573E-03,-4.73204E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 1.08767E-02, 0.00000E+00, 0.00000E+00/
      DATA PQ2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-8.08279E-03,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.21769E-04,
     * -2.27387E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.26769E-03,
     *  3.16901E-03, 4.60316E-04,-1.01431E-04, 1.02131E-03, 9.96601E-04,
     *  1.25707E+00, 2.50114E+01, 0.00000E+00, 4.24472E-01,-2.77655E+01,
     *  3.44625E-01, 2.75412E+01, 0.00000E+00, 7.94251E+02, 0.00000E+00,
     *  2.45835E-03, 1.38871E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.00000E+00/
C          TN2(4) TN3(1)
      DATA PR1/
     *  1.01890E+00,-2.46603E-02, 1.00078E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-6.70977E-02,
     * -4.02286E-02, 0.00000E+00, 0.00000E+00,-2.29466E+01,-7.47019E-03,
     *  2.26580E-03, 2.63931E-02, 3.72625E+01,-6.39041E-03, 0.00000E+00,
     *  9.58383E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.85291E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 1.39717E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 9.19771E-03,-3.69121E+02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.57067E-02, 0.00000E+00, 0.00000E+00/
      DATA PR2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-7.07265E-03,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.92953E-03,
     * -2.77739E-03,-4.40092E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.47280E-03,
     *  2.95035E-04,-1.81246E-03, 2.81945E-03, 4.27296E-03, 9.78863E-04,
     *  1.40545E+00,-6.19173E+00, 0.00000E+00, 0.00000E+00,-7.93632E+01,
     *  4.44643E-01,-4.03085E+02, 0.00000E+00, 1.15603E+01, 0.00000E+00,
     *  2.25068E-03, 8.48557E-04,-2.98493E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.00000E+00/
C          TN3(2)
      DATA PS1/
     *  9.75801E-01, 3.80680E-02,-3.05198E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.85575E-02,
     *  5.04057E-02, 0.00000E+00, 0.00000E+00,-1.76046E+02, 1.44594E-02,
     * -1.48297E-03,-3.68560E-03, 3.02185E+01,-3.23338E-03, 0.00000E+00,
     *  1.53569E-02, 0.00000E+00,-1.15558E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 4.89620E-03, 0.00000E+00, 0.00000E+00,-1.00616E-02,
     * -8.21324E-03,-1.57757E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 6.63564E-03, 4.58410E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-2.51280E-02, 0.00000E+00, 0.00000E+00/
      DATA PS2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 9.91215E-03,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-8.73148E-04,
     * -1.29648E-03,-7.32026E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-4.68110E-03,
     * -4.66003E-03,-1.31567E-03,-7.39390E-04, 6.32499E-04,-4.65588E-04,
     * -1.29785E+00,-1.57139E+02, 0.00000E+00, 2.58350E-01,-3.69453E+01,
     *  4.10672E-01, 9.78196E+00,-1.52064E-01,-3.85084E+03, 0.00000E+00,
     * -8.52706E-04,-1.40945E-03,-7.26786E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.00000E+00/
C          TN3(3)
      DATA PU1/
     *  9.60722E-01, 7.03757E-02,-3.00266E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.22671E-02,
     *  4.10423E-02, 0.00000E+00, 0.00000E+00,-1.63070E+02, 1.06073E-02,
     *  5.40747E-04, 7.79481E-03, 1.44908E+02, 1.51484E-04, 0.00000E+00,
     *  1.97547E-02, 0.00000E+00,-1.41844E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 5.77884E-03, 0.00000E+00, 0.00000E+00, 9.74319E-03,
     *  0.00000E+00,-2.88015E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-4.44902E-03,-2.92760E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 2.34419E-02, 0.00000E+00, 0.00000E+00/
      DATA PU2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.36685E-03,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-4.65325E-04,
     * -5.50628E-04, 3.31465E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.06179E-03,
     * -3.08575E-03,-7.93589E-04,-1.08629E-04, 5.95511E-04,-9.05050E-04,
     *  1.18997E+00, 4.15924E+01, 0.00000E+00,-4.72064E-01,-9.47150E+02,
     *  3.98723E-01, 1.98304E+01, 0.00000E+00, 3.73219E+03, 0.00000E+00,
     * -1.50040E-03,-1.14933E-03,-1.56769E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.00000E+00/
C          TN3(4)
      DATA PV1/
     *  1.03123E+00,-7.05124E-02, 8.71615E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-3.82621E-02,
     * -9.80975E-03, 0.00000E+00, 0.00000E+00, 2.89286E+01, 9.57341E-03,
     *  0.00000E+00, 0.00000E+00, 8.66153E+01, 7.91938E-04, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 4.68917E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 7.86638E-03, 0.00000E+00, 0.00000E+00, 9.90827E-03,
     *  0.00000E+00, 6.55573E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-4.00200E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 7.07457E-03, 0.00000E+00, 0.00000E+00/
      DATA PV2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.72268E-03,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.04970E-04,
     *  1.21560E-03,-8.05579E-06, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.49941E-03,
     * -4.57256E-04,-1.59311E-04, 2.96481E-04,-1.77318E-03,-6.37918E-04,
     *  1.02395E+00, 1.28172E+01, 0.00000E+00, 1.49903E-01,-2.63818E+01,
     *  0.00000E+00, 4.70628E+01,-2.22139E-01, 4.82292E-02, 0.00000E+00,
     * -8.67075E-04,-5.86479E-04, 5.32462E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.00000E+00/
C          TN3(5) SURFACE TEMP TSL
      DATA PW1/
     *  1.00828E+00,-9.10404E-02,-2.26549E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.32420E-02,
     * -9.08925E-03, 0.00000E+00, 0.00000E+00, 3.36105E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.24957E+01,-5.87939E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.79765E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.01237E+03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.75553E-02, 0.00000E+00, 0.00000E+00/
      DATA PW2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.29699E-03,
     *  1.26659E-03, 2.68402E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.17894E-03,
     *  1.48746E-03, 1.06478E-04, 1.34743E-04,-2.20939E-03,-6.23523E-04,
     *  6.36539E-01, 1.13621E+01, 0.00000E+00,-3.93777E-01, 2.38687E+03,
     *  0.00000E+00, 6.61865E+02,-1.21434E-01, 9.27608E+00, 0.00000E+00,
     *  1.68478E-04, 1.24892E-03, 1.71345E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.00000E+00/
C          TGN3(2) SURFACE GRAD TSLG
      DATA PX1/
     *  1.57293E+00,-6.78400E-01, 6.47500E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-7.62974E-02,
     * -3.60423E-01, 0.00000E+00, 0.00000E+00, 1.28358E+02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 4.68038E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.67898E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.90994E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 3.15706E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PX2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.00000E+00/
C          TGN2(1) TGN1(2)
      DATA PY1/
     *  8.60028E-01, 3.77052E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.17570E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 7.77757E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 1.01024E+02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 6.54251E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PY2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.56959E-02,
     *  1.91001E-02, 3.15971E-02, 1.00982E-02,-6.71565E-03, 2.57693E-03,
     *  1.38692E+00, 2.82132E-01, 0.00000E+00, 0.00000E+00, 3.81511E+02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.00000E+00/
C          TGN3(1) TGN2(2)
      DATA PZ1/
     *  1.06029E+00,-5.25231E-02, 3.73034E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.31072E-02,
     * -3.88409E-01, 0.00000E+00, 0.00000E+00,-1.65295E+02,-2.13801E-01,
     * -4.38916E-02,-3.22716E-01,-8.82393E+01, 1.18458E-01, 0.00000E+00,
     * -4.35863E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.19782E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.62229E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-5.37443E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-4.55788E-01, 0.00000E+00, 0.00000E+00/
      DATA PZ2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.84009E-02,
     *  3.96733E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.05494E-02,
     *  7.39617E-02, 1.92200E-02,-8.46151E-03,-1.34244E-02, 1.96338E-02,
     *  1.50421E+00, 1.88368E+01, 0.00000E+00, 0.00000E+00,-5.13114E+01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  5.11923E-02, 3.61225E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.00000E+00/
C          SEMIANNUAL MULT SAM
      DATA PAA1/
     *  1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00,
     *  1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00,
     *  1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00,
     *  1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00,
     *  1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00,
     *  1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00,
     *  1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00,
     *  1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00,
     *  1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00,
     *  1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00/
      DATA PAA2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         MIDDLE ATMOSPHERE AVERAGES
      DATA PAVGM/
     M  2.61000E+02, 2.64000E+02, 2.29000E+02, 2.17000E+02, 2.17000E+02,
     M  2.23000E+02, 2.86760E+02,-2.93940E+00, 2.50000E+00, 0.00000E+00/
      END

