C      TEST DRIVER FOR GWS5 (HWM93 HORIZONTAL WIND MODEL)
      DIMENSION W(2,20)
      DIMENSION IDAY(20),UT(20),ALT(20),XLAT(20),XLONG(20),XLST(20),
     $ F107A(20),F107(20),AP(20)    
      COMMON/HWMC/WBT(2),WCT(2)
      DATA IDAY/172,81,8*172,3*81.,7*182./
      DATA UT/29000.,29000.,75000.,17*29000./
      DATA ALT/400.,400.,400.,200.,6*400.,5*100.,80.,60.,40.,20.,0/
      DATA XLAT/4*60.,0.,5*60.,4*45.,0,45.,45.,-45.,45.,45./
      DATA XLONG/5*-70.,0.,4*-70.,5*0,90.,90.,0,-90.,-90./
      DATA XLST/6*16.,4.,3*16.,0.,6.,9.,12.,0,5*0/
      DATA F107A/7*150.,70.,150.,150.,10*150./
      DATA F107/8*150.,180.,150.,10*150./
      DATA AP/9*4.,40.,4.,40.,8*4./
      
      DIMENSION OULOULOUW(2, 10), OULOULOUALT(10)
      DATA OULOULOUALT/0.,10.,25.,50.,75.,100.,200.,300.,400.,500./
      
      WRITE(6, *) 'ALT MW ZW'
      DO I=1,10
        CALL GWS5(0, 43200., OULOULOUALT(I), 0., 0., 12., 150., 150.,
     &            4., OULOULOUW(1, I))
        WRITE(6, *) OULOULOUALT(I), (OULOULOUW(J, I), J=1,2)
        WRITE(6, *) ''
      ENDDO
      
      DO 10 I=1,20
         CALL GWS5(IDAY(I),UT(I),ALT(I),XLAT(I),XLONG(I),XLST(I),
     $             F107A(I),F107(I),AP(I),W(1,I))
         WRITE(6,100) W(1,I),WBT(1),WCT(1),W(2,I),WBT(2),WCT(2)
   10 CONTINUE
      WRITE(6,200) (IDAY(I),I=1,5)
      WRITE(6,201) (UT(I),I=1,5)
      WRITE(6,202) (ALT(I),I=1,5)
      WRITE(6,203) (XLAT(I),I=1,5)
      WRITE(6,204) (XLONG(I),I=1,5)
      WRITE(6,205) (XLST(I),I=1,5)
      WRITE(6,206) (F107A(I),I=1,5)
      WRITE(6,207) (F107(I),I=1,5)
      WRITE(6,210) (AP(I),I=1,5)
      WRITE(6,208) (W(1,I),I=1,5)
      WRITE(6,209) (W(2,I),I=1,5)
      WRITE(6,200) (IDAY(I),I=6,10)
      WRITE(6,201) (UT(I),I=6,10)
      WRITE(6,202) (ALT(I),I=6,10)
      WRITE(6,203) (XLAT(I),I=6,10)
      WRITE(6,204) (XLONG(I),I=6,10)
      WRITE(6,205) (XLST(I),I=6,10)
      WRITE(6,206) (F107A(I),I=6,10)
      WRITE(6,207) (F107(I),I=6,10)
      WRITE(6,210) (AP(I),I=6,10)
      WRITE(6,208) (W(1,I),I=6,10)
      WRITE(6,209) (W(2,I),I=6,10)
      WRITE(6,200) (IDAY(I),I=11,15)
      WRITE(6,201) (UT(I),I=11,15)
      WRITE(6,202) (ALT(I),I=11,15)
      WRITE(6,203) (XLAT(I),I=11,15)
      WRITE(6,204) (XLONG(I),I=11,15)
      WRITE(6,205) (XLST(I),I=11,15)
      WRITE(6,206) (F107A(I),I=11,15)
      WRITE(6,207) (F107(I),I=11,15)
      WRITE(6,210) (AP(I),I=11,15)
      WRITE(6,208) (W(1,I),I=11,15)
      WRITE(6,209) (W(2,I),I=11,15)
      WRITE(6,200) (IDAY(I),I=16,20)
      WRITE(6,201) (UT(I),I=16,20)
      WRITE(6,202) (ALT(I),I=16,20)
      WRITE(6,203) (XLAT(I),I=16,20)
      WRITE(6,204) (XLONG(I),I=16,20)
      WRITE(6,205) (XLST(I),I=16,20)
      WRITE(6,206) (F107A(I),I=16,20)
      WRITE(6,207) (F107(I),I=16,20)
      WRITE(6,210) (AP(I),I=16,20)
      WRITE(6,208) (W(1,I),I=16,20)
      WRITE(6,209) (W(2,I),I=16,20)
  100 FORMAT(1X,6F10.2)
  200 FORMAT(//' DAY  ',5I12)
  201 FORMAT(' UT   ',5F12.0)
  202 FORMAT(' ALT  ',5F12.0)
  203 FORMAT(' LAT  ',5F12.0)
  204 FORMAT(' LONG ',5F12.0)
  205 FORMAT(' LST  ',5F12.0)
  206 FORMAT(' F107A',5F12.0)
  207 FORMAT(' F107 ',5F12.0)
  210 FORMAT(' AP   ',5F12.0)
  208 FORMAT(/' MERID',5F12.2)
  209 FORMAT(' ZONAL',5F12.2)
      STOP
      END
