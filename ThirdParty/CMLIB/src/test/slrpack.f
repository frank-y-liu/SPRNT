C
C   DRIVER FOR TESTING CMLIB ROUTINES
C     DLINFS  LINFS   RYBAR
C     RGM     RWILL   RYARNG
C     RYORK   WWILL   WYORK
C
C    ONE INPUT DATA CARD IS REQUIRED
C         READ(LIN,1) KPRINT,TIMES
C    1    FORMAT(I1,E10.0)
C
C     KPRINT = 0   NO PRINTING
C              1   NO PRINTING FOR PASSED TESTS, SHORT MESSAGE
C                  FOR FAILED TESTS
C              2   PRINT SHORT MESSAGE FOR PASSED TESTS, FULLER
C                  INFORMATION FOR FAILED TESTS
C              3   PRINT COMPLETE QUICK-CHECK RESULTS
C
C                ***IMPORTANT NOTE***
C         ALL QUICK CHECKS USE ROUTINES R2MACH AND D2MACH
C         TO SET THE ERROR TOLERANCES.
C     TIMES IS A CONSTANT MULTIPLIER THAT CAN BE USED TO SCALE THE
C     VALUES OF R1MACH AND D1MACH SO THAT
C               R2MACH(I) = R1MACH(I) * TIMES   FOR I=3,4,5
C               D2MACH(I) = D1MACH(I) * TIMES   FOR I=3,4,5
C     THIS MAKES IT EASILY POSSIBLE TO CHANGE THE ERROR TOLERANCES
C     USED IN THE QUICK CHECKS.
C     IF TIMES .LE. 0.0 THEN TIMES IS DEFAULTED TO 1.0
C
C              ***END NOTE***
C
      COMMON/UNIT/LUN
      COMMON/MSG/ICNT,JTEST(38)
      COMMON/XXMULT/TIMES
      LUN=I1MACH(2)
      LIN=I1MACH(1)
      ITEST=1
C
C     READ KPRINT,TIMES PARAMETERS FROM DATA CARD..
C
      READ(LIN,1) KPRINT,TIMES
1     FORMAT(I1,E10.0)
      IF(TIMES.LE.0.) TIMES=1.
      CALL XSETUN(LUN)
      CALL XSETF(1)
      CALL XERMAX(1000)
C   TEST RGM
      CALL SLRQX1(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C   TEST RYORK
      CALL SLRQX2(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C   TEST RWILL
      CALL SLRQX3(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C   TEST LINFS
      CALL SLRQX4(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C   TEST DLINFS
      CALL SLRQX5(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C
      IF(KPRINT.GE.1.AND.ITEST.NE.1) WRITE(LUN,2)
2     FORMAT(/' ***** WARNING -- AT LEAST ONE TEST FOR SUBLIBRARY SLRPACK
     1 HAS FAILED ***** ')
      IF(KPRINT.GE.1.AND.ITEST.EQ.1) WRITE(LUN,3)
3     FORMAT(/' ----- SUBLIBRARY SLRPACK PASSED ALL TESTS ----- ')
      END
      DOUBLE PRECISION FUNCTION D2MACH(I)
      DOUBLE PRECISION D1MACH
      COMMON/XXMULT/TIMES
      D2MACH=D1MACH(I)
      IF(I.EQ.1.OR. I.EQ.2) RETURN
      D2MACH = D2MACH * DBLE(TIMES)
      RETURN
      END
      REAL FUNCTION R2MACH(I)
      COMMON/XXMULT/TIMES
      R2MACH=R1MACH(I)
      IF(I.EQ.1.OR. I.EQ.2) RETURN
      R2MACH = R2MACH * TIMES
      RETURN
      END
      SUBROUTINE SLRQX1(KPRINT,IPASS)
      COMMON/UNIT/LUN
      INTEGER KPRINT,IPASS
      DIMENSION X(10), Y(10), OUTPUT(7)
      DATA X/0.0,0.9,1.8,2.6,3.3,4.4,5.2,6.1,6.5,7.4/
      DATA Y/5.9,5.4,4.4,4.6,3.5,3.7,2.8,2.8,2.4,1.5/
      IPASS = 1
          N = 10
        EPS = R2MACH(4)
         BT =  -.552576514441537
         AT =  5.810842285166643
      XBART =  3.819999999999979
      YBART =  3.699999999999983
       SDBT =   .037679130341852
       SDAT =   .881503409843887
         RT =  -.976475222674530
C-----
C----- COMPUTE GEOMETRIC MEAN SLOPE, INTERCEPT, STANDARD DEVIATIONS
C-----
      CALL RGM (N,X,Y,OUTPUT,IER)
      ERR = SQRT(EPS)
      IF (IER .EQ. 0                               .AND.
     1    ABS(OUTPUT(2) - AT) / AT        .LE. ERR .AND.
     2    ABS(OUTPUT(1) - BT) / BT        .LE. ERR .AND.
     3    ABS(OUTPUT(4) - SDAT) / SDAT    .LE. ERR .AND.
     4    ABS(OUTPUT(3) - SDBT) / SDBT    .LE. ERR .AND.
     5    ABS(OUTPUT(5) - XBART) / XBART  .LE. ERR .AND.
     6    ABS(OUTPUT(6) - YBART) / YBART  .LE. ERR .AND.
     7    ABS(OUTPUT(7) - RT) / RT        .LE. ERR) GO TO 10
C-----
C----- EXECUTION ERROR
C-----
      IF(KPRINT.GE.1) WRITE(LUN,900)
  900 FORMAT(/,' ** TEST 1 FOR SUBLIBRARY SLRPACK FAILED ** ')
      IF(KPRINT.GE.2) WRITE (LUN,901)
  901 FORMAT (' EXECUTION ERROR - ',/
     1        ' SUBROUTINE RGM RESULTS DO NOT AGREE WITH KNOWN RESULTS')
      IF(KPRINT.GE.2) WRITE(LUN,903) OUTPUT(2), AT, OUTPUT(1), BT,
     1 OUTPUT(4), SDAT, OUTPUT(3), SDBT, OUTPUT(5), XBART, OUTPUT(6),
     2 YBART, OUTPUT(7), RT
 903  FORMAT(8X, 'RGM RESULTS', 24X, 'EXPECTED RESULTS', /,
     1 7(3X, F20.15, 17X, F20.15, /))
      IF(KPRINT.GE.2) WRITE(LUN, 904) IER
      IPASS=2
 904  FORMAT(10X, 'RGM ERROR CONDITION ', 5X, I3)
      GO TO 999
C-----
C----- PRINT CORRECT EXECUTION MESSAGE
C-----
   10 IF(KPRINT.GE.2) WRITE (LUN,902)
  902 FORMAT (' ** TEST 1 FOR SUBLIBRARY SLRPACK PASSED ** ')
      GO TO 999
  999 RETURN
      END
      SUBROUTINE SLRQX2(KPRINT,IPASS)
      COMMON/UNIT/LUN
      INTEGER KPRINT,IPASS
      DIMENSION X(10), Y(10), WX(10), WY(10), RES(10,2), WORK(109),
     1          OUTPUT(7)
      DATA X/0.0,0.9,1.8,2.6,3.3,4.4,5.2,6.1,6.5,7.4/
      DATA Y/5.9,5.4,4.4,4.6,3.5,3.7,2.8,2.8,2.4,1.5/
      DATA WX/1000.,1000.,500.,800.,200.,80.,60.,20.,1.8,1.0/
      DATA WY/1.0,1.8,4.0,8.0,20.,20.,70.,70.,100.,500./
       IPASS = 1
           N = 10
        IRES = 1
        MAXD = 10
      MXITER = 10
      BSTART = -0.546
      DELMAX = 0.0001
         EPS = R2MACH(4)
          BT = -.480534835891692
        SDBT =  .071006409434287
          AT = 5.4799176942544
        SDAT = 0.361872433019652
      XBAR5T = 4.910991714976291
      YBAR5T = 3.120015096432942
C-----
C----- COMPUTE SLOPE, Y-INTERCEPT, VARIANCES, AND RESIDUALS (IF REQUESTED)
C----- BY YORK PROCEDURE WITH INITIAL SLOPE VALUE BSTART
C-----
  300 CALL RYORK (N, X, Y, WX, WY, BSTART, MXITER, DELMAX, IRES,
     1            MAXD, WORK, OUTPUT, RES, IER)
C-----
C----- CHECK CORRECTNESS OF PROGRAM EXECUTION
C-----
       ERR = SQRT(EPS)
       IF (IER .EQ. 0                                  .AND.
     1     INT(OUTPUT(7)) .EQ. 5                       .AND.
     2     ABS(OUTPUT(2) - AT) / AT         .LE. ERR   .AND.
     3     ABS(OUTPUT(1) - BT) / BT         .LE. ERR   .AND.
     4     ABS(OUTPUT(4) - SDAT) / SDAT     .LE. ERR   .AND.
     5     ABS(OUTPUT(3) - SDBT) / SDBT     .LE. ERR   .AND.
     6     ABS(OUTPUT(5) - XBAR5T) / XBAR5T .LE. ERR   .AND.
     7     ABS(OUTPUT(6) - YBAR5T) / YBAR5T .LE. ERR)  GO TO 10
C-----
C----- EXECUTION ERROR
C-----
      IF(KPRINT.GE.1) WRITE(LUN,900)
 900  FORMAT(' ** TEST 2 FOR SUBLIBRARY SLRPACK FAILED ** ')
      IF(KPRINT.GE.2) WRITE (LUN,901)
  901 FORMAT (' EXECUTION ERROR - SUBROUTINE RYORK  ')
      IF(KPRINT.GE.2) WRITE (LUN,902)
  902 FORMAT (' RESULTS DO NOT AGREE WITH KNOWN RESULTS')
      IF(KPRINT.GT.2) WRITE(LUN,904) OUTPUT(2), AT, OUTPUT(1), BT,
     1OUTPUT(4),SDAT,OUTPUT(3), SDBT, OUTPUT(5), XBAR5T, OUTPUT(6), 
     2YBAR5T
 904  FORMAT(8X, 'RYORK RESULTS', 22X, 'EXPECTED RESULTS', /,
     1  6(3X, F20.15, 17X, F20.15, /))
      IF(KPRINT.GT.2) WRITE(LUN,905) INT(OUTPUT(7)), IER
 905  FORMAT(3X, 'NUMBER OF ITERATIONS', 5X, I3, 5X,
     1      'RYORK ERROR CONDITION', 5X, I3)
      IPASS=2
      GO TO 999
C-----
C----- PRINT CORRECT EXECUTION MESSAGE
C-----
   10 IF(KPRINT.GE.2) WRITE (LUN,903)
  903 FORMAT (' ** TEST 2 FOR SUBLIBRARY SLRPACK PASSED ** ')
  999 RETURN
      END
      SUBROUTINE SLRQX3(KPRINT,IPASS)
      COMMON/UNIT/LUN
      INTEGER KPRINT,IPASS
      DIMENSION OUTPUT(11), WORK(50), X(10), Y(10), U(10), V(10),
     1          XRES(10), YRES(10)
      DATA X/0.0,0.9,1.8,2.6,3.3,4.4,5.2,6.1,6.5,7.4/
      DATA Y/5.9,5.4,4.4,4.6,3.5,3.7,2.8,2.8,2.4,1.5/
      DATA U/.001,.001,.002,.00125,.005,.0125,.01667,.05,.55556,1.0/
      DATA V/1.0,.55556,.250,.125,.05,.05,.01429,.01429,.01,.002/
       IPASS = 1
           N = 10
      MXITER = 10
      BSTART = -0.546
      DELMAX = 0.0001
         EPS = R2MACH(4)
          BT = -.480532799091556
      AJSDBT =  .057619766382573
        SDBT =  .070172013548005
          AT = 5.479914011835746
      AJSDAT =  .291943155534277
        SDAT =  .355541862654874
      XBAR4T = 4.910865781673408
      YBAR4T = 3.120081931805274
       RMSXT =  .288752957347914
       RMSYT =  .277613150086520
C-----
C----- COMPUTE SLOPE, Y-INTERCEPT, VARIANCES, AND RESIDUALS
C----- BY WILLIAMSON PROCEDURE WITH INITIAL SLOPE VALUE BSTART
C-----
      CALL RWILL (N, X, Y, U, V, BSTART, MXITER, DELMAX,
     1            OUTPUT, XRES, YRES, WORK, IER)
C-----
C----- CHECK CORRECTNESS OF PROGRAM EXECUTION
C-----
       ERR = SQRT(EPS)
       IF (IER .EQ. 0                                  .AND.
     *     INT(OUTPUT(11)) .EQ. 4                      .AND.
     *     ABS(OUTPUT(2) - AT) / AT          .LE. ERR  .AND.
     *     ABS(OUTPUT(1) - BT) / BT          .LE. ERR  .AND.
     *     ABS(OUTPUT(6) - AJSDAT) / AJSDAT  .LE. ERR  .AND.
     *     ABS(OUTPUT(5) - AJSDBT) / AJSDBT  .LE. ERR  .AND.
     *     ABS(OUTPUT(4) - SDAT) / SDAT      .LE. ERR  .AND.
     *     ABS(OUTPUT(3) - SDBT) / SDBT      .LE. ERR  .AND.
     *     ABS(OUTPUT(7) - XBAR4T) / XBAR4T  .LE. ERR  .AND.
     *     ABS(OUTPUT(8) - YBAR4T) / YBAR4T  .LE. ERR  .AND.
     *     ABS(OUTPUT(9) - RMSXT) / RMSXT    .LE. ERR  .AND.
     *     ABS(OUTPUT(10) - RMSYT) / RMSYT   .LE. ERR) GO TO 10
C-----
C----- EXECUTION ERROR
C-----
      IF(KPRINT.GE.1) WRITE(LUN,900)
  900 FORMAT(' ** TEST 3 FOR SUBLIBRARY SLRPACK FAILED ** ')
      IF(KPRINT.GE.2) WRITE (LUN,901)
  901 FORMAT (' EXECUTION ERROR - SUBROUTINE RWILL  ')
      IF(KPRINT.GT.2) WRITE (LUN,902)
  902 FORMAT (' RESULTS DO NOT AGREE WITH KNOWN RESULTS')
      IF(KPRINT.GT.2) WRITE(LUN,904) OUTPUT(2), AT, OUTPUT(1),
     1 BT, OUTPUT(6), AJSDAT, OUTPUT(5), AJSDBT, OUTPUT(4), SDAT,
     2 OUTPUT(3), SDBT, OUTPUT(7), XBAR4T, OUTPUT(8), YBAR4T,
     3 OUTPUT(9), RMSXT, OUTPUT(10), RMSYT
 904  FORMAT(8X, 'RWILL RESULTS', 22X, 'EXPECTED RESULTS', /,
     1 10(3X, F20.15, 17X, F20.15, /))
      IF(KPRINT.GT.2) WRITE (LUN,905) INT(OUTPUT(11)), IER
 905  FORMAT(5X, 'NUMBER OF ITERATIONS', 5X, I3, 5X,
     1      'RWILL ERROR CONDITION', 5X, I3)
      IPASS=2
      GO TO 999
C-----
C----- PRINT CORRECT EXECUTION MESSAGE
C-----
 10   IF(KPRINT.GE.2) WRITE (LUN,903)
  903 FORMAT (' ** TEST 3 FOR SUBLIBRARY SLRPACK PASSED ** ')
  999 RETURN
      END
      subroutine slrqx4(kprint,ipass)
      common/unit/lun
      dimension x(7),y(7)
      integer p, q, r, tp, tq, tr, kprint, ipass
      real lambda
c
      data x /9.4,7.5,5.9,0.0,1.5,2.5,4.3/
      data y /2.8,3.3,3.1,1.0,1.8,1.3,2.2/
      ipass = 1
      n = 7
      tbeta1 = 1.2869565217391
      tbeta2 =  .2173913043478
      tlambd =  .5304347826087
      tp = 3
      tq = 1
      tr = 6
      iter = 1
      eps = sqrt(r2mach(4))
c
c     perform simple linear regression with l(infinity)-norm
c
      call linfs(n, x, y, beta1, beta2, lambda, kount, p, q, r, ifault)
c
c     check correctness of program execution for given data
c
      if (ifault    .eq. 0      .and.
     *    kount     .eq. iter   .and.
     *    p         .eq. tp     .and.
     *    q         .eq. tq     .and.
     *    r         .eq. tr     .and.
     *    abs((beta1  - tbeta1) / tbeta1) .le. eps .and.
     *    abs((beta2  - tbeta2) / tbeta2) .le. eps .and.
     *    abs((lambda - tlambd) / tlambd) .le. eps) go to 10
c-----
c----- execution error
c-----
      if(kprint.ge.1) write(lun,900)
  900 format(' ** TEST 4 FOR SUBLIBRARY SLRPACK FAILED ** ')
      if(kprint.ge.2) write (lun,901)
  901 format (' EXECUTION ERROR - ',/
     1 ' SUBROUTINE LINFS RESULTS DO NOT AGREE WITH KNOWN RESULTS')
      if(kprint.gt.2) write(lun,903) tbeta1, beta1, tbeta2, beta2,
     1    tlambd, lambda, tp, p, tq, q, tr, r, iter, kount, ifault
 903  format(1x,//,5x,'EXPECTED RESULTS',5x,'CALCULATED RESULTS',
     1       3(/,f20.13,3x,f20.13), 4(/,2x,i11,13x,i11), /,
     2       /, 3x, 'IFAULT = ', i3)
      ipass=2
      go to 999
c-----
c----- print correct execution message
c-----
   10 if(kprint.ge.2) write (lun,902)
  902 format (' ** TEST 4 FOR SUBLIBRARY SLRPACK PASSED ** ')
  999 return
      end
      subroutine slrqx5(kprint,ipass)
      common/unit/lun
      double precision x(7), y(7), tbeta1, tbeta2, tlambd
      double precision beta1, beta2, lambda, d2mach
      double precision eps
      double precision reler1, reler2, reler3
      integer p, q, r, tp, tq, tr, kprint, ipass
c
      data tbeta1 / 1.28695652173912254493756171D0 /
      data tbeta2 /  .217391304347826803319708391D0 /
      data tlambd /  .530434782608692395407775728D0 /
      data x /9.4,7.5,5.9,0.0,1.5,2.5,4.3/
      data y /2.8,3.3,3.1,1.0,1.8,1.3,2.2/
      ipass=1
      n = 7
      iter = 1
      tp = 3
      tq = 1
      tr = 6
      eps = max(10.0d0*sqrt(d2mach(4)),1.0d-14)
c
c     perform simple linear regression with l(infinity)-norm
c
      call dlinfs(n, x, y, beta1, beta2, lambda, kount, p, q, r, ifault)
c
c     check correctness of program execution for given data
c
      reler1 = dabs((beta1  - tbeta1) / tbeta1)
      reler2 = dabs((beta2  - tbeta2) / tbeta2)
      reler3 = dabs((lambda - tlambd) / tlambd)

      if (ifault    .eq. 0      .and.
     *    p         .eq. tp     .and.
     *    q         .eq. tq     .and.
     *    r         .eq. tr     .and.
     *    reler1    .le. eps    .and.
     *    reler2    .le. eps    .and.
     *    reler3    .le. eps) go to 10
c-----
c----- execution error
c-----
      if(kprint.ge.1) write(lun,900)
  900 format(' ** TEST 5 FOR SUBLIBRARY SLRPACK FAILED ** ')
      if(kprint.ge.2) write (lun,901)
  901 format (' EXECUTION ERROR - ',/
     1 ' SUBROUTINE DLINFS RESULTS DO NOT AGREE WITH KNOWN RESULTS')
c      if(kprint.gt.2) write(lun,903) tbeta1, beta1, tbeta2, beta2,
c     1   tlambd, lambda, tp, p, tq, q, tr, r, iter, kount, ifault
c 903  format(1x,//,3x,'EXPECTED RESULTS',6x,'CALCULATED RESULTS',
c     1       3(/,d20.13,3x,d20.13), 4(/,2x,i9,11x,i11), /,
c     2       3x, 'IFAULT = ', i3)
      if(kprint.gt.2) write(lun,903) eps, tbeta1, beta1, reler1, tbeta2, 
     1   beta2, reler2, tlambd, lambda, reler3, tp, p, tq, q, tr, r,  
     2   iter, kount, ifault
 903  format(1x,//,3x,'TOLERANCE = ',d14.6,
     +       //,3x,'EXPECTED RESULTS',8x,'CALCULATED RESULTS',
     +       8x,'REL ERROR',
     1       3(/,d22.15,3x,d22.15,3x,d14.6), 4(/,2x,i9,11x,i11), /,
     2       3x, 'IFAULT = ', i3)
      ipass=2
      go to 999
c-----
c----- print correct execution message
c-----
   10 if(kprint.ge.2) write (lun,902)
  902 format (' ** TEST 5 FOR SUBLIBRARY SLRPACK PASSED ** ')
  999 return
      end
