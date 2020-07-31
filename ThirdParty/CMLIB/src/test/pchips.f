C
C   DRIVER FOR TESTING SLATEC FORTRAN '77 PROGRAMS
C     PCHIP
C
      COMMON/UNIT/LUN
      COMMON/MSG/ICNT,JTEST(38)
      LUN=I1MACH(2)
      LIN=I1MACH(1)
      ITEST=1
C
C     READ KPRINT PARAMETER FROM DATA CARD..
C     KPRINT = 0   NO PRINTING
C              1   NO PRINTING FOR PASSED TESTS, SHORT MESSAGE
C                  FOR FAILED TESTS
C              2   PRINT SHORT MESSAGE FOR PASSED TESTS, FULLER
C                  INFORMATION FOR FAILED TESTS
C              3   PRINT COMPLETE QUICK-CHECK RESULTS
C
      READ(LIN,1) KPRINT
1     FORMAT(I1)
      CALL XSETUN(LUN)
      CALL XSETF(1)
      CALL XERMAX(1000)
C   TEST PCHIP PACKAGE
      CALL PCHQK1(KPRINT,IPAS)
      IPASS=0
      IF(IPAS.EQ.0) IPASS=1
      ITEST=ITEST*IPASS
      CALL PCHQK2(KPRINT,IPAS)
      IPASS=0
      IF(IPAS.EQ.0) IPASS=1
      ITEST=ITEST*IPASS
C
      IF(KPRINT.GE.1.AND.ITEST.NE.1) WRITE(LUN,2)
2     FORMAT(/' ***** WARNING -- AT LEAST ONE TEST FOR THE PCHIP,
     1 SUBLIBRARY FAILED ***** ')
      IF(KPRINT.GE.1.AND.ITEST.EQ.1) WRITE(LUN,3)
3     FORMAT(/' ----- THE PCHIP SUBLIBRARY PASSED ALL TESTS ----- ')
      END
      SUBROUTINE PCHQK1 ( KPRINT, IPASS)
      COMMON/UNIT/LUN
      INTEGER  LUN, KPRINT, IPASS
C
C-----------------------------------------------------------------------
C
C              PCHIP QUICKCHECK NUMBER 1
C
C     TESTS THE EVALUATORS:  CHFDV, CHFEV, PCHFD, PCHFE.
C
C-----------------------------------------------------------------------
C
C     RETURNED VALUES:
C        IPASS = 0  IF ALL TESTS PASSED.
C        IPASS BETWEEN 1 AND 7 IS THE SUM OF:
C           IPASS = 1  IF SINGLE CUBIC TEST FAILED. (SEE EVCHCK OUTPUT.)
C           IPASS = 2  IF PCHFD/PCHFE  TEST FAILED. (SEE EVPCCK OUTPUT.)
C           IPASS = 4  IF ERROR RETURN TEST FAILED. (SEE EVERCK OUTPUT.)
C
C-----------------------------------------------------------------------
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I1, I2, I3, I4, I5, I6, I7, I8, I9, NPTS
      REAL  WORK (4000)
      LOGICAL  FAIL
C
      IPASS = 0
C
C  TEST CHFDV AND CHFEV.
C
      NPTS = 1000
      I1 = 1  + NPTS
      I2 = I1 + NPTS
      I3 = I2 + NPTS
      CALL EVCHCK (LUN, KPRINT, NPTS, WORK(1), WORK(I1), WORK(I2),
     *                                          WORK(I3), FAIL)
      IF (FAIL)  IPASS = IPASS + 1
C
C  TEST PCHFD AND PCHFE.
C
      I1 = 1  +  10
      I2 = I1 +  10
      I3 = I2 + 100
      I4 = I3 + 100
      I5 = I4 + 100
      I6 = I5 +  51
      I7 = I6 +  51
      I8 = I7 +  51
      I9 = I8 +  51
      CALL EVPCCK (LUN, KPRINT, WORK(1), WORK(I1), WORK(I2), WORK(I3),
     *             WORK(I4), WORK(I5), WORK(I6), WORK(I7), WORK(I8),
     *             WORK(I9), FAIL)
      IF (FAIL)  IPASS = IPASS + 2
C
C  TEST ERROR RETURNS.
C
C   WHEN KPRINT .LE. 1
C  BYPASS ERROR CHECKS THAT ARE FATAL ON A VAX 11/780
C
      FAIL=.FALSE.
      IF(KPRINT.GT.1)
     1 CALL EVERCK (LUN, KPRINT, FAIL)
      IF (FAIL)  IPASS = IPASS + 4
C
      IF(KPRINT.GE.1.AND.IPASS.NE.0) WRITE(LUN,99999)
      IF(KPRINT.GE.2.AND.IPASS.EQ.0) WRITE(LUN,99998)
99998 FORMAT(/' ***** THE PCHIP SUBLIBRARY PASSED ALL TESTS ***** ')
99999 FORMAT(/' ***** THE PCHIP SUBLIBRARY FAILED SOME TESTS ***** ')
      RETURN
      END
      SUBROUTINE  EVCHCK (LUN, KPRINT, NPTS, XEV, FEV, DEV, FEV2, FAIL)
      INTEGER  LUN, KPRINT, NPTS
      REAL  XEV(NPTS), FEV(NPTS), DEV(NPTS), FEV2(NPTS)
      LOGICAL  FAIL
C
C -------- CODE TO TEST EVALUATION ACCURACY OF CHFDV AND CHFEV --------
C
C     USING FUNCTION AND DERIVATIVE VALUES FROM A CUBIC (COMPUTED IN
C     DOUBLE PRECISION) AT NINT DIFFERENT (X1,X2) PAIRS:
C     1. CHECKS THAT CHFDV AND CHFEV BOTH REPRODUCE ENDPOINT VALUES.
C     2. EVALUATES AT NPTS POINTS, 10 OF WHICH ARE OUTSIDE THE INTERVAL
C        AND:
C        A. CHECKS ACCURACY OF CHFDV FUNCTION AND DERIVATIVE VALUES
C           AGAINST EXACT VALUES.
C        B. CHECKS THAT RETURNED VALUES OF NEXT SUM TO 10.
C        C. CHECKS THAT FUNCTION VALUES FROM CHFEV AGREE WITH THOSE
C           FROM CHFDV.
C
C
C     FORTRAN INTRINSICS USED:  ABS, AMAX1, AMIN1, FLOAT.
C     FORTRAN LIBRARY ROUTINES USED:  SQRT, (READ), (WRITE).
C     SLATEC LIBRARY ROUTINES USED:  CHFDV, CHFEV, R1MACH, RAND.
C     OTHER ROUTINES USED:  FDTRUE.
C
C-----------------------------------------------------------------------
C
C     PROGRAMMED BY:  F. N. FRITSCH, LLNL.
C
C     CHANGE RECORD:
C       82-06-24/9 CONVERTED TO QUICKCHECK FOR SLATEC LIBRARY.
C       82-06-30   1. MODIFIED DEFINITIONS OF RELATIVE ERROR AND TEST
C                     TOLERANCES.
C                  2. VARIOUS IMPROVEMENTS TO OUTPUT FORMATS.
C       82-07-16   1. SET MACHEP VIA A CALL TO R1MACH.
C                  2. CHANGED FROM FORTLIB'S RANF TO SLATEC'S RAND.
C
C-----------------------------------------------------------------------
C
C  DECLARATIONS.
C
      EXTERNAL RAND
      INTEGER  I, IERR, IINT, NEXT(2), NEXT2(2), NINT
      REAL  AED, AED2, AEDMAX, AEDMIN, AEF, AEF2, AEFMAX, AEFMIN,
     *      CHECK(2), CHECKF(2), CHECKD(2), D1, D2, DERMAX, DTRUE, DX,
     *      EPS1, EPS2, F1, F2, FACT, FERMAX, FLOORD, FLOORF, FOUR,
     *      FTRUE, LEFT(3), MACHEP,
     *      ONE, RED, RED2, REDMAX, REDMIN, REF, REF2, REFMAX,
     *      REFMIN, RIGHT(3), SMALL, TEN, TOL1, TOL2,
     *      X1, X2, XADMAX, XADMIN, XAFMAX, XAFMIN, XRDMAX,
     *      XRDMIN, XRFMAX, XRFMIN, ZERO
      LOGICAL  FAILOC, FAILNX
      REAL  RAND
C
C  DEFINE RELATIVE ERROR WITH FLOOR.
C
      RERR(ERR,VALUE,FLOOR) = ERR / AMAX1(ABS(VALUE), FLOOR)
C
C  INITIALIZE.
C
      DATA  ZERO /0./,  ONE /1./,  FOUR /4./,  TEN /10./
      DATA  SMALL  /1.0E-10/
      DATA  NINT /3/
      DATA   LEFT /-1.5, 2.0E-10, 1.0   /
      DATA  RIGHT / 2.5, 3.0E-10, 1.0E+8/
C
      MACHEP = R1MACH(4)
      EPS1 = FOUR*MACHEP
      EPS2 = TEN*MACHEP
C
      FAIL = .FALSE.
C
      IF (KPRINT .EQ. 2)  WRITE (LUN, 3000)
 3000 FORMAT ('1'//10X,'EVCHCK RESULTS'/10X,'--------------')
C
C  CYCLE OVER INTERVALS.
C
      DO 90  IINT = 1, NINT
      X1 =  LEFT(IINT)
      X2 = RIGHT(IINT)
C
      FACT = AMAX1(SQRT(X2-X1), ONE)
      TOL1 = EPS1 * FACT
      TOL2 = EPS2 * FACT
C
C  COMPUTE AND PRINT ENDPOINT VALUES.
C
      CALL FDTRUE (X1, F1, D1)
      CALL FDTRUE (X2, F2, D2)
C
      IF (KPRINT .EQ. 3)  THEN
         WRITE (LUN, 2000)
 2000    FORMAT ('1'/30X,'CHFDV ACCURACY TEST'///)
         WRITE (LUN, 2001)  'X1', X1, 'X2', X2
         WRITE (LUN, 2001)  'F1', F1, 'F2', F2
         WRITE (LUN, 2001)  'D1', D1, 'D2', D2
 2001    FORMAT (/10X,A2,' =',E18.10,5X,A2,' =',E18.10)
      ENDIF
C
      IF (KPRINT .EQ. 2)  WRITE (LUN, 3001)  X1, X2
 3001 FORMAT (//10X,'INTERVAL = (',E12.5,',',E12.5,' ):' )
C
C  COMPUTE FLOORS FOR RELATIVE ERRORS.
C
      FLOORF = AMAX1( AMIN1(ABS(F1),ABS(F2)), SMALL)
      FLOORD = AMAX1( AMIN1(ABS(D1),ABS(D2)), SMALL)
C
C  CHECK REPRODUCTION OF ENDPOINT VALUES.
C
      XEV(1) = X1
      XEV(2) = X2
C     -----------------------------------------------------------
      CALL CHFDV (X1, X2, F1, F2, D1, D2, 2, XEV, CHECKF, CHECKD,
     *            NEXT, IERR)
C     -----------------------------------------------------------
      AEF  = CHECKF(1)-F1
      REF  = RERR(AEF , F1, FLOORF)
      AEF2 = CHECKF(2)-F2
      REF2 = RERR(AEF2, F2, FLOORF)
      AED  = CHECKD(1)-D1
      RED  = RERR(AED , D1, FLOORD)
      AED2 = CHECKD(2)-D2
      RED2 = RERR(AED2, D2, FLOORD)
C
      FAILOC = AMAX1(ABS(REF),ABS(REF2),ABS(RED),ABS(RED2)) .GT. TOL1
      FAIL = FAIL .OR. FAILOC
C
      IF (KPRINT .EQ. 3)  THEN
         WRITE (LUN, 2002)  NEXT, AEF, AEF2, AED, AED2
 2002    FORMAT (//5X,'ERRORS AT ENDPOINTS:',40X,'(NEXT =',2I3,')'
     *           // 7X,'F1:',E13.5,5X,'F2:',E13.5,
     *              5X,'D1:',E13.5,5X,'D2:',E13.5)
         WRITE (LUN, 2003)  REF, REF2, RED, RED2
 2003    FORMAT (2X,4(8X,E13.5))
      ENDIF
C
      IF (FAILOC .AND. (KPRINT.GE.2))  WRITE (LUN, 3002)
 3002 FORMAT (/' ***** CHFDV FAILED TO REPRODUCE ENDPOINT VALUES.')
C
C  CHFEV SHOULD AGREE EXACTLY WITH CHFDV.
C                     -------
C     --------------------------------------------------------------
      CALL CHFEV (X1, X2, F1, F2, D1, D2, 2, XEV, CHECK, NEXT, IERR)
C     --------------------------------------------------------------
      FAILOC = (CHECK(1).NE.CHECKF(1)) .OR. (CHECK(2).NE.CHECKF(2))
      FAIL = FAIL .OR. FAILOC
C
      IF (FAILOC .AND. (KPRINT.GE.2))  WRITE (LUN, 3003)
 3003 FORMAT (/' ***** CHFEV DOES NOT AGREE WITH CHFDV AT ENDPOINTS.')
C
C  EVALUATE AT NPTS 'UNIFORMLY RANDOM' POINTS IN (X1,X2).
C     THIS VERSION EXTENDS EVALUATION DOMAIN BY ADDING 4 SUBINTERVALS
C     TO LEFT AND 6 TO RIGHT OF [X1,X2].
C
      DX = (X2-X1)/FLOAT(NPTS-10)
      DO 20  I = 1, NPTS
         XEV(I) = (X1 + (I-5)*DX) + DX*RAND(ZERO)
   20 CONTINUE
C     --------------------------------------------------------
      CALL CHFDV (X1, X2, F1, F2, D1, D2, NPTS, XEV, FEV, DEV,
     *            NEXT, IERR)
C     --------------------------------------------------------
      IF (IERR .NE. 0)  THEN
         FAILOC = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LUN, 4003)  IERR
 4003    FORMAT (/' ***** ERROR ***** CHFDV RETURNED IERR =',I5)
      ELSE
C
C     CUMULATE LARGEST AND SMALLEST ERRORS FOR SUMMARY.
C
      DO 30  I = 1, NPTS
         CALL  FDTRUE (XEV(I), FTRUE, DTRUE)
         AEF = FEV(I) - FTRUE
         REF = RERR(AEF, FTRUE, FLOORF)
         AED = DEV(I) - DTRUE
         RED = RERR(AED, DTRUE, FLOORD)
C
         IF (I .EQ. 1)  THEN
C            INITIALIZE.
            AEFMIN = AEF
            AEFMAX = AEF
            AEDMIN = AED
            AEDMAX = AED
            REFMIN = REF
            REFMAX = REF
            REDMIN = RED
            REDMAX = RED
            XAFMIN = XEV(1)
            XAFMAX = XEV(1)
            XADMIN = XEV(1)
            XADMAX = XEV(1)
            XRFMIN = XEV(1)
            XRFMAX = XEV(1)
            XRDMIN = XEV(1)
            XRDMAX = XEV(1)
         ELSE
C            SELECT.
            IF (AEF .LT. AEFMIN)  THEN
               AEFMIN = AEF
               XAFMIN = XEV(I)
            ELSE IF (AEF .GT. AEFMAX)  THEN
               AEFMAX = AEF
               XAFMAX = XEV(I)
            ENDIF
            IF (AED .LT. AEDMIN)  THEN
               AEDMIN = AED
               XADMIN = XEV(I)
            ELSE IF (AED .GT. AEDMAX)  THEN
               AEDMAX = AED
               XADMAX = XEV(I)
            ENDIF
            IF (REF .LT. REFMIN)  THEN
               REFMIN = REF
               XRFMIN = XEV(I)
            ELSE IF (REF .GT. REFMAX)  THEN
               REFMAX = REF
               XRFMAX = XEV(I)
            ENDIF
            IF (RED .LT. REDMIN)  THEN
               REDMIN = RED
               XRDMIN = XEV(I)
            ELSE IF (RED .GT. REDMAX)  THEN
               REDMAX = RED
               XRDMAX = XEV(I)
            ENDIF
         ENDIF
   30    CONTINUE
C
         FERMAX = AMAX1 (ABS(REFMAX), ABS(REFMIN))
         DERMAX = AMAX1 (ABS(REDMAX), ABS(REDMIN))
C
         FAILNX = (NEXT(1) + NEXT(2)) .NE. 10
         FAILOC = FAILNX .OR. (AMAX1(FERMAX, DERMAX) .GT. TOL2)
      ENDIF
      FAIL = FAIL .OR. FAILOC
C
C  PRINT SUMMARY.
C
      IF (KPRINT .EQ. 3)  THEN
         WRITE (LUN, 2004)  NPTS-10, NEXT
 2004    FORMAT (//5X,'ERRORS AT ',I5,' INTERIOR POINTS + 10 OUTSIDE:',
     *                   15X,'(NEXT =',2I3,')'
     *           //35X,'FUNCTION',17X,'DERIVATIVE'
     *            /20X,2(11X,'ABS',9X,'REL') )
C
         WRITE (LUN, 2005)  'MIN', AEFMIN, REFMIN, AEDMIN, REDMIN
 2005    FORMAT (/10X,A3,'IMUM ERROR:  ',2E12.4,2X,2E12.4)
         WRITE (LUN, 2006) XAFMIN, XRFMIN, XADMIN, XRDMIN
 2006    FORMAT ( 10X,'LOCATED AT X =  ',2E12.4,2X,2E12.4)
         WRITE (LUN, 2005)  'MAX', AEFMAX, REFMAX, AEDMAX, REDMAX
         WRITE (LUN, 2006) XAFMAX, XRFMAX, XADMAX, XRDMAX
      ENDIF
C
      IF (KPRINT .GE. 2)  THEN
         IF (FAILOC) THEN
            IF (FERMAX .GT. TOL2)  WRITE (LUN, 3006) 'F', FERMAX, TOL2
 3006       FORMAT (/' ***** MAXIMUM RELATIVE ERROR IN ',A1,' =',
     *                 E12.5,', WHICH EXCEEDS TOLERANCE ',E12.5)
            IF (DERMAX .GT. TOL2)  WRITE (LUN, 3006) 'D', DERMAX, TOL2
            IF (FAILNX)  WRITE (LUN, 4006)  NEXT
 4006       FORMAT (/' ***** REPORTED NEXT =',2I5,
     *                   '   RATHER THAN    4    6')
         ELSE
            WRITE (LUN, 5006)
 5006       FORMAT (/' CHFDV RESULTS OK.')
         ENDIF
      ENDIF
C
C  CHECK THAT CHFEV AGREES WITH CHFDV.
C
C     -----------------------------------------------------------------
      CALL CHFEV (X1, X2, F1, F2, D1, D2, NPTS, XEV, FEV2, NEXT2, IERR)
C     -----------------------------------------------------------------
      IF (IERR .NE. 0)  THEN
         FAILOC = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LUN, 3007)  IERR
 3007    FORMAT (/' ***** ERROR ***** CHFEV RETURNED IERR =',I5)
      ELSE
         AEFMAX = ABS(FEV2(1) - FEV(1))
         XAFMAX = XEV(1)
         DO 40  I = 2, NPTS
            AEF = ABS(FEV2(I) - FEV(I))
            IF (AEF .GT. AEFMAX)  THEN
               AEFMAX = AEF
               XAFMAX = XEV(I)
            ENDIF
   40    CONTINUE
         FAILNX = (NEXT2(1).NE.NEXT(1)) .OR. (NEXT2(2).NE.NEXT(2))
         FAILOC = FAILNX .OR. (AEFMAX.NE.ZERO)
         IF (KPRINT .GE. 2)  THEN
            IF (FAILOC)  THEN
               WRITE (LUN, 3008)
 3008          FORMAT (/' ***** CHFEV DID NOT AGREE WITH CHFDV:')
               IF (AEFMAX.NE.ZERO)  WRITE (LUN, 3009)  AEFMAX, XAFMAX
 3009          FORMAT ( 7X,'MAXIMUM DIFFERENCE OF ',E12.5,
     *                           ' OCCURRED AT X =',E12.5)
               IF (FAILNX)  WRITE (LUN, 4009)  NEXT2, NEXT
 4009          FORMAT ( 7X,'REPORTED NEXT =',2I3,'   RATHER THAN ',2I3)
            ELSE
               WRITE (LUN, 5009)
 5009          FORMAT (/' CHFEV AGREES WITH CHFDV.')
            ENDIF
         ENDIF
      ENDIF
C
      FAIL = FAIL .OR. FAILOC
C
C  GO BACK FOR ANOTHER INTERVAL.
C
   90 CONTINUE
C
      RETURN
C
      END
      SUBROUTINE  FDTRUE (X, F, D)
      REAL  X, F, D
C
C        COMPUTE EXACT FUNCTION VALUES IN DOUBLE PRECISION.
C
C                   F(X) = X*(X+1)*(X-2)
C
      DOUBLE PRECISION  FACT1, FACT2, ONE, TWO, XX
      DATA  ONE /1.0D0/,  TWO /2.0D0/
C
      XX = X
      FACT1 = XX + ONE
      FACT2 = XX - TWO
      F = XX * FACT1 * FACT2
      D = FACT1*FACT2 + XX*(FACT1 + FACT2)
C
      RETURN
      END
      SUBROUTINE  EVPCCK (LUN, KPRINT, X, Y, F, FX, FY,
     *                    XE, YE, FE, DE, FE2, FAIL)
      INTEGER  LUN, KPRINT
      LOGICAL  FAIL
      REAL  X(10), Y(10), F(10,10), FX(10,10), FY(10,10),
     *      XE(51), YE(51), FE(51), DE(51), FE2(51)
C
C ---- CODE TO TEST USAGE OF INCREMENT ARGUMENT IN PCHFD AND PCHFE ----
C
C     EVALUATES A BICUBIC FUNCTION AND ITS FIRST PARTRIAL DERIVATIVES
C     ON A 4X6 MESH CONTAINED IN A 10X10 ARRAY.
C
C     INTERPOLATION OF THESE DATA ALONG MESH LINES IN EITHER DIMENSION
C     SHOULD AGRRE WITH CORRECT FUNCTION WITHIN ROUNDOFF ERROR.
C
C     ARRAYS ARE ARGUMENTS ONLY TO ALLOW SHARING STORAGE WITH OTHER
C     TEST ROUTINES.
C
C
C     FORTRAN INTRINSICS USED:  ABS, FLOAT.
C     FORTRAN LIBRARY ROUTINES USED:  (WRITE).
C     SLATEC LIBRARY ROUTINES USED:  PCHFD, PCHFE, R1MACH.
C
C-----------------------------------------------------------------------
C
C     PROGRAMMED BY:  F. N. FRITSCH, LLNL.
C
C     CHANGE RECORD:
C       82-06-30/
C          07-14   CONVERTED TO QUICKCHECK FOR SLATEC LIBRARY.
C       82-07-15   1. CORRECTED SOME FORMATS.
C                  2. ADDED CALL TO R1MACH TO SET MACHEP.
C
C-----------------------------------------------------------------------
C
C  DECLARATIONS.
C
      INTEGER  I, IER2, IERR, INC, J, K, NE, NERR, NMAX, NX, NY
      LOGICAL  FAILD, FAILE, FAILOC, SKIP
      REAL  DERMAX, DERR, DTRUE, DX, FDIFF, FDIFMX, FERMAX, FERR, FTRUE,
     *      MACHEP, TOL, PDERMX, PDIFMX, PFERMX, ZERO
C
C  DEFINE TEST FUNCTION AND DERIVATIVES.
C
      FCN (XDUM,YDUM) =  XDUM*(YDUM*YDUM)*(XDUM*XDUM + 1.)
      DFDX(XDUM,YDUM) = (YDUM*YDUM)*(3.*XDUM*XDUM + 1.)
      DFDY(XDUM,YDUM) =   2.*XDUM*YDUM*(XDUM*XDUM + 1.)
C
C  INITALIZE.
C
      DATA  NMAX /10/,  NX /4/,  NY /6/
      DATA  NE /51/
C
      DATA  ZERO /0./
C
      MACHEP = R1MACH(4)
      TOL = 10.*MACHEP
C
      FAIL = .FALSE.
C
C  SET UP 4-BY-6 MESH IN A 10-BY-10 ARRAY:
C     X =  0.25(0.25)1.   ;
C     Y = -0.75(0.5 )1.75 .
C
      DO 1  I = 1, NX
         X(I) = 0.25*FLOAT(I)
    1 CONTINUE
      DO 5  J = 1, NY
         Y(J) = 0.5*FLOAT(J) - 1.25
         DO 4  I = 1, NX
             F(I,J) = FCN (X(I), Y(J))
            FX(I,J) = DFDX(X(I), Y(J))
            FY(I,J) = DFDY(X(I), Y(J))
    4    CONTINUE
    5 CONTINUE
C
C  SET UP EVALUATION POINTS:
C     XE =  0.(0.02)1. ;
C     YE = -2.(0.08)2. .
C
      DX = 1./FLOAT(NE-1)
      DO 8  K = 1, NE
         XE(K) = DX*FLOAT(K-1)
         YE(K) = 4.*XE(K) - 2.
    8 CONTINUE
C
      IF (KPRINT .EQ. 2)  WRITE (LUN, 2000)
C
C  EVALUATE ON HORIZONTAL MESH LINES (Y FIXED, X RUNNING) ..............
C
      NERR = 0
      INC = 1
      SKIP = .FALSE.
      DO 20  J = 1, NY
C        --------------------------------------------------------------
         CALL PCHFD (NX, X, F(1,J), FX(1,J), INC, SKIP, NE, XE, FE, DE,
     *               IERR)
C        --------------------------------------------------------------
         IF (KPRINT .EQ. 3)
     *       WRITE (LUN, 2001)  INC, 'J', J, 'Y', Y(J), IERR, 'X'
         IF (IERR .LT. 0)  GO TO 15
C
C        PCHFE SHOULD AGREE EXACTLY WITH PCHFD.
C
C        -----------------------------------------------------------
         CALL PCHFE (NX, X, F(1,J), FX(1,J), INC, SKIP, NE, XE, FE2,
     *               IER2)
C        -----------------------------------------------------------
C
         DO 10  K = 1, NE
            FTRUE =  FCN(XE(K), Y(J))
            FERR = FE(K) - FTRUE
            DTRUE = DFDX(XE(K), Y(J))
            DERR = DE(K) - DTRUE
            IF (KPRINT .EQ. 3)
     *         WRITE (LUN, 2002)  XE(K), FTRUE, FE(K), FERR,
     *                                    DTRUE, DE(K), DERR
            IF (K .EQ. 1)  THEN
C              INITIALIZE.
               FERMAX = ABS(FERR)
               PFERMX = XE(1)
               DERMAX = ABS(DERR)
               PDERMX = XE(1)
               FDIFMX = ABS(FE2(1) - FE(1))
               PDIFMX = XE(1)
            ELSE
C              SELECT.
               FERR = ABS(FERR)
               IF (FERR .GT. FERMAX)  THEN
                  FERMAX = FERR
                  PFERMX = XE(K)
               ENDIF
               DERR = ABS(DERR)
               IF (DERR .GT. DERMAX)  THEN
                  DERMAX = DERR
                  PDERMX = XE(K)
               ENDIF
               FDIFF = ABS(FE2(K) - FE(K))
               IF (FDIFF .GT. FDIFMX)  THEN
                  FDIFMX = FDIFF
                  PDIFMX = XE(K)
               ENDIF
            ENDIF
   10    CONTINUE
C
         FAILD = (FERMAX.GT.TOL) .OR. (DERMAX.GT.TOL)
         FAILE = FDIFMX .NE. ZERO
         FAILOC = FAILD .OR. FAILE .OR. (IERR.NE.13) .OR. (IER2.NE.IERR)
C
         IF (FAILOC .AND. (KPRINT.EQ.2))
     *      WRITE (LUN, 2003)  'J', J, 'Y', Y(J)
C
         IF ((KPRINT.EQ.3) .OR. (FAILD.AND.(KPRINT.EQ.2)) )
     *      WRITE (LUN, 2004)  FERMAX, PFERMX, DERMAX, PDERMX
C
         IF ((KPRINT.EQ.3) .OR. (FAILE.AND.(KPRINT.EQ.2)) )
     *      WRITE (LUN, 2005)  FDIFMX, PDIFMX
C
         IF ((IERR.NE.13) .AND. (KPRINT.GE.2))
     *      WRITE (LUN, 2006)  'D', IERR, 13
C
         IF ((IER2.NE.IERR) .AND. (KPRINT.GE.2))
     *      WRITE (LUN, 2006)  'E', IER2, IERR
         GO TO 19
C
   15    CONTINUE
         FAILOC = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LUN, 3000) IERR
C
   19    CONTINUE
         IF (FAILOC)  NERR = NERR + 1
         FAIL = FAIL .OR. FAILOC
   20 CONTINUE
C
      IF (KPRINT .GE. 2)  THEN
         IF (NERR .GT. 0)  THEN
            WRITE (LUN, 3001)  NERR, 'J'
         ELSE
            WRITE (LUN, 4000)  'J'
         ENDIF
      ENDIF
C
C  EVALUATE ON VERTICAL MESH LINES (X FIXED, Y RUNNING) ................
C
      NERR = 0
      INC = NMAX
      SKIP = .FALSE.
      DO 40  I = 1, NX
C        --------------------------------------------------------------
         CALL PCHFD (NY, Y, F(I,1), FY(I,1), INC, SKIP, NE, YE, FE, DE,
     *               IERR)
C        --------------------------------------------------------------
         IF (KPRINT .EQ. 3)
     *       WRITE (LUN, 2001)  INC, 'I', I, 'X', X(I), IERR, 'Y'
         IF (IERR .LT. 0)  GO TO 35
C
C        PCHFE SHOULD AGREE EXACTLY WITH PCHFD.
C
C        -----------------------------------------------------------
         CALL PCHFE (NY, Y, F(I,1), FY(I,1), INC, SKIP, NE, YE, FE2,
     *               IER2)
C        -----------------------------------------------------------
C
         DO 30  K = 1, NE
            FTRUE =  FCN(X(I), YE(K))
            FERR = FE(K) - FTRUE
            DTRUE = DFDY(X(I), YE(K))
            DERR = DE(K) - DTRUE
            IF (KPRINT .EQ. 3)
     *         WRITE (LUN, 2002)  YE(K), FTRUE, FE(K), FERR,
     *                                    DTRUE, DE(K), DERR
            IF (K .EQ. 1)  THEN
C              INITIALIZE.
               FERMAX = ABS(FERR)
               PFERMX = YE(1)
               DERMAX = ABS(DERR)
               PDERMX = YE(1)
               FDIFMX = ABS(FE2(1) - FE(1))
               PDIFMX = YE(1)
            ELSE
C              SELECT.
               FERR = ABS(FERR)
               IF (FERR .GT. FERMAX)  THEN
                  FERMAX = FERR
                  PFERMX = YE(K)
               ENDIF
               DERR = ABS(DERR)
               IF (DERR .GT. DERMAX)  THEN
                  DERMAX = DERR
                  PDERMX = YE(K)
               ENDIF
               FDIFF = ABS(FE2(K) - FE(K))
               IF (FDIFF .GT. FDIFMX)  THEN
                  FDIFMX = FDIFF
                  PDIFMX = YE(K)
               ENDIF
            ENDIF
   30    CONTINUE
C
         FAILD = (FERMAX.GT.TOL) .OR. (DERMAX.GT.TOL)
         FAILE = FDIFMX .NE. ZERO
         FAILOC = FAILD .OR. FAILE .OR. (IERR.NE.20) .OR. (IER2.NE.IERR)
C
         IF (FAILOC .AND. (KPRINT.EQ.2))
     *      WRITE (LUN, 2003)  'I', I, 'X', X(I)
C
         IF ((KPRINT.EQ.3) .OR. (FAILD.AND.(KPRINT.EQ.2)) )
     *      WRITE (LUN, 2004)  FERMAX, PFERMX, DERMAX, PDERMX
C
         IF ((KPRINT.EQ.3) .OR. (FAILE.AND.(KPRINT.EQ.2)) )
     *      WRITE (LUN, 2005)  FDIFMX, PDIFMX
C
         IF ((IERR.NE.20) .AND. (KPRINT.GE.2))
     *      WRITE (LUN, 2006)  'D', IERR, 20
C
         IF ((IER2.NE.IERR) .AND. (KPRINT.GE.2))
     *      WRITE (LUN, 2006)  'E', IER2, IERR
         GO TO 39
C
   35    CONTINUE
         FAILOC = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LUN, 3000) IERR
C
   39    CONTINUE
         IF (FAILOC)  NERR = NERR + 1
         FAIL = FAIL .OR. FAILOC
   40 CONTINUE
C
      IF (KPRINT .GE. 2)  THEN
         IF (NERR .GT. 0)  THEN
            WRITE (LUN, 3001)  NERR, 'I'
         ELSE
            WRITE (LUN, 4000)  'I'
         ENDIF
      ENDIF
C
C  TERMINATE.
C
      RETURN
C
C  FORMATS.
C
 2000 FORMAT ('1'//10X,'EVPCCK RESULTS'/10X,'--------------')
 2001 FORMAT ('1',20X,'PCHFD INCREMENT TEST -- INCFD = ',I2
     *        //15X,'ON ',A1,'-LINE ',I2,',  ',A1,' =',F8.4,
     *           '  --  IERR =',I3
     *        //3X,A1,'E',10X,'F',8X,'FE',9X,'DIFF',
     *                    13X,'D',8X,'DE',9X,'DIFF')
 2002 FORMAT (F7.2,2(2X,2F10.5,E15.5))
 2003 FORMAT (/' ***** PCHFD AND/OR PCHFE FAILED ON ',A1,'-LINE ',I1,
     *                             ',  ',A1,' =',F8.4)
 2004 FORMAT (/'  MAXIMUM ERROR IN FUNCTION =',E13.5,' (AT',F6.2,
     *                    '), IN DERIVATIVE =',E13.5,' (AT',F6.2,').' )
 2005 FORMAT (/'  MAXIMUM DIFFERENCE BETWEEN PCHFE AND PCHFD =',
     *                                         E13.5,' (AT',F6.2,').' )
 2006 FORMAT (/'  PCHF',A1,' RETURNED IERR = ',I2,' INSTEAD OF ',I2)
 3000 FORMAT (//' ***** ERROR ***** PCHFD RETURNED IERR =',I5//)
 3001 FORMAT (//' ***** ERROR ***** PCHFD AND/OR PCHFE FAILED ON',I2,1X,
     *                                   A1,'-LINES.'//)
 4000 FORMAT (/' PCHFD AND PCHFE OK ON ',A1,'-LINES.')
C
      END
      SUBROUTINE  EVERCK (LUN, KPRINT, FAIL)
      INTEGER  LUN, KPRINT
      LOGICAL  FAIL
C
C --------- CODE TO TEST ERROR RETURNS FROM PCHIP EVALUATORS. ---------
C
C
C     FORTRAN LIBRARY ROUTINES USED:  (WRITE).
C     SLATEC LIBRARY ROUTINES USED:  CHFDV, CHFEV, PCHFD, PCHFE.
C
C-----------------------------------------------------------------------
C
C     PROGRAMMED BY:  F. N. FRITSCH, LLNL.
C
C     CHANGE RECORD:
C       82-07-15   CONVERTED TO QUICKCHECK FOR SLATEC LIBRARY.
C
C-----------------------------------------------------------------------
C
C  DECLARATIONS.
C
      INTEGER  I, IERR, N, NERR, NEXT(2)
      REAL  D(10), DUM, F(10), TEMP, X(10)
      LOGICAL  COMP, SKIP
C
C  INITIALIZE.
C
      DATA  N /10/
      NERR = 0
C
      IF (KPRINT .GE. 2)  WRITE (LUN, 5000)
C
C  FIRST, TEST CHFEV AND CHFDV.
C
      IF (KPRINT .EQ. 3)  WRITE (LUN, 5001)  (-1)
      CALL CHFEV (0., 1., 3., 7., 3., 6., 0, DUM, DUM, NEXT, IERR)
      IF (.NOT. COMP (IERR, -1, LUN, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .EQ. 3)  WRITE (LUN, 5001)  (-2)
      CALL CHFEV (1., 1., 3., 7., 3., 6., 1, DUM, DUM, NEXT, IERR)
      IF (.NOT. COMP (IERR, -2, LUN, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .EQ. 3)  WRITE (LUN, 5001)  (-1)
      CALL CHFDV (0., 1., 3., 7., 3., 6., 0, DUM, DUM, DUM, NEXT, IERR)
      IF (.NOT. COMP (IERR, -1, LUN, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .EQ. 3)  WRITE (LUN, 5001)  (-2)
      CALL CHFDV (1., 1., 3., 7., 3., 6., 1, DUM, DUM, DUM, NEXT, IERR)
      IF (.NOT. COMP (IERR, -2, LUN, KPRINT) )  NERR = NERR + 1
C
C  SET UP PCH DEFINITION.
C
 9999 CONTINUE
      DO 10  I = 1, N
         X(I) = I
         F(I) = I + 2
         D(I) = 1.
   10 CONTINUE
C
C  SWAP POINTS 4 AND 7, SO X-ARRAY IS OUT OF ORDER.
C
      TEMP = X(4)
      X(4) = X(7)
      X(7) = TEMP
C
C  NOW, TEST PCHFE AND PCHFD.
C
      IF (KPRINT .EQ. 3)  WRITE (LUN, 5001)  (-1)
      SKIP = .FALSE.
      CALL PCHFE (1, X, F, D, 0, SKIP, 0, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -1, LUN, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .EQ. 3)  WRITE (LUN, 5001)  (-2)
      SKIP = .FALSE.
      CALL PCHFE (4, X, F, D, 0, SKIP, 0, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -2, LUN, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .EQ. 3)  WRITE (LUN, 5001)  (-3)
      SKIP = .FALSE.
      CALL PCHFE (N, X, F, D, 1, SKIP, 0, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -3, LUN, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .EQ. 3)  WRITE (LUN, 5001)  (-4)
      SKIP = .TRUE.
      CALL PCHFE (N, X, F, D, 1, SKIP, 0, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -4, LUN, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .EQ. 3)  WRITE (LUN, 5001)  (-1)
      SKIP = .FALSE.
      CALL PCHFD (1, X, F, D, 0, SKIP, 0, DUM, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -1, LUN, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .EQ. 3)  WRITE (LUN, 5001)  (-2)
      SKIP = .FALSE.
      CALL PCHFD (4, X, F, D, 0, SKIP, 0, DUM, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -2, LUN, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .EQ. 3)  WRITE (LUN, 5001)  (-3)
      SKIP = .FALSE.
      CALL PCHFD (N, X, F, D, 1, SKIP, 0, DUM, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -3, LUN, KPRINT) )  NERR = NERR + 1
C
      IF (KPRINT .EQ. 3)  WRITE (LUN, 5001)  (-4)
      SKIP = .TRUE.
      CALL PCHFD (N, X, F, D, 1, SKIP, 0, DUM, DUM, DUM, IERR)
      IF (.NOT. COMP (IERR, -4, LUN, KPRINT) )  NERR = NERR + 1
C
C  SUMMARIZE RESULTS.
C
      IF (NERR .EQ. 0)  THEN
         FAIL = .FALSE.
         IF (KPRINT .GE. 2)  WRITE (LUN, 5002)
      ELSE
         FAIL = .TRUE.
         IF (KPRINT .GE. 2)  WRITE (LUN, 5003)  NERR
      ENDIF
C
C  TERMINATE.
C
      RETURN
C
C  FORMATS.
C
 5000 FORMAT ('1'//10X,'EVERCK RESULTS'/10X,'--------------')
 5001 FORMAT (/' THIS CALL SHOULD RETURN IERR =',I3)
 5002 FORMAT (///' ALL ERROR RETURNS OK.')
 5003 FORMAT (///' ***** TROUBLE IN EVERCK *****'
     *         //5X,I5,' TESTS FAILED TO GIVE EXPECTED RESULTS.')
C
      END
      LOGICAL FUNCTION  COMP (IERACT, IEREXP, LUN, KPRINT)
      INTEGER  IERACT, IEREXP, LUN, KPRINT
C
C     COMPARE ACTUAL VALUE OF IERR WITH EXPECTED VALUE.
C        PRINT ERROR MESSAGE IF THEY DON'T AGREE.
C
      IF (IERACT .EQ. IEREXP)  THEN
         COMP = .TRUE.
         IF (KPRINT .EQ. 3)  WRITE (LUN, 5010)
 5010    FORMAT (' OK.')
      ELSE
         COMP = .FALSE.
         IF (KPRINT .EQ. 3)  WRITE (LUN, 5020)  IERACT
 5020    FORMAT (' *** COMPARE FAILED -- IERR =',I5)
      ENDIF
C
      RETURN
      END
      SUBROUTINE PCHQK2 (KPRINT, IPASS)
      COMMON/UNIT/LUN
      INTEGER  LUN, KPRINT, IPASS
C
C-----------------------------------------------------------------------
C
C              PCHIP QUICKCHECK NUMBER 2
C
C     TESTS THE INTEGRATORS:  PCHIA, PCHID.
C
C-----------------------------------------------------------------------
C
C     RETURNED VALUES:
C        IPASS = 0  IF ALL TESTS PASSED.
C        IPASS = K > 0  IF K TEST INTEGRALS FAILED.
C
C-----------------------------------------------------------------------
C
C  DECLARE VARIABLES.
C
      INTEGER  I, IEREXP(17), IERR, N, NPAIRS
      REAL  A(17), B(17), CALC, D(7), ERRMAX, ERROR, F(7), MACHEP,
     *      ONE, THREE, THRQTR, TOL, TRUE, TWO, X(7)
      LOGICAL  FAIL, SKIP
      REAL  PCHIA
C
C  INITIALIZE.
C
      DATA  THRQTR /0.75/,  ONE /1./,  TWO /2./,  THREE /3./
C
      DATA  N /7/
      DATA  X /-4., -2., -0.9, 0., 0.9, 2., 4./
C
      DATA  NPAIRS /17/
      DATA  A /-3.0, 3.0,-0.5,-0.5,-0.5,-4.0,-4.0, 3.0,-5.0,-5.0,-6.0,
     *          6.0,-1.5,-1.5,-3.0, 3.0, 0.5/
      DATA  B / 3.0,-3.0, 1.0, 2.0, 5.0,-0.5, 4.0, 5.0,-3.0, 5.0,-5.0,
     *          5.0,-0.5,-1.0,-2.5, 3.5, 0.5/
      DATA  IEREXP /0,0,0,0,2,0,0,2,1,3,3,3,0,0,0,0,0/
C
C  DEFINE TEST FUNCTIONS.
C
      FCN(XDUM) = THREE*XDUM*XDUM*(XDUM-TWO)
      DERIV(XDUM) = THREE*XDUM*(TWO*(XDUM-TWO) + XDUM)
      ANTDER(XDUM) = XDUM**3 * (THRQTR*XDUM - TWO)
C
C  SET PASS/FAIL TOLERANCE.
C
      MACHEP = R1MACH(4)
      TOL = 100.*MACHEP
C
C  SET UP PCH FUNCTION DEFINITION.
C
      DO 10  I = 1, N
         F(I) =   FCN(X(I))
         D(I) = DERIV(X(I))
   10 CONTINUE
C
      IF (KPRINT .EQ. 3)  WRITE (LUN, 1000)  (X(I), F(I), D(I), I=1,N)
      IF (KPRINT .EQ. 2)  WRITE (LUN, 1001)
C
C  LOOP OVER (A,B)-PAIRS.
C
      IF (KPRINT .EQ. 3)  WRITE (LUN, 2000)
C
      IPASS = 0
C
      SKIP = .FALSE.
      DO 20  I = 1, NPAIRS
C               ---------------------------------------------
         CALC = PCHIA (N, X, F, D, 1, SKIP, A(I), B(I), IERR)
C               ---------------------------------------------
         IF (IERR .GE. 0)  THEN
            FAIL = IERR .NE. IEREXP(I)
            TRUE = ANTDER(B(I)) - ANTDER(A(I))
            ERROR = CALC - TRUE
            IF (KPRINT .EQ. 3)  THEN
               IF (FAIL)  THEN
                 WRITE (LUN, 2001) A(I), B(I), IERR, TRUE, CALC, ERROR,
     *                                          IEREXP(I)
               ELSE
                 WRITE (LUN, 2002) A(I), B(I), IERR, TRUE, CALC, ERROR
               ENDIF
            ENDIF
C
            ERROR = ABS(ERROR) / AMAX1(ONE, ABS(TRUE))
            IF (FAIL .OR. (ERROR.GT.TOL))  IPASS = IPASS + 1
            IF (I .EQ. 1)  THEN
               ERRMAX = ERROR
            ELSE
               ERRMAX = AMAX1(ERRMAX, ERROR)
            ENDIF
         ELSE
            IF (KPRINT .EQ. 3)  WRITE (LUN, 2002)  A(I), B(I), IERR
            IPASS = IPASS + 1
         ENDIF
   20 CONTINUE
C
      IF (KPRINT .GE. 2)  THEN
         WRITE (LUN, 2003)  ERRMAX, TOL
         IF (IPASS .EQ. 0)  THEN
            WRITE (LUN, 3000)
         ELSE
            WRITE (LUN ,3001)  IPASS
         ENDIF
      ENDIF
C
C  TERMINATE.
C
      IF(KPRINT.GE.1.AND.IPASS.NE.0) WRITE(LUN,99999)
      IF(KPRINT.GE.2.AND.IPASS.EQ.0) WRITE(LUN,99998)
99998 FORMAT(/'   ************PCHIP PASSED  ALL TESTS ****************')
99999 FORMAT(/'   ************PCHIP FAILED SOME TESTS ****************')
      RETURN
C
C  FORMATS.
C
 1000 FORMAT ('1'//10X,'CODE TO TEST PCHIP INTEGRATORS'
     *           // 5X,'DATA:' //11X,'X',9X,'F',9X,'D'
     *            /(5X,3F10.3) )
 1001 FORMAT ('1'//10X,'PCHQK2 RESULTS'/10X,'--------------')
 2000 FORMAT (// 5X,'TEST RESULTS:'
     *        //'    A     B    ERR     TRUE',16X,'CALC',15X,'ERROR')
 2001 FORMAT (/2F6.1,I5,2E20.10,E15.5,'  (',I1,') *****' )
 2002 FORMAT (/2F6.1,I5,2E20.10,E15.5)
 2003 FORMAT (/'  MAXIMUM RELATIVE ERROR IS:',E15.5,
     *                        5X,'TOLERANCE:',E15.5)
 3000 FORMAT (/' ALL PCHIP INTEGRATOR TESTS PASSED OK.')
 3001 FORMAT (/' *** TROUBLE *** ',I5,' PCHIP INTEGRATOR TESTS FAILED.')
C
      END
