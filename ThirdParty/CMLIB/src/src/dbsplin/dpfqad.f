      SUBROUTINE DPFQAD(F,LDC,C,XI,LXI,K,ID,X1,X2,TOL,QUAD,IERR)
C***BEGIN PROLOGUE  DPFQAD
C***DATE WRITTEN   800901   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  H2A2A1,E3,K6
C***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
C             QUADRATURE,SPLINE
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  Computes the integral on (X1,X2) of a product of a
C            function F and the ID-th derivative of a B-spline,
C            (PP-representation).
C***DESCRIPTION
C
C     Written by D. E. Amos, June, 1979.
C
C     Reference SAND-79-1825
C
C     Abstract    **** a double precision routine ****
C
C         DPFQAD computes the integral on (X1,X2) of a product of a
C         function F and the ID-th derivative of a B-spline, using the
C         PP-representation (C,XI,LXI,K).  (X1,X2) is normally a sub-
C         interval of XI(1) .LE. X .LE. XI(LXI+1).  An integration
C         routine, DPPGQ8 (a modification of GAUS8), integrates the
C         product on subintervals of (X1,X2) formed by the included
C         break points.  Integration outside of (XI(1),XI(LXI+1)) is
C         permitted provided F is defined.
C
C         The maximum number of significant digits obtainable in
C         DBSQAD is the smaller of 18 and the number of digits
C         carried in double precision arithmetic.
C
C         DPFQAD calls DINTRV, DPPVAL, DPPGQ8, D1MACH, XERROR
C
C     Description of arguments
C         Input      F,C,XI,X1,X2,TOL are double precision
C           F      - external function of one argument for the
C                    integrand PF(X)=F(X)*DPPVAL(LDC,C,XI,LXI,K,ID,X,
C                    INPPV)
C           LDC    - leading dimension of matrix C, LDC .GE. K
C           C(I,J) - right Taylor derivatives at XI(J), I=1,K , J=1,LXI
C           XI(*)  - break point array of length LXI+1
C           LXI    - number of polynomial pieces
C           K      - order of B-spline, K .GE. 1
C           ID     - order of the spline derivative, 0 .LE. ID .LE. K-1
C                    ID=0 gives the spline function
C           X1,X2  - end points of quadrature interval, normally in
C                    XI(1) .LE. X .LE. XI(LXI+1)
C           TOL    - desired accuracy for the quadrature, suggest
C                    10.*DTOL .LT. TOL .LE. 0.1 where DTOL is the
C                    maximum of 1.0D-18 and double precision unit
C                    roundoff for the machine = D1MACH(4)
C
C         Output     QUAD is double precision
C           QUAD   - integral of PF(X) on (X1,X2)
C           IERR   - a status code
C                    IERR=1 normal return
C                         2 some quadrature does not meet the
C                           requested tolerance
C
C     Error Conditions
C         Improper input is a fatal error.
C         Some quadrature does not meet the requested tolerance.
C***REFERENCES  D.E. AMOS, *QUADRATURE SUBROUTINES FOR SPLINES AND
C                 B-SPLINES*, SAND79-1825, SANDIA LABORATORIES,
C                 DECEMBER 1979.
C***ROUTINES CALLED  D1MACH,DINTRV,DPPGQ8,XERROR
C***END PROLOGUE  DPFQAD
C
C
      INTEGER ID,IERR,IFLG,ILO,IL1,IL2,INPPV,K,LDC,LEFT,LXI,MF1,MF2
      DOUBLE PRECISION A,AA,ANS,B,BB,C,Q,QUAD,TA,TB,TOL,WTOL,XI,X1,X2
      DOUBLE PRECISION D1MACH, F
      DIMENSION XI(*), C(LDC,*)
      EXTERNAL F
C
C***FIRST EXECUTABLE STATEMENT  DPFQAD
      IERR = 1
      QUAD = 0.0D0
      IF(K.LT.1) GO TO 100
      IF(LDC.LT.K) GO TO 105
      IF(ID.LT.0 .OR. ID.GE.K) GO TO 110
      IF(LXI.LT.1) GO TO 115
      WTOL = D1MACH(4)
      WTOL = DMAX1(WTOL,1.0D-18)
      IF (TOL.LT.WTOL .OR. TOL.GT.0.1D0) GO TO 20
      AA = DMIN1(X1,X2)
      BB = DMAX1(X1,X2)
      IF (AA.EQ.BB) RETURN
      ILO = 1
      CALL DINTRV(XI, LXI, AA, ILO, IL1, MF1)
      CALL DINTRV(XI, LXI, BB, ILO, IL2, MF2)
      Q = 0.0D0
      INPPV = 1
      DO 10 LEFT=IL1,IL2
        TA = XI(LEFT)
        A = DMAX1(AA,TA)
        IF (LEFT.EQ.1) A = AA
        TB = BB
        IF (LEFT.LT.LXI) TB = XI(LEFT+1)
        B = DMIN1(BB,TB)
        CALL DPPGQ8(F,LDC,C,XI,LXI,K,ID,A,B,INPPV,TOL,ANS,IFLG)
        IF (IFLG.GT.1) IERR = 2
        Q = Q + ANS
   10 CONTINUE
      IF (X1.GT.X2) Q = -Q
      QUAD = Q
      RETURN
C
   20 CONTINUE
      CALL XERROR( ' DPFQAD,  TOL IS LESS DTOL OR GREATER THAN 0.1',
     1 46, 2, 1)
      RETURN
  100 CONTINUE
      CALL XERROR( ' DPFQAD,  K DOES NOT SATISFY K.GE.1', 35, 2, 1)
      RETURN
  105 CONTINUE
      CALL XERROR( ' DPFQAD,  LDC DOES NOT SATISFY LDC.GE.K', 39, 2, 1)
      RETURN
  110 CONTINUE
      CALL XERROR( ' DPFQAD,  ID DOES NOT SATISFY 0.LE.ID.LT.K', 42,
     1 2, 1)
      RETURN
  115 CONTINUE
      CALL XERROR( ' DPFQAD,  LXI DOES NOT SATISFY LXI.GE.1', 39, 2, 1)
      RETURN
      END
