      INTEGER FUNCTION DCHFMC(D1,D2,DELTA)
C***BEGIN PROLOGUE  DCHFMC
C***REFER TO  DPCHMC
C***ROUTINES CALLED  D1MACH
C***REVISION DATE  870707   (YYMMDD)
C***DESCRIPTION
C
C          DCHFMC:  Cubic Hermite Function Monotonicity Checker.
C
C    Called by  DPCHMC  to determine the monotonicity properties of the
C    cubic with boundary derivative values D1,D2 and chord slope DELTA.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        DOUBLE PRECISION  D1, D2, DELTA
C        INTEGER  ISMON, DCHFMC
C
C        ISMON = DCHFMC (D1, D2, DELTA)
C
C   Parameters:
C
C     D1,D2 -- (input) derivative values at the ends of an interval.
C
C     DELTA -- (input) data slope over that interval.
C
C     ISMON -- (output) integer function value, indicating the monoto-
C           nicity of the cubic segment:
C             ISMON = -1  if function is strictly decreasing;
C             ISMON =  0  if function is constant;
C             ISMON =  1  if function is strictly increasing;
C             ISMON =  2  if function is non-monotonic;
C             ISMON =  3  if unable to determine.
C
C  Fortran intrinsics used:  DSIGN.
C  Other routines used:  D1MACH.
C
C***END PROLOGUE  DCHFMC
C
C ----------------------------------------------------------------------
C
C  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
C                  Mathematics and Statistics Division,
C                  Lawrence Livermore National Laboratory.
C
C  Change record:
C     82-08-05   Converted to SLATEC library version.
C     83-12-01   Changed from  ISIGN  to SIGN  to correct bug that
C                produced wrong sign when -1 .LT. DELTA .LT. 0 .
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DCHFMC to CHFMC wherever it occurs,
C        b. Change the double precision declarations to real, and
C        c. Change the constants ZERO, ONE, ... to single precision.
C
C  DECLARE ARGUMENTS.
C
      DOUBLE PRECISION  D1, D2, DELTA, D1MACH
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER ISMON, ITRUE
      DOUBLE PRECISION  A, B, EPS, FOUR, ONE, PHI, TEN, THREE, TWO,
     * ZERO
C
C  INITIALIZE.
C
      DATA ZERO /0.D0/, ONE/1.D0/, TWO/2.D0/, THREE/3.D0/, FOUR/4.D0/,
     1      TEN /10.D0/
C
C        MACHINE-DEPENDENT PARAMETER -- SHOULD BE ABOUT 10*UROUND.
C***FIRST EXECUTABLE STATEMENT  DCHFMC
      EPS = TEN*D1MACH(4)
C
C  MAKE THE CHECK.
C
      IF (DELTA .EQ. ZERO)  THEN
C        CASE OF CONSTANT DATA.
         IF ((D1.EQ.ZERO) .AND. (D2.EQ.ZERO))  THEN
            ISMON = 0
         ELSE
            ISMON = 2
         ENDIF
      ELSE
C        DATA IS NOT CONSTANT -- PICK UP SIGN.
         ITRUE = DSIGN (ONE, DELTA)
         A = D1/DELTA
         B = D2/DELTA
         IF ((A.LT.ZERO) .OR. (B.LT.ZERO))  THEN
            ISMON = 2
         ELSE IF ((A.LE.THREE-EPS) .AND. (B.LE.THREE-EPS))  THEN
C           INSIDE SQUARE (0,3)X(0,3)  IMPLIES   OK.
            ISMON = ITRUE
         ELSE IF ((A.GT.FOUR+EPS) .AND. (B.GT.FOUR+EPS))  THEN
C           OUTSIDE SQUARE (0,4)X(0,4)  IMPLIES   NONMONOTONIC.
            ISMON = 2
         ELSE
C           MUST CHECK AGAINST BOUNDARY OF ELLIPSE.
            A = A - TWO
            B = B - TWO
            PHI = ((A*A + B*B) + A*B) - THREE
            IF (PHI .LT. -EPS)  THEN
               ISMON = ITRUE
            ELSE IF (PHI .GT. EPS)  THEN
               ISMON = 2
            ELSE
C              TO CLOSE TO BOUNDARY TO TELL,
C                  IN THE PRESENCE OF ROUND-OFF ERRORS.
               ISMON = 3
            ENDIF
         ENDIF
      ENDIF
C
C  RETURN VALUE.
C
      DCHFMC = ISMON
      RETURN
C------------- LAST LINE OF DCHFMC FOLLOWS -----------------------------
      END
