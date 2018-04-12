      INTEGER FUNCTION CHFMC(D1,D2,DELTA)
C***BEGIN PROLOGUE  CHFMC
C***REFER TO  PCHMC
C***ROUTINES CALLED  R1MACH
C***REVISION DATE  870707   (YYMMDD)
C***DESCRIPTION
C
C          CHFMC:  Cubic Hermite Function Monotonicity Checker.
C
C    Called by  PCHMC  to determine the monotonicity properties of the
C    cubic with boundary derivative values D1,D2 and chord slope DELTA.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        REAL  D1, D2, DELTA
C        INTEGER  ISMON, CHFMC
C
C        ISMON = CHFMC (D1, D2, DELTA)
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
C  Fortran intrinsics used:  SIGN.
C  Other routines used:  R1MACH.
C
C***END PROLOGUE  CHFMC
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
C     To produce a double precision version, simply:
C        a. Change CHFMC to DCHFMC wherever it occurs,
C        b. Change the real declarations to double precision, and
C        c. Change the constants ZERO, ONE, ... to double precision.
C
C  DECLARE ARGUMENTS.
C
      REAL  D1, D2, DELTA
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  ISMON, ITRUE
      REAL  A, B, EPS, FOUR, ONE, PHI, TEN, THREE, TWO, ZERO
C
C  INITIALIZE.
C
      DATA  ZERO /0./,  ONE /1.0/,  TWO /2./,  THREE /3./,  FOUR /4./,
     1      TEN /10./
C
C        MACHINE-DEPENDENT PARAMETER -- SHOULD BE ABOUT 10*UROUND.
C***FIRST EXECUTABLE STATEMENT  CHFMC
      EPS = TEN*R1MACH(4)
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
         ITRUE = SIGN (ONE, DELTA)
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
      CHFMC = ISMON
      RETURN
C------------- LAST LINE OF CHFMC FOLLOWS ------------------------------
      END
