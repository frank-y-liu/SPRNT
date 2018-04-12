      DOUBLE PRECISION FUNCTION DCHFIV(X1,X2,F1,F2,D1,D2,A,B,IERR)
C***BEGIN PROLOGUE  DCHFIV
C***REFER TO  DPCHIA
C***ROUTINES CALLED  XERROR
C***REVISION DATE  870707   (YYMMDD)
C***DESCRIPTION
C
C          DCHFIV:  Cubic Hermite Function Integral Evaluator.
C
C     Called by  DPCHIA  to evaluate the integral of a single cubic (in
C     Hermite form) over an arbitrary interval (A,B).
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        INTEGER  IERR
C        DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, A, B
C        DOUBLE PRECISION  VALUE, DCHFIV
C
C        VALUE = DCHFIV (X1, X2, F1, F2, D1, D2, A, B, IERR)
C
C   Parameters:
C
C     VALUE -- (output) value of the requested integral.
C
C     X1,X2 -- (input) endpoints if interval of definition of cubic.
C           (Must be distinct.  Error return if not.)
C
C     F1,F2 -- (input) function values at the ends of the interval.
C
C     D1,D2 -- (input) derivative values at the ends of the interval.
C
C     A,B -- (input) endpoints of interval of integration.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0 (no errors).
C           "Recoverable errors":
C              IERR = -1  if X1.EQ.X2 .
C                (VALUE has not been set in this case.)
C
C***END PROLOGUE  DCHFIV
C
C ----------------------------------------------------------------------
C
C  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
C                  Mathematics and Statistics Division,
C                  Lawrence Livermore National Laboratory.
C
C  Change record:
C     82-08-05   Converted to SLATEC library version.
C     87-07-07   Corrected XERROR calls for d.p. name(s).
C                The name was changed from DCHIV to DCHFIV.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DCHFIV to CHFIV wherever it occurs,
C        b. Change the double precision declarations to real, and
C        c. Change the constants HALF, TWO, ... to single precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER IERR
      DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, A, B
C
C  DECLARE LOCAL VARIABLES.
C
      DOUBLE PRECISION  DTERM, FOUR, FTERM, H, HALF, PHIA1, PHIA2,
     *      PHIB1, PHIB2, PSIA1, PSIA2, PSIB1, PSIB2, TA1, TA2, TB1,
     *      TB2, THREE, TWO, UA1, UA2, UB1, UB2
C
C  INITIALIZE.
C
      DATA  HALF/.5D0/, TWO/2.D0/, THREE/3.D0/, FOUR/4.D0/, SIX/6.D0/
C
C  VALIDITY CHECK INPUT.
C
C***FIRST EXECUTABLE STATEMENT  DCHFIV
      IF (X1 .EQ. X2)  GO TO 5001
      IERR = 0
C
C  COMPUTE INTEGRAL.
C
      H = X2 - X1
      TA1 = (A - X1) / H
      TA2 = (X2 - A) / H
      TB1 = (B - X1) / H
      TB2 = (X2 - B) / H
C
      UA1 = TA1**3
      PHIA1 = UA1 * (TWO - TA1)
      PSIA1 = UA1 * (THREE*TA1 - FOUR)
      UA2 = TA2**3
      PHIA2 =  UA2 * (TWO - TA2)
      PSIA2 = -UA2 * (THREE*TA2 - FOUR)
C
      UB1 = TB1**3
      PHIB1 = UB1 * (TWO - TB1)
      PSIB1 = UB1 * (THREE*TB1 - FOUR)
      UB2 = TB2**3
      PHIB2 =  UB2 * (TWO - TB2)
      PSIB2 = -UB2 * (THREE*TB2 - FOUR)
C
      FTERM =   F1*(PHIA2 - PHIB2) + F2*(PHIB1 - PHIA1)
      DTERM = ( D1*(PSIA2 - PSIB2) + D2*(PSIB1 - PSIA1) )*(H/SIX)
C
C  RETURN VALUE.
C
      DCHFIV = (HALF*H) * (FTERM + DTERM)
      RETURN
C
C  ERROR RETURN.
C
 5001 CONTINUE
      IERR = -1
      CALL XERROR ('DCHFIV -- X1 EQUAL TO X2'
     *           , 0, IERR, 1)
      RETURN
C------------- LAST LINE OF DCHFIV FOLLOWS -----------------------------
      END
