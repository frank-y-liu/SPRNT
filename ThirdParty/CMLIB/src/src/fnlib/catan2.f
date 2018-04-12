      COMPLEX FUNCTION CATAN2(CSN,CCS)
C***BEGIN PROLOGUE  CATAN2
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C4A
C***KEYWORDS  ARC TANGENT,COMPLEX,ELEMENTARY FUNCTION,QUADRANT
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the complex arc Tangent in the proper quadrant.
C***DESCRIPTION
C
C CATAN2(CSN,CCS) calculates the complex trigonometric arc
C tangent of the ratio CSN/CCS and returns a result whose real
C part is in the correct quadrant (within a multiple of 2*PI).  The
C result is in units of radians and the real part is between -PI
C and +PI.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CATAN,XERROR
C***END PROLOGUE  CATAN2
      COMPLEX CSN, CCS, CATAN
      DATA PI / 3.1415926535 8979323846E0 /
C***FIRST EXECUTABLE STATEMENT  CATAN2
      IF (CABS(CCS).EQ.0.) GO TO 10
C
      CATAN2 = CATAN (CSN/CCS)
      IF (REAL(CCS).LT.0.) CATAN2 = CATAN2 + PI
      IF (REAL(CATAN2).GT.PI) CATAN2 = CATAN2 - 2.0*PI
      RETURN
C
 10   IF (CABS(CSN).EQ.0.) CALL XERROR (  'CATAN2  CALLED WITH BOTH ARGS
     1 ZERO', 34, 1, 2)
C
      CATAN2 = CMPLX (SIGN(0.5*PI,REAL(CSN)), 0.0)
C
      RETURN
      END
