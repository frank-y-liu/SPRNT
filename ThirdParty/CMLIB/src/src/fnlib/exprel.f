      FUNCTION EXPREL(X)
C***BEGIN PROLOGUE  EXPREL
C***DATE WRITTEN   770801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C4B
C***KEYWORDS  ELEMENTARY FUNCTION,EXPONENTIAL,RELATIVE ERROR
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Evaluates EXPREL(X) = (EXP(X)-1)/X)
C***DESCRIPTION
C
C Evaluate  EXPREL(X) = (EXP(X) - 1.0) / X.   For small ABS(X) the
C Taylor series is used.  If X is negative, the reflection formula
C         EXPREL(X) = EXP(X) * EXPREL(ABS(X))
C may be used.  This reflection formula will be of use when the
C evaluation for small ABS(X) is done by Chebyshev series rather than
C Taylor series.  EXPREL and X are single precision.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH
C***END PROLOGUE  EXPREL
      DATA NTERMS, XBND / 0, 0.0 /
C***FIRST EXECUTABLE STATEMENT  EXPREL
      IF (NTERMS.NE.0) GO TO 10
      ALNEPS = ALOG(R1MACH(3))
      XN = 3.72 - 0.3*ALNEPS
      XLN = ALOG((XN+1.0)/1.36)
      NTERMS = XN - (XN*XLN+ALNEPS)/(XLN+1.36) + 1.5
      XBND = R1MACH(3)
C
 10   ABSX = ABS(X)
      IF (ABSX.GT.0.5) EXPREL = (EXP(X) - 1.0) / X
      IF (ABSX.GT.0.5) RETURN
C
      EXPREL = 1.0
      IF (ABSX.LT.XBND) RETURN
C
      EXPREL = 0.0
      DO 20 I=1,NTERMS
        EXPREL = 1.0 + EXPREL*X/FLOAT(NTERMS+2-I)
 20   CONTINUE
C
      RETURN
      END
