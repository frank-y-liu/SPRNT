      DOUBLE PRECISION FUNCTION DEXPRL(X)
C***BEGIN PROLOGUE  DEXPRL
C***DATE WRITTEN   770801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C4B
C***KEYWORDS  DOUBLE PRECISION,ELEMENTARY FUNCTION,EXPONENTIAL,
C             FIRST ORDER,RELATIVE ACCURACY
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Calculates the d.p. relative error exponential
C            (DEXP(X)-1)/X.
C***DESCRIPTION
C
C Evaluate  EXPREL(X) = (EXP(X) - 1.0) / X.   For small ABS(X) the
C Taylor series is used.  If X is negative the reflection formula
C         EXPREL(X) = EXP(X) * EXPREL(ABS(X))
C may be used.  This reflection formula will be of use when the
C evaluation for small ABS(X) is done by Chebyshev series rather than
C Taylor series.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH
C***END PROLOGUE  DEXPRL
      DOUBLE PRECISION X, ABSX, ALNEPS, XBND, XLN, XN,  D1MACH
      DATA NTERMS, XBND / 0, 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  DEXPRL
      IF (NTERMS.NE.0) GO TO 10
      ALNEPS = DLOG(D1MACH(3))
      XN = 3.72D0 - 0.3D0*ALNEPS
      XLN = DLOG((XN+1.0D0)/1.36D0)
      NTERMS = XN - (XN*XLN+ALNEPS)/(XLN+1.36D0) + 1.5D0
      XBND = D1MACH(3)
C
 10   ABSX = DABS(X)
      IF (ABSX.GT.0.5D0) DEXPRL = (DEXP(X)-1.0D0)/X
      IF (ABSX.GT.0.5D0) RETURN
C
      DEXPRL = 1.0D0
      IF (ABSX.LT.XBND) RETURN
C
      DEXPRL = 0.0D0
      DO 20 I=1,NTERMS
        DEXPRL = 1.0D0 + DEXPRL*X/DBLE(FLOAT(NTERMS+2-I))
 20   CONTINUE
C
      RETURN
      END
