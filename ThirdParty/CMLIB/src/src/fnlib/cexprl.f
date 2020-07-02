      COMPLEX FUNCTION CEXPRL(Z)
C***BEGIN PROLOGUE  CEXPRL
C***DATE WRITTEN   770801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY (YYMMDD)
C   000601  Changed CEXP to generic EXP
C***CATEGORY NO.  C4B
C***KEYWORDS  COMPLEX,ELEMENTARY FUNCTION,EXPONENTIAL,FIRST ORDER,
C             RELATIVE ERROR
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Evaluates (CEXP(Z)-1)/Z so that EXP(Z) = 1 + Z*CEXPRL(Z)
C***DESCRIPTION
C
C Evaluate  (CEXP(Z)-1)/Z .  For small CABS(Z), we use the Taylor
C series.  We could instead use the expression
C        CEXPRL(Z) = (EXP(X)*EXP(I*Y)-1)/Z
C                  = (X*EXPREL(X) * (1 - 2*SIN(Y/2)**2) - 2*SIN(Y/2)**2
C                                    + I*SIN(Y)*(1+X*EXPREL(X))) / Z
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH
C***END PROLOGUE  CEXPRL
      COMPLEX Z
      DATA NTERMS, RBND / 0, 0.0 /
C***FIRST EXECUTABLE STATEMENT  CEXPRL
      IF (NTERMS.NE.0) GO TO 10
      ALNEPS = ALOG(R1MACH(3))
      XN = 3.72 - 0.3*ALNEPS
      XLN = ALOG((XN+1.0)/1.36)
      NTERMS = XN - (XN*XLN+ALNEPS)/(XLN+1.36) + 1.5
      RBND = R1MACH(3)
C
 10   R = CABS(Z)
      IF (R.GT.0.5) CEXPRL = (EXP(Z) - 1.0) / Z
      IF (R.GT.0.5) RETURN
C
      CEXPRL = (1.0, 0.0)
      IF (R.LT.RBND) RETURN
C
      CEXPRL = (0.0, 0.0)
      DO 20 I=1,NTERMS
        CEXPRL = 1.0 + CEXPRL*Z/FLOAT(NTERMS+2-I)
 20   CONTINUE
C
      RETURN
      END
