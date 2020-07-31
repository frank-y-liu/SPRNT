      REAL FUNCTION RELDST(P, D, X, X0)
C
C  ***  COMPUTE AND RETURN RELATIVE DIFFERENCE BETWEEN X AND X0  ***
C  ***  NL2SOL VERSION 2.2  ***
C
      INTEGER P
      REAL D(P), X(P), X0(P)
C/+
      REAL  ABS
C/
      INTEGER I
      REAL EMAX, T, XMAX, ZERO
C/6
      DATA ZERO/0.E+0/
C/7
C     PARAMETER (ZERO=0.E+0)
C/
C
      EMAX = ZERO
      XMAX = ZERO
      DO 10 I = 1, P
         T =  ABS(D(I) * (X(I) - X0(I)))
         IF (EMAX .LT. T) EMAX = T
         T = D(I) * ( ABS(X(I)) +  ABS(X0(I)))
         IF (XMAX .LT. T) XMAX = T
 10      CONTINUE
      RELDST = ZERO
      IF (XMAX .GT. ZERO) RELDST = EMAX / XMAX
 999  RETURN
C  ***  LAST CARD OF RELDST FOLLOWS  ***
      END
