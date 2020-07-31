      DOUBLE PRECISION FUNCTION D9LGIT(A,X,ALGAP1)
C***BEGIN PROLOGUE  D9LGIT
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  DOUBLE PRECISION,GAMMA,INCOMPLETE GAMMA FUNCTION,
C             LOGARITHM,SPECIAL FUNCTION,TRICOMI
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes  the log of Tricomi's incomplete Gamma function
C            with Perron's continued fraction for large X and A .GE. X.
C***DESCRIPTION
C
C Compute the log of Tricomi's incomplete gamma function with Perron's
C continued fraction for large X and for A .GE. X.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,XERROR
C***END PROLOGUE  D9LGIT
      DOUBLE PRECISION A, X, ALGAP1, AX, A1X, EPS, FK, HSTAR, P, R, S,
     1  SQEPS, T, D1MACH
      DATA EPS, SQEPS / 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  D9LGIT
      IF (EPS.NE.0.D0) GO TO 10
      EPS = 0.5D0*D1MACH(3)
      SQEPS = DSQRT (D1MACH(4))
C
 10   IF (X.LE.0.D0 .OR. A.LT.X) CALL XERROR ( 'D9LGIT  X SHOULD BE GT 0
     1.0 AND LE A', 35, 2, 2)
C
      AX = A + X
      A1X = AX + 1.0D0
      R = 0.D0
      P = 1.D0
      S = P
      DO 20 K=1,200
        FK = K
        T = (A+FK)*X*(1.D0+R)
        R = T/((AX+FK)*(A1X+FK)-T)
        P = R*P
        S = S + P
        IF (DABS(P).LT.EPS*S) GO TO 30
 20   CONTINUE
      CALL XERROR ( 'D9LGIT  NO CONVERGENCE IN 200 TERMS OF CONTINUED FR
     1ACTION', 57, 3, 2)
C
 30   HSTAR = 1.0D0 - X*S/A1X
      IF (HSTAR.LT.SQEPS) CALL XERROR ( 'D9LGIT  RESULT LESS THAN HALF P
     1RECISION', 39, 1, 1)
C
      D9LGIT = -X - ALGAP1 - DLOG(HSTAR)
      RETURN
C
      END
