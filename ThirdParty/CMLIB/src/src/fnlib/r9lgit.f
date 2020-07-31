      FUNCTION R9LGIT(A,X,ALGAP1)
C***BEGIN PROLOGUE  R9LGIT
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  INCOMPLETE GAMMA FUNCTION,LOGARITHM,
C             PERRON-S CONTINUED FRACTION,SPECIAL FUNCTION,TRICOMI-S
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the log of Tricomi's incomplete Gamma function
C            with Perron's continued fraction for large X and A .GE. X.
C***DESCRIPTION
C
C Compute the log of Tricomi's incomplete gamma function with Perron's
C continued fraction for large X and for A .GE. X.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH,XERROR
C***END PROLOGUE  R9LGIT
      DATA EPS, SQEPS / 2*0.0 /
C***FIRST EXECUTABLE STATEMENT  R9LGIT
      IF (EPS.EQ.0.0) EPS = 0.5*R1MACH(3)
      IF (SQEPS.EQ.0.0) SQEPS = SQRT(R1MACH(4))
C
      IF (X.LE.0.0 .OR. A.LT.X) CALL XERROR ( 'R9LGIT  X SHOULD BE GT 0.
     10 AND LE A', 35, 2, 2)
C
      AX = A + X
      A1X = AX + 1.0
      R = 0.0
      P = 1.0
      S = P
      DO 20 K=1,200
        FK = K
        T = (A+FK)*X*(1.0+R)
        R = T/((AX+FK)*(A1X+FK)-T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 30
 20   CONTINUE
      CALL XERROR ( 'R9LGIT  NO CONVERGENCE IN 200 TERMS OF CONTINUED FR
     1ACTION', 57, 3, 2)
C
 30   HSTAR = 1.0 - X*S/A1X
      IF (HSTAR.LT.SQEPS) CALL XERROR ( 'R9LGIT  RESULT LESS THAN HALF P
     1RECISION', 39, 1, 1)
C
      R9LGIT = -X - ALGAP1 - ALOG(HSTAR)
C
      RETURN
      END
