      FUNCTION R9LGIC(A,X,ALX)
C***BEGIN PROLOGUE  R9LGIC
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUMNCTION,LARGE X,
C             LOGARITHM,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the log complementary incomplete Gamma function
C            for large X and for A .LE. X.
C***DESCRIPTION
C
C Compute the log complementary incomplete gamma function for large X
C and for A .LE. X.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH,XERROR
C***END PROLOGUE  R9LGIC
      DATA EPS / 0.0 /
C***FIRST EXECUTABLE STATEMENT  R9LGIC
      IF (EPS.EQ.0.0) EPS = 0.5*R1MACH(3)
C
      XPA = X + 1.0 - A
      XMA = X - 1.0 - A
C
      R = 0.0
      P = 1.0
      S = P
      DO 10 K=1,200
        FK = K
        T = FK*(A-FK)*(1.0+R)
        R = -T/((XMA+2.0*FK)*(XPA+2.0*FK)+T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 20
 10   CONTINUE
      CALL XERROR ( 'R9LGIC  NO CONVERGENCE IN 200 TERMS OF CONTINUED FR
     1ACTION', 57, 1, 2)
C
 20   R9LGIC = A*ALX - X + ALOG(S/XPA)
C
      RETURN
      END
