      DOUBLE PRECISION FUNCTION DGAMI(A,X)
C***BEGIN PROLOGUE  DGAMI
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  DOUBLE PRECISION,GAMMA,GAMMA FUNCTION,
C             INCOMPLETE GAMMA FUNCTION,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Calculates the d.p. incomplete Gamma function.
C***DESCRIPTION
C
C Evaluate the incomplete gamma function defined by
C
C DGAMI = integral from T = 0 to X of EXP(-T) * T**(A-1.0) .
C
C DGAMI is evaluated for positive values of A and non-negative values
C of X.  A slight deterioration of 2 or 3 digits accuracy will occur
C when DGAMI is very large or very small, because logarithmic variables
C are used.  The function and both arguments are double precision.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DGAMIT,DLNGAM,XERROR
C***END PROLOGUE  DGAMI
      DOUBLE PRECISION A, X, FACTOR, DLNGAM, DGAMIT
C***FIRST EXECUTABLE STATEMENT  DGAMI
      IF (A.LE.0.D0) CALL XERROR ( 'DGAMI   A MUST BE GT ZERO', 25, 1,2)
      IF (X.LT.0.D0) CALL XERROR ( 'DGAMI   X MUST BE GE ZERO', 25, 2,2)
C
      DGAMI = 0.D0
      IF (X.EQ.0.0D0) RETURN
C
C THE ONLY ERROR POSSIBLE IN THE EXPRESSION BELOW IS A FATAL OVERFLOW.
      FACTOR = DEXP (DLNGAM(A) + A*DLOG(X))
C
      DGAMI = FACTOR * DGAMIT (A, X)
C
      RETURN
      END
