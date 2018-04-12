      FUNCTION GAMI(A,X)
C***BEGIN PROLOGUE  GAMI
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  GAMMA FUNCTION,INCOMPLETE GAMMA FUNCTION,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the incomplete Gamma function.
C***DESCRIPTION
C
C Evaluate the incomplete gamma function defined by
C
C GAMI = integral from T = 0 to X of EXP(-T) * T**(A-1.0) .
C
C GAMI is evaluated for positive values of A and non-negative values
C of X.  A slight deterioration of 2 or 3 digits accuracy will occur
C when GAMI is very large or very small, because logarithmic variables
C are used.  GAMI, A, and X are single precision.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALNGAM,GAMIT,XERROR
C***END PROLOGUE  GAMI
C***FIRST EXECUTABLE STATEMENT  GAMI
      IF (A.LE.0.0) CALL XERROR ( 'GAMI    A MUST BE GT ZERO', 25, 1, 2)
      IF (X.LT.0.0) CALL XERROR ( 'GAMI    X MUST BE GE ZERO', 25, 2, 2)
C
      GAMI = 0.0
      IF (X.EQ.0.0) RETURN
C
C THE ONLY ERROR POSSIBLE IN THE EXPRESSION BELOW IS A FATAL OVERFLOW.
      FACTOR = EXP (ALNGAM(A) + A*ALOG(X) )
C
      GAMI = FACTOR * GAMIT(A, X)
C
      RETURN
      END
