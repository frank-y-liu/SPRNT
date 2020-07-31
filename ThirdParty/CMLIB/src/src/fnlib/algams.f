      SUBROUTINE ALGAMS(X,ALGAM,SGNGAM)
C***BEGIN PROLOGUE  ALGAMS
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7A
C***KEYWORDS  GAMMA FUNCTION,LOGARITHM,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes log(ABS(GAMMA(X))).
C***DESCRIPTION
C
C Evaluates the logarithm of the absolute value of the gamma
C function.
C     X           - input argument
C     ALGAM       - result
C     SGNGAM      - is set to the sign of GAMMA(X) and will
C                   be returned at +1.0 or -1.0.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALNGAM
C***END PROLOGUE  ALGAMS
C***FIRST EXECUTABLE STATEMENT  ALGAMS
      ALGAM = ALNGAM(X)
      SGNGAM = 1.0
      IF (X.GT.0.0) RETURN
C
      INT = AMOD (-AINT(X), 2.0) + 0.1
      IF (INT.EQ.0) SGNGAM = -1.0
C
      RETURN
      END
