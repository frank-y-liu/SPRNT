      DOUBLE PRECISION FUNCTION DEI(X)
C***BEGIN PROLOGUE  DEI
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C5
C***KEYWORDS  DOUBLE PRECISION,EXPONENTIAL INTEGRAL,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. exponential integral EI(X).
C***DESCRIPTION
C
C DEI(X) calculates the double precision exponential integral Ei(X)
C for positive double precision argument X and the Cauchy principal
C value for negative X.  If principal values are used everywhere,
C then for all X
C
C        Ei(X) = -Ei(-X)   or
C        Ei(X) = -Ei(-X).
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DE1
C***END PROLOGUE  DEI
      DOUBLE PRECISION X, DE1
C***FIRST EXECUTABLE STATEMENT  DEI
      DEI = -DE1(-X)
C
      RETURN
      END
