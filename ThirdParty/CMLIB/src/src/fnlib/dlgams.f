      SUBROUTINE DLGAMS(X,DLGAM,SGNGAM)
C***BEGIN PROLOGUE  DLGAMS
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7A
C***KEYWORDS  ABSOLUTE VALUE,DOUBLE PRECISION,GAMMA FUNCTION,LOGARITHM,
C             SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Calculates the log of the absolute value of the Gamma
C            function
C***DESCRIPTION
C
C DLGAMS(X,DLGAM,SGNGAM) calculates the double precision natural
C logarithm of the absolute value of the gamma function for
C double precision argument X and stores the result in double
C precision argument DLGAM.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DINT,DLNGAM
C***END PROLOGUE  DLGAMS
      DOUBLE PRECISION X, DLGAM, SGNGAM, DINT, DLNGAM
C***FIRST EXECUTABLE STATEMENT  DLGAMS
      DLGAM = DLNGAM(X)
      SGNGAM = 1.0D0
      IF (X.GT.0.D0) RETURN
C
      INT = DMOD (-DINT(X), 2.0D0) + 0.1D0
      IF (INT.EQ.0) SGNGAM = -1.0D0
C
      RETURN
      END
