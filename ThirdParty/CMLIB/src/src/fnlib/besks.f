      SUBROUTINE BESKS(XNU,X,NIN,BK)
C***BEGIN PROLOGUE  BESKS
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  C10B3
C***KEYWORDS  BESSEL FUNCTION,FRACTIONAL ORDER,MODIFIED BESSEL FUNCTION,
C             SEQUENCE,THIRD KIND
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes a sequence of modified Bessel functions of the
C            third kind of fractional order.
C***DESCRIPTION
C
C BESKS computes a sequence of modified Bessel functions of the third
C kind of order XNU + I at X, where X .GT. 0, XNU lies in (-1,1),
C and I = 0, 1, ... , NIN - 1, if NIN is positive and I = 0, 1, ... ,
C NIN + 1, if NIN is negative.  On return, the vector BK(.) Contains
C the results at X for order starting at XNU.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  BESKES,R1MACH,XERROR
C***END PROLOGUE  BESKS
      DIMENSION BK(*)
      DATA XMAX / 0.0 /
C***FIRST EXECUTABLE STATEMENT  BESKS
      IF (XMAX.EQ.0.0) XMAX = -ALOG (R1MACH(1))
C
      IF (X.GT.XMAX) CALL XERROR ( 'BESKS   X SO BIG BESSEL K UNDERFLOWS
     1', 36, 1, 2)
C
      CALL BESKES (XNU, X, NIN, BK)
C
      EXPXI = EXP (-X)
      N = IABS (NIN)
      DO 20 I=1,N
        BK(I) = EXPXI * BK(I)
 20   CONTINUE
C
      RETURN
      END
