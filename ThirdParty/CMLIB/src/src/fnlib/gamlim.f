      SUBROUTINE GAMLIM(XMIN,XMAX)
C***BEGIN PROLOGUE  GAMLIM
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7A,R2
C***KEYWORDS  GAMMA FUNCTION,LIMITS,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the minimum and maximum bounds for X in GAMMA(X).
C***DESCRIPTION
C
C Calculate the minimum and maximum legal bounds for X in GAMMA(X).
C XMIN and XMAX are not the only bounds, but they are the only non-
C trivial ones to calculate.
C
C             Output Arguments --
C XMIN   minimum legal value of X in GAMMA(X).  Any smaller value of
C        X might result in underflow.
C XMAX   maximum legal value of X in GAMMA(X).  Any larger value will
C        cause overflow.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH,XERROR
C***END PROLOGUE  GAMLIM
C***FIRST EXECUTABLE STATEMENT  GAMLIM
      ALNSML = ALOG(R1MACH(1))
      XMIN = -ALNSML
      DO 10 I=1,10
        XOLD = XMIN
        XLN = ALOG(XMIN)
        XMIN = XMIN - XMIN*((XMIN+0.5)*XLN - XMIN - 0.2258 + ALNSML)
     1    / (XMIN*XLN + 0.5)
        IF (ABS(XMIN-XOLD).LT.0.005) GO TO 20
 10   CONTINUE
      CALL XERROR ( 'GAMLIM  UNABLE TO FIND XMIN', 27, 1, 2)
C
 20   XMIN = -XMIN + 0.01
C
      ALNBIG = ALOG(R1MACH(2))
      XMAX = ALNBIG
      DO 30 I=1,10
        XOLD = XMAX
        XLN = ALOG(XMAX)
        XMAX = XMAX - XMAX*((XMAX-0.5)*XLN - XMAX + 0.9189 - ALNBIG)
     1    / (XMAX*XLN - 0.5)
        IF (ABS(XMAX-XOLD).LT.0.005) GO TO 40
 30   CONTINUE
      CALL XERROR ( 'GAMLIM  UNABLE TO FIND XMAX', 27, 2, 2)
C
 40   XMAX = XMAX - 0.01
      XMIN = AMAX1 (XMIN, -XMAX+1.)
C
      RETURN
      END
