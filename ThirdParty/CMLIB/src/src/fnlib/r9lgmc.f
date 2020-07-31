      FUNCTION R9LGMC(X)
C***BEGIN PROLOGUE  R9LGMC
C***DATE WRITTEN   770801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  CORRECTION FACTOR,LOG GAMMA,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the log Gamma correction factor so that
C            ALOG(GAMMA(X)) = ALOG(SQRT(2*PI)) + (X-.5)*ALOG(X) - X
C            + R9LGMC(X)
C***DESCRIPTION
C
C Compute the log gamma correction factor for X .GE. 10.0 so that
C  ALOG (GAMMA(X)) = ALOG(SQRT(2*PI)) + (X-.5)*ALOG(X) - X + R9LGMC(X)
C
C Series for ALGM       on the interval  0.          to  1.00000D-02
C                                        with weighted error   3.40E-16
C                                         log weighted error  15.47
C                               significant figures required  14.39
C                                    decimal places required  15.86
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL,INITS,R1MACH,XERROR
C***END PROLOGUE  R9LGMC
      DIMENSION ALGMCS(6)
      DATA ALGMCS( 1) /    .1666389480 45186E0 /
      DATA ALGMCS( 2) /   -.0000138494 817606E0 /
      DATA ALGMCS( 3) /    .0000000098 108256E0 /
      DATA ALGMCS( 4) /   -.0000000000 180912E0 /
      DATA ALGMCS( 5) /    .0000000000 000622E0 /
      DATA ALGMCS( 6) /   -.0000000000 000003E0 /
      DATA NALGM, XBIG, XMAX / 0, 2*0.0 /
C***FIRST EXECUTABLE STATEMENT  R9LGMC
      IF (NALGM.NE.0) GO TO 10
      NALGM = INITS (ALGMCS, 6, R1MACH(3))
      XBIG = 1.0/SQRT(R1MACH(3))
      XMAX = EXP (AMIN1(ALOG(R1MACH(2)/12.0), -ALOG(12.0*R1MACH(1))) )
C
 10   IF (X.LT.10.0) CALL XERROR ( 'R9LGMC  X MUST BE GE 10', 23, 1, 2)
      IF (X.GE.XMAX) GO TO 20
C
      R9LGMC = 1.0/(12.0*X)
      IF (X.LT.XBIG) R9LGMC = CSEVL (2.0*(10./X)**2-1., ALGMCS, NALGM)/X
      RETURN
C
 20   R9LGMC = 0.0
      CALL XERROR ( 'R9LGMC  X SO BIG R9LGMC UNDERFLOWS', 34, 2, 1)
      RETURN
C
      END
