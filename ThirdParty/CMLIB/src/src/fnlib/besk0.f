      FUNCTION BESK0(X)
C***BEGIN PROLOGUE  BESK0
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C10B1
C***KEYWORDS  BESSEL FUNCTION,HYBERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION,ORDER ZERO,THIRD KIND
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the hyperbolic Bessel function of the third kind
C            of order zero
C***DESCRIPTION
C
C BESK0(X) calculates the modified (hyperbolic) Bessel function
C of the third kind of order zero for real argument X .GT. 0.0.
C
C Series for BK0        on the interval  0.          to  4.00000D+00
C                                        with weighted error   3.57E-19
C                                         log weighted error  18.45
C                               significant figures required  17.99
C                                    decimal places required  18.97
C***REFERENCES  (NONE)
C***ROUTINES CALLED  BESI0,BESK0E,CSEVL,INITS,R1MACH,XERROR
C***END PROLOGUE  BESK0
      DIMENSION BK0CS(11)
      DATA BK0 CS( 1) /   -.0353273932 3390276872E0 /
      DATA BK0 CS( 2) /    .3442898999 246284869E0 /
      DATA BK0 CS( 3) /    .0359799365 1536150163E0 /
      DATA BK0 CS( 4) /    .0012646154 1144692592E0 /
      DATA BK0 CS( 5) /    .0000228621 2103119451E0 /
      DATA BK0 CS( 6) /    .0000002534 7910790261E0 /
      DATA BK0 CS( 7) /    .0000000019 0451637722E0 /
      DATA BK0 CS( 8) /    .0000000000 1034969525E0 /
      DATA BK0 CS( 9) /    .0000000000 0004259816E0 /
      DATA BK0 CS(10) /    .0000000000 0000013744E0 /
      DATA BK0 CS(11) /    .0000000000 0000000035E0 /
      DATA NTK0, XSML, XMAX / 0, 0., 0. /
C***FIRST EXECUTABLE STATEMENT  BESK0
      IF (NTK0.NE.0) GO TO 10
      NTK0 = INITS (BK0CS, 11, 0.1*R1MACH(3))
      XSML = SQRT (4.0*R1MACH(3))
      XMAX = -ALOG(R1MACH(1))
      XMAX = XMAX - 0.5*XMAX*ALOG(XMAX)/(XMAX+0.5) - 0.01
C
 10   IF (X.LE.0.) CALL XERROR ( 'BESK0   X IS ZERO OR NEGATIVE', 29,
     1  2, 2)
      IF (X.GT.2.) GO TO 20
C
      Y = 0.
      IF (X.GT.XSML) Y = X*X
      BESK0 = -ALOG(0.5*X)*BESI0(X) - .25 + CSEVL (.5*Y-1., BK0CS, NTK0)
      RETURN
C
 20   BESK0 = 0.
      IF (X.GT.XMAX) CALL XERROR ( 'BESK0   X SO BIG K0 UNDERFLOWS', 30,
     1  1, 1)
      IF (X.GT.XMAX) RETURN
C
      BESK0 = EXP(-X) * BESK0E(X)
C
      RETURN
      END
