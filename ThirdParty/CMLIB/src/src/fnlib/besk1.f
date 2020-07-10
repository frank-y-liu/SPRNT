      FUNCTION BESK1(X)
C***BEGIN PROLOGUE  BESK1
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C10B1
C***KEYWORDS  BESSEL FUNCTION,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION,ORDER ONE,THIRD KIND
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the hyperbolic Bessel function of the third kind
C            of order one
C***DESCRIPTION
C
C BESK1(X) computes the modified (hyperbolic) Bessel function of third
C kind of order one for real argument X, where X .GT. 0.
C
C Series for BK1        on the interval  0.          to  4.00000D+00
C                                        with weighted error   7.02E-18
C                                         log weighted error  17.15
C                               significant figures required  16.73
C                                    decimal places required  17.67
C***REFERENCES  (NONE)
C***ROUTINES CALLED  BESI1,BESK1E,CSEVL,INITS,R1MACH,XERROR
C***END PROLOGUE  BESK1
      DIMENSION BK1CS(11)
      DATA BK1 CS( 1) /    .0253002273 389477705E0 /
      DATA BK1 CS( 2) /   -.3531559607 76544876E0 /
      DATA BK1 CS( 3) /   -.1226111808 22657148E0 /
      DATA BK1 CS( 4) /   -.0069757238 596398643E0 /
      DATA BK1 CS( 5) /   -.0001730288 957513052E0 /
      DATA BK1 CS( 6) /   -.0000024334 061415659E0 /
      DATA BK1 CS( 7) /   -.0000000221 338763073E0 /
      DATA BK1 CS( 8) /   -.0000000001 411488392E0 /
      DATA BK1 CS( 9) /   -.0000000000 006666901E0 /
      DATA BK1 CS(10) /   -.0000000000 000024274E0 /
      DATA BK1 CS(11) /   -.0000000000 000000070E0 /
      DATA NTK1, XMIN, XSML, XMAX /0, 0., 0., 0. /
C***FIRST EXECUTABLE STATEMENT  BESK1
      IF (NTK1.NE.0) GO TO 10
      NTK1 = INITS (BK1CS, 11, 0.1*R1MACH(3))
      XMIN = EXP (AMAX1(ALOG(R1MACH(1)), -ALOG(R1MACH(2))) + .01)
      XSML = SQRT (4.0*R1MACH(3))
      XMAX = -ALOG(R1MACH(1))
      XMAX = XMAX - 0.5*XMAX*ALOG(XMAX)/(XMAX+0.5)
C
 10   IF (X.LE.0.) CALL XERROR ( 'BESK1   X IS ZERO OR NEGATIVE', 29,
     1  2, 2)
      IF (X.GT.2.0) GO TO 20
C
      IF (X.LT.XMIN) CALL XERROR ( 'BESK1   X SO SMALL K1 OVERFLOWS',
     1  31, 3, 2)
      Y = 0.
      IF (X.GT.XSML) Y = X*X
      BESK1 = ALOG(0.5*X)*BESI1(X) +
     1  (0.75 + CSEVL (.5*Y-1., BK1CS, NTK1))/X
      RETURN
C
 20   BESK1 = 0.
      IF (X.GT.XMAX) CALL XERROR ( 'BESK1   X SO BIG K1 UNDERFLOWS',
     1  30, 1, 1)
      IF (X.GT.XMAX) RETURN
C
      BESK1 = EXP(-X) * BESK1E(X)
C
      RETURN
      END
