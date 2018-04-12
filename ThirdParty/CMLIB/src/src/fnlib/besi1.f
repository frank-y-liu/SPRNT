      FUNCTION BESI1(X)
C***BEGIN PROLOGUE  BESI1
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C10B1
C***KEYWORDS  BESSEL FUNCTION,FIRST KIND,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION,ORDER ONE
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the hyperbolic Bessel function of first kind of
C            order one
C***DESCRIPTION
C
C BESI1(X) calculates the modified (hyperbolic) Bessel function
C of the first kind of order one for real argument X.
C
C Series for BI1        on the interval  0.          to  9.00000D+00
C                                        with weighted error   2.40E-17
C                                         log weighted error  16.62
C                               significant figures required  16.23
C                                    decimal places required  17.14
C***REFERENCES  (NONE)
C***ROUTINES CALLED  BESI1E,CSEVL,INITS,R1MACH,XERROR
C***END PROLOGUE  BESI1
      DIMENSION BI1CS(11)
      DATA BI1 CS( 1) /   -.0019717132 61099859E0 /
      DATA BI1 CS( 2) /    .4073488766 7546481E0 /
      DATA BI1 CS( 3) /    .0348389942 99959456E0 /
      DATA BI1 CS( 4) /    .0015453945 56300123E0 /
      DATA BI1 CS( 5) /    .0000418885 21098377E0 /
      DATA BI1 CS( 6) /    .0000007649 02676483E0 /
      DATA BI1 CS( 7) /    .0000000100 42493924E0 /
      DATA BI1 CS( 8) /    .0000000000 99322077E0 /
      DATA BI1 CS( 9) /    .0000000000 00766380E0 /
      DATA BI1 CS(10) /    .0000000000 00004741E0 /
      DATA BI1 CS(11) /    .0000000000 00000024E0 /
      DATA NTI1, XMIN, XSML, XMAX / 0, 3*0. /
C***FIRST EXECUTABLE STATEMENT  BESI1
      IF (NTI1.NE.0) GO TO 10
      NTI1 = INITS (BI1CS, 11, 0.1*R1MACH(3))
      XMIN = 2.0*R1MACH(1)
      XSML = SQRT (8.0*R1MACH(3))
      XMAX = ALOG (R1MACH(2))
C
 10   Y = ABS(X)
      IF (Y.GT.3.0) GO TO 20
C
      BESI1 = 0.0
      IF (Y.EQ.0.0)  RETURN
C
      IF (Y.LT.XMIN) CALL XERROR ( 'BESI1   ABS(X) SO SMALL I1 UNDERFLOW
     1S', 37, 1, 1)
      IF (Y.GT.XMIN) BESI1 = 0.5*X
      IF (Y.GT.XSML) BESI1 = X * (.875 + CSEVL(Y*Y/4.5-1., BI1CS, NTI1))
      RETURN
C
 20   IF (Y.GT.XMAX) CALL XERROR ( 'BESI1   ABS(X) SO BIG I1 OVERFLOWS',
     1 34, 2, 2)
C
      BESI1 = EXP(Y) * BESI1E(X)
C
      RETURN
      END
