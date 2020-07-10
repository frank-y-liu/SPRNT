      FUNCTION AI(X)
C***BEGIN PROLOGUE  AI
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C10D
C***KEYWORDS  AIRY FUNCTION,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the Airy function.
C***DESCRIPTION
C
C AI(X) computes the Airy function Ai(X)
C Series for AIF        on the interval -1.00000D+00 to  1.00000D+00
C                                        with weighted error   1.09E-19
C                                         log weighted error  18.96
C                               significant figures required  17.76
C                                    decimal places required  19.44
C
C Series for AIG        on the interval -1.00000D+00 to  1.00000D+00
C                                        with weighted error   1.51E-17
C                                         log weighted error  16.82
C                               significant figures required  15.19
C                                    decimal places required  17.27
C***REFERENCES  (NONE)
C***ROUTINES CALLED  AIE,CSEVL,INITS,R1MACH,R9AIMP,XERROR
C***END PROLOGUE  AI
      DIMENSION AIFCS(9), AIGCS(8)
      DATA AIF CS( 1) /   -.0379713584 9666999750E0 /
      DATA AIF CS( 2) /    .0591918885 3726363857E0 /
      DATA AIF CS( 3) /    .0009862928 0577279975E0 /
      DATA AIF CS( 4) /    .0000068488 4381907656E0 /
      DATA AIF CS( 5) /    .0000000259 4202596219E0 /
      DATA AIF CS( 6) /    .0000000000 6176612774E0 /
      DATA AIF CS( 7) /    .0000000000 0010092454E0 /
      DATA AIF CS( 8) /    .0000000000 0000012014E0 /
      DATA AIF CS( 9) /    .0000000000 0000000010E0 /
      DATA AIG CS( 1) /    .0181523655 8116127E0 /
      DATA AIG CS( 2) /    .0215725631 6601076E0 /
      DATA AIG CS( 3) /    .0002567835 6987483E0 /
      DATA AIG CS( 4) /    .0000014265 2141197E0 /
      DATA AIG CS( 5) /    .0000000045 7211492E0 /
      DATA AIG CS( 6) /    .0000000000 0952517E0 /
      DATA AIG CS( 7) /    .0000000000 0001392E0 /
      DATA AIG CS( 8) /    .0000000000 0000001E0 /
      DATA NAIF, NAIG, X3SML, XMAX / 2*0, 2*0.0 /
C***FIRST EXECUTABLE STATEMENT  AI
      IF (NAIF.NE.0) GO TO 10
      NAIF = INITS (AIFCS, 9, 0.1*R1MACH(3))
      NAIG = INITS (AIGCS, 8, 0.1*R1MACH(3))
C
      X3SML = R1MACH(3)**0.3334
      XMAX = (-1.5*ALOG(R1MACH(1)))**0.6667
      XMAX = XMAX - XMAX*ALOG(XMAX)/(4.0*SQRT(XMAX)+1.0) - 0.01
C
 10   IF (X.GE.(-1.0)) GO TO 20
      CALL R9AIMP (X, XM, THETA)
      AI = XM * COS(THETA)
      RETURN
C
 20   IF (X.GT.1.0) GO TO 30
      Z = 0.0
      IF (ABS(X).GT.X3SML) Z = X**3
      AI = 0.375 + (CSEVL (Z, AIFCS, NAIF) - X*(0.25 +
     1  CSEVL (Z, AIGCS, NAIG)) )
      RETURN
C
 30   IF (X.GT.XMAX) GO TO 40
      AI = AIE(X) * EXP(-2.0*X*SQRT(X)/3.0)
      RETURN
C
 40   AI = 0.0
      CALL XERROR ( 'AI      X SO BIG AI UNDERFLOWS', 30, 1, 1)
      RETURN
C
      END
