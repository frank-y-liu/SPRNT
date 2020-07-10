      FUNCTION BI(X)
C***BEGIN PROLOGUE  BI
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C10D
C***KEYWORDS  BAIRY FUNCTION,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the Bairy function.
C***DESCRIPTION
C
C BI(X) calculates the Airy function of the second kind for real
C argument X.
C
C Series for BIF        on the interval -1.00000D+00 to  1.00000D+00
C                                        with weighted error   1.88E-19
C                                         log weighted error  18.72
C                               significant figures required  17.74
C                                    decimal places required  19.20
C
C Series for BIG        on the interval -1.00000D+00 to  1.00000D+00
C                                        with weighted error   2.61E-17
C                                         log weighted error  16.58
C                               significant figures required  15.17
C                                    decimal places required  17.03
C
C Series for BIF2       on the interval  1.00000D+00 to  8.00000D+00
C                                        with weighted error   1.11E-17
C                                         log weighted error  16.95
C                        approx significant figures required  16.5
C                                    decimal places required  17.45
C
C Series for BIG2       on the interval  1.00000D+00 to  8.00000D+00
C                                        with weighted error   1.19E-18
C                                         log weighted error  17.92
C                        approx significant figures required  17.2
C                                    decimal places required  18.42
C***REFERENCES  (NONE)
C***ROUTINES CALLED  BIE,CSEVL,INITS,R1MACH,R9AIMP,XERROR
C***END PROLOGUE  BI
      DIMENSION BIFCS(9), BIGCS(8), BIF2CS(10), BIG2CS(10)
      DATA BIF CS( 1) /   -.0167302164 7198664948E0 /
      DATA BIF CS( 2) /    .1025233583 424944561E0 /
      DATA BIF CS( 3) /    .0017083092 5073815165E0 /
      DATA BIF CS( 4) /    .0000118625 4546774468E0 /
      DATA BIF CS( 5) /    .0000000449 3290701779E0 /
      DATA BIF CS( 6) /    .0000000001 0698207143E0 /
      DATA BIF CS( 7) /    .0000000000 0017480643E0 /
      DATA BIF CS( 8) /    .0000000000 0000020810E0 /
      DATA BIF CS( 9) /    .0000000000 0000000018E0 /
      DATA BIG CS( 1) /    .0224662232 4857452E0 /
      DATA BIG CS( 2) /    .0373647754 5301955E0 /
      DATA BIG CS( 3) /    .0004447621 8957212E0 /
      DATA BIG CS( 4) /    .0000024708 0756363E0 /
      DATA BIG CS( 5) /    .0000000079 1913533E0 /
      DATA BIG CS( 6) /    .0000000000 1649807E0 /
      DATA BIG CS( 7) /    .0000000000 0002411E0 /
      DATA BIG CS( 8) /    .0000000000 0000002E0 /
      DATA BIF2CS( 1) /   0.0998457269 3816041E0 /
      DATA BIF2CS( 2) /    .4786249778 63005538E0 /
      DATA BIF2CS( 3) /    .0251552119 604330118E0 /
      DATA BIF2CS( 4) /    .0005820693 885232645E0 /
      DATA BIF2CS( 5) /    .0000074997 659644377E0 /
      DATA BIF2CS( 6) /    .0000000613 460287034E0 /
      DATA BIF2CS( 7) /    .0000000003 462753885E0 /
      DATA BIF2CS( 8) /    .0000000000 014288910E0 /
      DATA BIF2CS( 9) /    .0000000000 000044962E0 /
      DATA BIF2CS(10) /    .0000000000 000000111E0 /
      DATA BIG2CS( 1) /    .0333056621 45514340E0 /
      DATA BIG2CS( 2) /    .1613092151 23197068E0 /
      DATA BIG2CS( 3) /    .0063190073 096134286E0 /
      DATA BIG2CS( 4) /    .0001187904 568162517E0 /
      DATA BIG2CS( 5) /    .0000013045 345886200E0 /
      DATA BIG2CS( 6) /    .0000000093 741259955E0 /
      DATA BIG2CS( 7) /    .0000000000 474580188E0 /
      DATA BIG2CS( 8) /    .0000000000 001783107E0 /
      DATA BIG2CS( 9) /    .0000000000 000005167E0 /
      DATA BIG2CS(10) /    .0000000000 000000011E0 /
      DATA NBIF, NBIG, NBIF2, NBIG2, X3SML, XMAX / 4*0, 2*0.0 /
C***FIRST EXECUTABLE STATEMENT  BI
      IF (NBIF.NE.0) GO TO 10
      ETA = 0.1*R1MACH(3)
      NBIF  = INITS (BIFCS , 9, ETA)
      NBIG  = INITS (BIGCS , 8, ETA)
      NBIF2 = INITS (BIF2CS, 10, ETA)
      NBIG2 = INITS (BIG2CS, 10, ETA)
C
      X3SML = ETA**0.3333
      XMAX = (1.5*ALOG(R1MACH(2)))**0.6666
C
 10   IF (X.GE.(-1.0)) GO TO 20
      CALL R9AIMP (X, XM, THETA)
      BI = XM * SIN(THETA)
      RETURN
C
 20   IF (X.GT.1.0) GO TO 30
      Z = 0.0
      IF (ABS(X).GT.X3SML) Z = X**3
      BI = 0.625 + CSEVL (Z, BIFCS, NBIF) + X*(0.4375 +
     1  CSEVL (Z, BIGCS, NBIG))
      RETURN
C
 30   IF (X.GT.2.0) GO TO 40
      Z = (2.0*X**3 - 9.0) / 7.0
      BI = 1.125 + CSEVL (Z, BIF2CS, NBIF2) + X*(0.625 +
     1  CSEVL (Z, BIG2CS, NBIG2))
      RETURN
C
 40   IF (X.GT.XMAX) CALL XERROR ( 'BI      X SO BIG THAT BI OVERFLOWS',
     1 34, 1, 2)
C
      BI = BIE(X) * EXP(2.0*X*SQRT(X)/3.0)
      RETURN
C
      END
