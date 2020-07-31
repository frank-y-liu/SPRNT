      FUNCTION ALNREL(X)
C***BEGIN PROLOGUE  ALNREL
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C4B
C***KEYWORDS  ELEMENTARY FUNCTION,LOGARITHM,RELATIVE
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Evaluates ln(1+X) accurate in the sense of relative error.
C***DESCRIPTION
C
C ALNREL(X) evaluates ln(1+X) accurately in the sense of relative
C error when X is very small.  This routine must be used to
C maintain relative error accuracy whenever X is small and
C accurately known.
C
C Series for ALNR       on the interval -3.75000D-01 to  3.75000D-01
C                                        with weighted error   1.93E-17
C                                         log weighted error  16.72
C                               significant figures required  16.44
C                                    decimal places required  17.40
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL,INITS,R1MACH,XERROR
C***END PROLOGUE  ALNREL
      DIMENSION ALNRCS(23)
      DATA ALNRCS( 1) /   1.0378693562 743770E0 /
      DATA ALNRCS( 2) /   -.1336430150 4908918E0 /
      DATA ALNRCS( 3) /    .0194082491 35520563E0 /
      DATA ALNRCS( 4) /   -.0030107551 12753577E0 /
      DATA ALNRCS( 5) /    .0004869461 47971548E0 /
      DATA ALNRCS( 6) /   -.0000810548 81893175E0 /
      DATA ALNRCS( 7) /    .0000137788 47799559E0 /
      DATA ALNRCS( 8) /   -.0000023802 21089435E0 /
      DATA ALNRCS( 9) /    .0000004164 04162138E0 /
      DATA ALNRCS(10) /   -.0000000735 95828378E0 /
      DATA ALNRCS(11) /    .0000000131 17611876E0 /
      DATA ALNRCS(12) /   -.0000000023 54670931E0 /
      DATA ALNRCS(13) /    .0000000004 25227732E0 /
      DATA ALNRCS(14) /   -.0000000000 77190894E0 /
      DATA ALNRCS(15) /    .0000000000 14075746E0 /
      DATA ALNRCS(16) /   -.0000000000 02576907E0 /
      DATA ALNRCS(17) /    .0000000000 00473424E0 /
      DATA ALNRCS(18) /   -.0000000000 00087249E0 /
      DATA ALNRCS(19) /    .0000000000 00016124E0 /
      DATA ALNRCS(20) /   -.0000000000 00002987E0 /
      DATA ALNRCS(21) /    .0000000000 00000554E0 /
      DATA ALNRCS(22) /   -.0000000000 00000103E0 /
      DATA ALNRCS(23) /    .0000000000 00000019E0 /
      DATA NLNREL, XMIN /0, 0./
C***FIRST EXECUTABLE STATEMENT  ALNREL
      IF (NLNREL.NE.0) GO TO 10
      NLNREL = INITS (ALNRCS, 23, 0.1*R1MACH(3))
      XMIN = -1.0 + SQRT(R1MACH(4))
C
 10   IF (X.LE.(-1.0)) CALL XERROR ( 'ALNREL  X IS LE -1', 18, 2, 2)
      IF (X.LT.XMIN) CALL XERROR ( 'ALNREL  ANSWER LT HALF PRECISION BEC
     1AUSE X TOO NEAR -1', 54,    1, 1)
C
      IF (ABS(X).LE.0.375) ALNREL = X*(1. -
     1  X*CSEVL (X/.375, ALNRCS, NLNREL))
      IF (ABS(X).GT.0.375) ALNREL = ALOG (1.0+X)
C
      RETURN
      END
