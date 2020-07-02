      FUNCTION ERF(X)
C***BEGIN PROLOGUE  ERF
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C8A,L5A1E
C***KEYWORDS  ERF,ERROR FUNCTION,SPECIAL FUNCTION 
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the error function, ERF. 
C***DESCRIPTION
C
C ERF(X) calculates the single precision error function for
C single precision argument X.
C
C Series for ERF        on the interval  0.          to  1.00000D+00
C                                        with weighted error   7.10E-18
C                                         log weighted error  17.15
C                               significant figures required  16.31
C                                    decimal places required  17.71
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL,ERFC,INITS,R1MACH
C***END PROLOGUE  ERF
      DIMENSION ERFCS(13)
      EXTERNAL ERFC 
      DATA ERF CS( 1) /   -.0490461212 34691808E0 /
      DATA ERF CS( 2) /   -.1422612051 0371364E0 /
      DATA ERF CS( 3) /    .0100355821 87599796E0 /
      DATA ERF CS( 4) /   -.0005768764 69976748E0 /
      DATA ERF CS( 5) /    .0000274199 31252196E0 /
      DATA ERF CS( 6) /   -.0000011043 17550734E0 /
      DATA ERF CS( 7) /    .0000000384 88755420E0 /
      DATA ERF CS( 8) /   -.0000000011 80858253E0 /
      DATA ERF CS( 9) /    .0000000000 32334215E0 /
      DATA ERF CS(10) /   -.0000000000 00799101E0 /
      DATA ERF CS(11) /    .0000000000 00017990E0 /
      DATA ERF CS(12) /   -.0000000000 00000371E0 /
      DATA ERF CS(13) /    .0000000000 00000007E0 /
      DATA SQRTPI /1.772453850 9055160E0/
      DATA NTERF, XBIG, SQEPS / 0, 0., 0./
C***FIRST EXECUTABLE STATEMENT  ERF
      IF (NTERF.NE.0) GO TO 10
      NTERF = INITS (ERFCS, 13, 0.1*R1MACH(3))
      XBIG = SQRT(-ALOG(SQRTPI*R1MACH(3)))
      SQEPS = SQRT(2.0*R1MACH(3))
C
 10   Y = ABS(X)
      IF (Y.GT.1.) GO TO 20
C
C ERF(X) = 1. - ERFC(X) FOR -1. .LE. X .LE. 1.
C
      IF (Y.LE.SQEPS) ERF = 2.0*X/SQRTPI
      IF (Y.GT.SQEPS) ERF = X*(1.0 + CSEVL(2.*X**2-1., ERFCS, NTERF)) 
      RETURN
C
C ERF(X) = 1. - ERFC(X) FOR  ABS(X) .GT. 1.
C
 20   IF (Y.LE.XBIG) ERF = SIGN (1.0-ERFC(Y), X)
      IF (Y.GT.XBIG) ERF = SIGN (1.0, X)
C
      RETURN
      END 
