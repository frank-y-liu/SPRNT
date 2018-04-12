      FUNCTION PSI(X)
C***BEGIN PROLOGUE  PSI
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7C
C***KEYWORDS  DIGAMMA FUNCTION,PSI FUNCTION,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the Psi (or Digamma) function.
C***DESCRIPTION
C
C PSI(X) calculates the psi (or digamma) function for real argument X.
C PSI(X) is the logarithmic derivative of the gamma function of X.
C
C Series for PSI        on the interval  0.          to  1.00000D+00
C                                        with weighted error   2.03E-17
C                                         log weighted error  16.69
C                               significant figures required  16.39
C                                    decimal places required  17.37
C
C Series for APSI       on the interval  0.          to  2.50000D-01
C                                        with weighted error   5.54E-17
C                                         log weighted error  16.26
C                               significant figures required  14.42
C                                    decimal places required  16.86
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL,INITS,R1MACH,XERROR
C***END PROLOGUE  PSI
      DIMENSION PSICS(23), APSICS(16)
      DATA PSI CS( 1) /   -.0380570808 35217922E0 /
      DATA PSI CS( 2) /    .4914153930 2938713E0 /
      DATA PSI CS( 3) /   -.0568157478 21244730E0 /
      DATA PSI CS( 4) /    .0083578212 25914313E0 /
      DATA PSI CS( 5) /   -.0013332328 57994342E0 /
      DATA PSI CS( 6) /    .0002203132 87069308E0 /
      DATA PSI CS( 7) /   -.0000370402 38178456E0 /
      DATA PSI CS( 8) /    .0000062837 93654854E0 /
      DATA PSI CS( 9) /   -.0000010712 63908506E0 /
      DATA PSI CS(10) /    .0000001831 28394654E0 /
      DATA PSI CS(11) /   -.0000000313 53509361E0 /
      DATA PSI CS(12) /    .0000000053 72808776E0 /
      DATA PSI CS(13) /   -.0000000009 21168141E0 /
      DATA PSI CS(14) /    .0000000001 57981265E0 /
      DATA PSI CS(15) /   -.0000000000 27098646E0 /
      DATA PSI CS(16) /    .0000000000 04648722E0 /
      DATA PSI CS(17) /   -.0000000000 00797527E0 /
      DATA PSI CS(18) /    .0000000000 00136827E0 /
      DATA PSI CS(19) /   -.0000000000 00023475E0 /
      DATA PSI CS(20) /    .0000000000 00004027E0 /
      DATA PSI CS(21) /   -.0000000000 00000691E0 /
      DATA PSI CS(22) /    .0000000000 00000118E0 /
      DATA PSI CS(23) /   -.0000000000 00000020E0 /
      DATA APSICS( 1) /   -.0204749044 678185E0 /
      DATA APSICS( 2) /   -.0101801271 534859E0 /
      DATA APSICS( 3) /    .0000559718 725387E0 /
      DATA APSICS( 4) /   -.0000012917 176570E0 /
      DATA APSICS( 5) /    .0000000572 858606E0 /
      DATA APSICS( 6) /   -.0000000038 213539E0 /
      DATA APSICS( 7) /    .0000000003 397434E0 /
      DATA APSICS( 8) /   -.0000000000 374838E0 /
      DATA APSICS( 9) /    .0000000000 048990E0 /
      DATA APSICS(10) /   -.0000000000 007344E0 /
      DATA APSICS(11) /    .0000000000 001233E0 /
      DATA APSICS(12) /   -.0000000000 000228E0 /
      DATA APSICS(13) /    .0000000000 000045E0 /
      DATA APSICS(14) /   -.0000000000 000009E0 /
      DATA APSICS(15) /    .0000000000 000002E0 /
      DATA APSICS(16) /   -.0000000000 000000E0 /
      DATA PI     / 3.1415926535 8979324E0/
      DATA NTPSI, NTAPSI, XBIG, DXREL /0, 0, 0., 0./
C***FIRST EXECUTABLE STATEMENT  PSI
      IF (NTPSI.NE.0) GO TO 10
      NTPSI = INITS (PSICS, 23, 0.1*R1MACH(3))
      NTAPSI = INITS (APSICS, 16, 0.1*R1MACH(3))
C
      XBIG = 1.0/SQRT(R1MACH(3))
      DXREL = SQRT (R1MACH(4))
C
 10   Y = ABS(X)
      IF (Y.GE.2.0) GO TO 30
C
C PSI(X) FOR -2. .LT. X .LT. 2.
C
      N = X
      IF (X.LT.0.) N = N - 1
      Y = X - FLOAT(N)
      N = N - 1
      PSI = CSEVL (2.*Y-1., PSICS, NTPSI)
      IF (N.EQ.0) RETURN
C
      N = -N
      IF (X.EQ.0.) CALL XERROR ( 'PSI     X IS 0', 14, 2, 2)
      IF (X.LT.0. .AND. X+FLOAT(N-2).EQ.0.) CALL XERROR (  'PSI     X IS
     1 A NEGATIVE INTEGER', 31, 3, 2)
      IF (X.LT.(-0.5) .AND. ABS((X-AINT(X-0.5))/X).LT.DXREL) CALL XERROR
     1  (  'PSI     ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE
     2 INTEGER', 68, 1, 1)
C
      DO 20 I=1,N
        PSI = PSI - 1.0/(X+FLOAT(I-1))
 20   CONTINUE
      RETURN
C
C PSI(X) FOR ABS(X) .GE. 2.
C
 30   AUX = 0.
      IF (Y.LT.XBIG) AUX = CSEVL (8./Y**2-1., APSICS, NTAPSI)
      IF (X.LT.0.) PSI = ALOG(ABS(X)) - 0.5/X + AUX - PI*COT(PI*X)
      IF (X.GT.0.) PSI = ALOG(X) - 0.5/X + AUX
      RETURN
C
      END
