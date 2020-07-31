      DOUBLE PRECISION FUNCTION DBESI0(X)
C***BEGIN PROLOGUE  DBESI0
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C10B1
C***KEYWORDS  BESSEL FUNCTION,DOUBLE PRECISION,FIRST KIND,
C             MODIFIED BESSEL FUNCTION,ORDER ZERO,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. hyperbolic Bessel function of the first
C            kind of order zero.
C***DESCRIPTION
C
C DBESI0(X) calculates the double precision modified (hyperbolic)
C Bessel function of the first kind of order zero and double
C precision argument X.
C
C Series for BI0        on the interval  0.          to  9.00000E+00
C                                        with weighted error   9.51E-34
C                                         log weighted error  33.02
C                               significant figures required  33.31
C                                    decimal places required  33.65
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DBSI0E,DCSEVL,INITDS,XERROR
C***END PROLOGUE  DBESI0
      DOUBLE PRECISION X, BI0CS(18), XMAX, XSML, Y, D1MACH,
     1  DCSEVL, DBSI0E
      DATA BI0 CS(  1) / -.7660547252 8391449510 8189497624 3285 D-1   /
      DATA BI0 CS(  2) / +.1927337953 9938082699 5240875088 1196 D+1   /
      DATA BI0 CS(  3) / +.2282644586 9203013389 3702929233 0415 D+0   /
      DATA BI0 CS(  4) / +.1304891466 7072904280 7933421069 1888 D-1   /
      DATA BI0 CS(  5) / +.4344270900 8164874513 7868268102 6107 D-3   /
      DATA BI0 CS(  6) / +.9422657686 0019346639 2317174411 8766 D-5   /
      DATA BI0 CS(  7) / +.1434006289 5106910799 6209187817 9957 D-6   /
      DATA BI0 CS(  8) / +.1613849069 6617490699 1541971999 4611 D-8   /
      DATA BI0 CS(  9) / +.1396650044 5356696994 9509270814 2522 D-10  /
      DATA BI0 CS( 10) / +.9579451725 5054453446 2752317189 3333 D-13  /
      DATA BI0 CS( 11) / +.5333981859 8625021310 1510774400 0000 D-15  /
      DATA BI0 CS( 12) / +.2458716088 4374707746 9678591999 9999 D-17  /
      DATA BI0 CS( 13) / +.9535680890 2487700269 4434133333 3333 D-20  /
      DATA BI0 CS( 14) / +.3154382039 7214273367 8933333333 3333 D-22  /
      DATA BI0 CS( 15) / +.9004564101 0946374314 6666666666 6666 D-25  /
      DATA BI0 CS( 16) / +.2240647369 1236700160 0000000000 0000 D-27  /
      DATA BI0 CS( 17) / +.4903034603 2428373333 3333333333 3333 D-30  /
      DATA BI0 CS( 18) / +.9508172606 1226666666 6666666666 6666 D-33  /
      DATA NTI0, XSML, XMAX / 0, 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DBESI0
      IF (NTI0.NE.0) GO TO 10
      NTI0 = INITDS (BI0CS, 18, 0.1*SNGL(D1MACH(3)))
      XSML = DSQRT (8.0D0*D1MACH(3))
      XMAX = DLOG (D1MACH(2))
C
 10   Y = DABS(X)
      IF (Y.GT.3.0D0) GO TO 20
C
      DBESI0 = 1.0D0
      IF (Y.GT.XSML) DBESI0 = 2.75D0 + DCSEVL (Y*Y/4.5D0-1.D0, BI0CS,
     1  NTI0)
      RETURN
C
 20   IF (Y.GT.XMAX) CALL XERROR ( 'DBESI0  DABS(X) SO BIG I0 OVERFLOWS'
     1, 35, 2, 2)
C
      DBESI0 = DEXP(Y) * DBSI0E(X)
C
      RETURN
      END