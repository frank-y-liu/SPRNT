      FUNCTION POCH1(A,X)
C***BEGIN PROLOGUE  POCH1
C***DATE WRITTEN   770801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C1,C7A
C***KEYWORDS  FIRST ORDER,POCHHAMMER,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes Pochhammer's symbol from first order.
C***DESCRIPTION
C
C Evaluate a generalization of Pochhammer's symbol for special
C situations that require especially accurate values when X is small in
C        POCH1(A,X) = (POCH(A,X)-1)/X
C                   = (GAMMA(A+X)/GAMMA(A) - 1.0)/X .
C This specification is particularly suited for stably computing
C expressions such as
C        (GAMMA(A+X)/GAMMA(A) - GAMMA(B+X)/GAMMA(B))/X
C             = POCH1(A,X) - POCH1(B,X)
C Note that POCH1(A,0.0) = PSI(A)
C
C When ABS(X) is so small that substantial cancellation will occur if
C the straightforward formula is used, we  use an expansion due
C to Fields and discussed by Y. L. Luke, The Special Functions and Their
C Approximations, Vol. 1, Academic Press, 1969, page 34.
C
C The ratio POCH(A,X) = GAMMA(A+X)/GAMMA(A) is written by Luke as
C        (A+(X-1)/2)**X * polynomial in (A+(X-1)/2)**(-2) .
C In order to maintain significance in POCH1, we write for positive A
C        (A+(X-1)/2)**X = EXP(X*ALOG(A+(X-1)/2)) = EXP(Q)
C                       = 1.0 + Q*EXPREL(Q) .
C Likewise the polynomial is written
C        POLY = 1.0 + X*POLY1(A,X) .
C Thus,
C        POCH1(A,X) = (POCH(A,X) - 1) / X
C                   = EXPREL(Q)*(Q/X + Q*POLY1(A,X)) + POLY1(A,X)
C***REFERENCES  (NONE)
C***ROUTINES CALLED  EXPREL,POCH,PSI,R1MACH,XERROR
C***END PROLOGUE  POCH1
      DIMENSION BERN(9), GBERN(10)
      DATA BERN( 1) /   .8333333333 3333333E-01 /
      DATA BERN( 2) /  -.1388888888 8888889E-02 /
      DATA BERN( 3) /   .3306878306 8783069E-04 /
      DATA BERN( 4) /  -.8267195767 1957672E-06 /
      DATA BERN( 5) /   .2087675698 7868099E-07 /
      DATA BERN( 6) /  -.5284190138 6874932E-09 /
      DATA BERN( 7) /   .1338253653 0684679E-10 /
      DATA BERN( 8) /  -.3389680296 3225829E-12 /
      DATA BERN( 9) /   .8586062056 2778446E-14 /
      DATA PI / 3.1415926535 8979324 E0 /
      DATA SQTBIG, ALNEPS / 2*0.0 /
C***FIRST EXECUTABLE STATEMENT  POCH1
      IF (SQTBIG.NE.0.0) GO TO 10
      SQTBIG = 1.0/SQRT(24.0*R1MACH(1))
      ALNEPS = ALOG(R1MACH(3))
C
 10   IF (X.EQ.0.0) POCH1 = PSI(A)
      IF (X.EQ.0.0) RETURN
C
      ABSX = ABS(X)
      ABSA = ABS(A)
      IF (ABSX.GT.0.1*ABSA) GO TO 70
      IF (ABSX*ALOG(AMAX1(ABSA,2.0)).GT.0.1) GO TO 70
C
      BP = A
      IF (A.LT.(-0.5)) BP = 1.0 - A - X
      INCR = 0
      IF (BP.LT.10.0) INCR = 11.0 - BP
      B = BP + FLOAT(INCR)
C
      VAR = B + 0.5*(X-1.0)
      ALNVAR = ALOG(VAR)
      Q = X*ALNVAR
C
      POLY1 = 0.0
      IF (VAR.GE.SQTBIG) GO TO 40
      VAR2 = (1.0/VAR)**2
C
      RHO = 0.5*(X+1.0)
      GBERN(1) = 1.0
      GBERN(2) = -RHO/12.0
      TERM = VAR2
      POLY1 = GBERN(2)*TERM
C
      NTERMS = -0.5*ALNEPS/ALNVAR + 1.0
      IF (NTERMS.GT.9) CALL XERROR ( 'POCH1   NTERMS IS TOO BIG, MAYBE R
     11MACH(3) IS BAD', 49, 1, 2)
      IF (NTERMS.LT.2) GO TO 40
C
      DO 30 K=2,NTERMS
        GBK = 0.0
        DO 20 J=1,K
          NDX = K - J + 1
          GBK = GBK + BERN(NDX)*GBERN(J)
 20     CONTINUE
        GBERN(K+1) = -RHO*GBK/FLOAT(K)
C
        TERM = TERM * (FLOAT(2*K-2)-X)*(FLOAT(2*K-1)-X)*VAR2
        POLY1 = POLY1 + GBERN(K+1)*TERM
 30   CONTINUE
C
 40   POLY1 = (X-1.0)*POLY1
      POCH1 = EXPREL(Q)*(ALNVAR + Q*POLY1) + POLY1
C
      IF (INCR.EQ.0) GO TO 60
C
C WE HAVE POCH1(B,X).  BUT BP IS SMALL, SO WE USE BACKWARDS RECURSION
C TO OBTAIN POCH1(BP,X).
C
      DO 50 II=1,INCR
        I = INCR - II
        BINV = 1.0/(BP+FLOAT(I))
        POCH1 = (POCH1-BINV)/(1.0+X*BINV)
 50   CONTINUE
C
 60   IF (BP.EQ.A) RETURN
C
C WE HAVE POCH1(BP,X), BUT A IS LT -0.5.  WE THEREFORE USE A REFLECTION
C FORMULA TO OBTAIN POCH1(A,X).
C
      SINPXX = SIN(PI*X)/X
      SINPX2 = SIN(0.5*PI*X)
      TRIG = SINPXX*COT(PI*B) - 2.0*SINPX2*(SINPX2/X)
C
      POCH1 = TRIG + (1.0 + X*TRIG) * POCH1
      RETURN
C
 70   POCH1 = (POCH(A,X) - 1.0) / X
      RETURN
C
      END
