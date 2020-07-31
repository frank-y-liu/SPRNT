      SUBROUTINE XDLEGF(DNU1,NUDIFF,MU1,MU2,THETA,ID,PQA,IPQA)
C***BEGIN PROLOGUE   XDLEGF
C***DATE WRITTEN   820728   (YYMMDD)
C***REVISION DATE  871020   (YYMMDD)
C***CATEGORY NO.  C3a2,C9
C***KEYWORDS  LEGENDRE FUNCTIONS
C***AUTHOR  SMITH, JOHN M. (NBS AND GEORGE MASON UNIVERSITY)
C***PURPOSE  TO COMPUTE THE VALUES OF LEGENDRE FUNCTIONS
C***DESCRIPTION
C
C   XDLEGF: Extended-range Double-precision Legendre Functions
C
C   A feature of the XDLEGF subroutine for Legendre functions is
C the use of extended-range arithmetic, a software extension of
C ordinary floating-point arithmetic that greatly increases the
C exponent range of the representable numbers. This avoids the
C need for scaling the solutions to lie within the exponent range
C of the most restrictive manufacturer's hardware. The increased
C exponent range is achieved by allocating an integer storage
C location together with each floating-point storage location.
C
C   The interpretation of the pair (X,I) where X is floating-point
C and I is integer is X*(IR**I) where IR is the internal radix of
C the computer arithmetic.
C
C   This subroutine computes one of the following vectors:
C
C 1. Legendre function of the first kind of negative order, either
C    a. P(-MU1,NU,X), P(-MU1-1,NU,X), ..., P(-MU2,NU,X) or
C    b. P(-MU,NU1,X), P(-MU,NU1+1,X), ..., P(-MU,NU2,X)
C 2. Legendre function of the second kind, either
C    a. Q(MU1,NU,X), Q(MU1+1,NU,X), ..., Q(MU2,NU,X) or
C    b. Q(MU,NU1,X), Q(MU,NU1+1,X), ..., Q(MU,NU2,X)
C 3. Legendre function of the first kind of positive order, either
C    a. P(MU1,NU,X), P(MU1+1,NU,X), ..., P(MU2,NU,X) or
C    b. P(MU,NU1,X), P(MU,NU1+1,X), ..., P(MU,NU2,X)
C 4. Normalized Legendre polynomials, either
C    a. PN(MU1,NU,X), PN(MU1+1,NU,X), ..., PN(MU2,NU,X) or
C    b. PN(MU,NU1,X), PN(MU,NU1+1,X), ..., PN(MU,NU2,X)
C
C where X = COS(THETA).
C
C   The input values to XDLEGF are DNU1, NUDIFF, MU1, MU2, THETA,
C and ID. These must satisfy
C
C    DNU1 is DOUBLE PRECISION and greater than or equal to -0.5;
C    NUDIFF is INTEGER and non-negative;
C    MU1 is INTEGER and non-negative;
C    MU2 is INTEGER and greater than or equal to MU1;
C    THETA is DOUBLE PRECISION and in the half-open interval (0,PI/2];
C    ID is INTEGER and equal to 1, 2, 3 or 4;
C
C and  additionally either NUDIFF = 0 or MU2 = MU1.
C
C   If ID=1 and NUDIFF=0, a vector of type 1a above is computed
C with NU=DNU1.
C
C   If ID=1 and MU1=MU2, a vector of type 1b above is computed
C with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1.
C
C   If ID=2 and NUDIFF=0, a vector of type 2a above is computed
C with NU=DNU1.
C
C   If ID=2 and MU1=MU2, a vector of type 2b above is computed
C with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1.
C
C   If ID=3 and NUDIFF=0, a vector of type 3a above is computed
C with NU=DNU1.
C
C   If ID=3 and MU1=MU2, a vector of type 3b above is computed
C with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1.
C
C   If ID=4 and NUDIFF=0, a vector of type 4a above is computed
C with NU=DNU1.
C
C   If ID=4 and MU1=MU2, a vector of type 4b above is computed
C with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1.
C
C   In each case the vector of computed Legendre function values
C is returned in the extended-range vector (PQA(I),IPQA(I)). The
C length of this vector is either MU2-MU1+1 or NUDIFF+1.
C
C Where possible, XDLEGF returns IPQA(I) as zero. In this case the
C value of the Legendre function is contained entirely in PQA(I),
C so it can be used in subsequent computations without further
C consideration of extended-range arithmetic. If IPQA(I) is nonzero,
C then the value of the Legendre function is not representable in
C floating-point because of underflow or overflow. The program that
C calls XDLEGF must test IPQA(I) to ensure correct usage.
C
C***REFERENCES  OLVER AND SMITH,J.COMPUT.PHYSICS,51(1983),N0.3,502-518.
C***ROUTINES CALLED  XERROR, XDPMU, XDPMUP, XDPNRM, XDPQNU, XDQMU, XDQNU
C                    XDRED, XDSET
C***END PROLOGUE  XDLEGF
      DOUBLE PRECISION PQA,DNU1,DNU2,SX,THETA,X,PI2
      DIMENSION PQA(*),IPQA(*)
C
C***FIRST EXECUTABLE STATEMENT   XDLEGF
      CALL XDSET (0, 0, 0.0D0, 0)
      PI2=2.D0*DATAN(1.D0)
C
C        ZERO OUTPUT ARRAYS
C
      L=(MU2-MU1)+NUDIFF+1
      DO 290 I=1,L
      PQA(I)=0.
  290 IPQA(I)=0
C
C        CHECK FOR VALID INPUT VALUES
C
      IF(NUDIFF.LT.0) GO TO 400
      IF(DNU1.LT.-.5D0) GO TO 400
      IF(MU2.LT.MU1) GO TO 400
      IF(MU1.LT.0) GO TO 400
      IF(THETA.LE.0.D0.OR.THETA.GT.PI2) GO TO 420
      IF(ID.LT.1.OR.ID.GT.4) GO TO 400
      IF((MU1.NE.MU2).AND.(NUDIFF.GT.0)) GO TO 400
C
C        IF DNU1 IS NOT AN INTEGER, NORMALIZED P(MU,DNU,X)
C        CANNOT BE CALCULATED.  IF DNU1 IS AN INTEGER AND
C        MU1.GT.DNU2 THEN ALL VALUES OF P(+MU,DNU,X) AND
C        NORMALIZED P(MU,NU,X) WILL BE ZERO.
C
      DNU2=DNU1+NUDIFF
      IF((ID.EQ.3).AND.(DMOD(DNU1,1.D0).NE.0.)) GO TO 295
      IF((ID.EQ.4).AND.(DMOD(DNU1,1.D0).NE.0.)) GO TO 400
      IF((ID.EQ.3.OR.ID.EQ.4).AND.DBLE(FLOAT(MU1)).GT.DNU2) RETURN
  295 CONTINUE
C
      X=DCOS(THETA)
      SX=1.D0/DSIN(THETA)
      IF(ID.EQ.2) GO TO 300
      IF(MU2-MU1.LE.0) GO TO 360
C
C        FIXED NU, VARIABLE MU
C        CALL XDPMU TO CALCULATE P(-MU1,NU,X),....,P(-MU2,NU,X)
C
      CALL XDPMU(DNU1,DNU2,MU1,MU2,THETA,X,SX,ID,PQA,IPQA)
      GO TO 380
C
  300 IF(MU2.EQ.MU1) GO TO 320
C
C        FIXED NU, VARIABLE MU
C        CALL XDQMU TO CALCULATE Q(MU1,NU,X),....,Q(MU2,NU,X)
C
      CALL XDQMU(DNU1,DNU2,MU1,MU2,THETA,X,SX,ID,PQA,IPQA)
      GO TO 390
C
C        FIXED MU, VARIABLE NU
C        CALL XDQNU TO CALCULATE Q(MU,DNU1,X),....,Q(MU,DNU2,X)
C
  320 CALL XDQNU(DNU1,DNU2,MU1,THETA,X,SX,ID,PQA,IPQA)
      GO TO 390
C
C        FIXED MU, VARIABLE NU
C        CALL XDPQNU TO CALCULATE P(-MU,DNU1,X),....,P(-MU,DNU2,X)
C
  360 CALL XDPQNU(DNU1,DNU2,MU1,THETA,ID,PQA,IPQA)
C
C        IF ID = 3, TRANSFORM P(-MU,NU,X) VECTOR INTO
C        P(MU,NU,X) VECTOR.
C
  380 IF(ID.EQ.3) CALL XDPMUP(DNU1,DNU2,MU1,MU2,PQA,IPQA)
C
C        IF ID = 4, TRANSFORM P(-MU,NU,X) VECTOR INTO
C        NORMALIZED P(MU,NU,X) VECTOR.
C
      IF(ID.EQ.4) CALL XDPNRM(DNU1,DNU2,MU1,MU2,PQA,IPQA)
C
C        PLACE RESULTS IN REDUCED FORM IF POSSIBLE
C        AND RETURN TO MAIN PROGRAM.
C
  390 DO 395 I=1,L
      CALL XDRED(PQA(I),IPQA(I))
  395 CONTINUE
      RETURN
C
C        *****     ERROR TERMINATION     *****
C
  400 CALL XERROR('XDLEGF: DNU1,NUDIFF,MU1,MU2, or ID not valid',44,1,1)
      GO TO 430
  420 CALL XERROR('XDLEGF: THETA out of range',26,2,1)
  430 RETURN
      END
