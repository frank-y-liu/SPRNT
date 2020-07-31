      SUBROUTINE DBSGQ8(FUN,XT,BC,N,KK,ID,A,B,INBV,ERR,ANS,IERR,WORK)
C***BEGIN PROLOGUE  DBSGQ8
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***REFER TO  DBFQAD
C
C     Written by R.E. Jones and modified by D.E. Amos
C
C     Abstract    **** A DOUBLE PRECISION routine ****
C
C        DBSGQ8, a modification of GAUS8, integrates the
C        product of FUN(X) by the ID-th derivative of a spline
C        DBVALU(XT,BC,N,KK,ID,X,INBV,WORK)  between limits A and B.
C
C        DBSGQ8 calls DBVALU, DINTRV, I1MACH, D1MACH, XERROR
C
C     Description of Arguments
C
C        INPUT-- FUN,XT,BC,A,B,ERR are DOUBLE PRECISION
C        FUN - Name of external function of one argument which
C              multiplies DBVALU.
C        XT  - Knot array for DBVALU
C        BC  - B-coefficient array for DBVALU
C        N   - Number of B-coefficients for DBVALU
C        KK  - Order of the spline, KK.GE.1
C        ID  - Order of the spline derivative, 0.LE.ID.LE.KK-1
C        A   - Lower limit of integral
C        B   - Upper limit of integral (may be less than A)
C        INBV- Initialization parameter for DBVALU
C        ERR - Is a requested pseudorelative error tolerance.  Normally
C              pick a value of DABS(ERR).LT.1D-3.  ANS will normally
C              have no more error than DABS(ERR) times the integral of
C              the absolute value of FUN(X)*DBVALU(XT,BC,N,KK,X,ID,
C              INBV,WORK).
C
C
C        OUTPUT-- ERR,ANS,WORK are DOUBLE PRECISION
C        ERR - Will be an estimate of the absolute error in ANS if the
C              input value of ERR was negative.  (ERR is unchanged if
C              the input value of ERR was nonnegative.)  The estimated
C              error is solely for information to the user and should
C              not be used as a correction to the computed integral.
C        ANS - Computed value of integral
C        IERR- A status code
C            --Normal Codes
C               1 ANS most likely meets requested error tolerance,
C                 or A=B.
C              -1 A and B are too nearly equal to allow normal
C                 integration.  ANS is set to zero.
C            --Abnormal Code
C               2 ANS probably does not meet requested error tolerance.
C        WORK- Work vector of length 3*K for DBVALU
C***ROUTINES CALLED  D1MACH,DBVALU,I1MACH,XERROR
C***END PROLOGUE  DBSGQ8
C
      INTEGER ICALL,ID,IERR,INBV,K, KK, KML, KMX, L, LMN, LMX, LR, MXL,
     1 N, NBITS, NIB, NLMN, NLMX
      INTEGER I1MACH
      DOUBLE PRECISION A,AA,AE,ANIB,ANS,AREA,B,BC,C,CE,EE,EF,EPS,ERR,
     1 EST,GL,GLR,GR,HH,SQ2,TOL,VL,VR,WORK,W1, W2, W3, W4, XT, X1,
     2 X2, X3, X4, X, H
      DOUBLE PRECISION D1MACH, DBVALU, G8, FUN
      DIMENSION XT(*), BC(*), WORK(*)
      DIMENSION AA(60), HH(60), LR(60), VL(60), GR(60)
      DATA X1, X2, X3, X4/
     1     1.83434642495649805D-01,     5.25532409916328986D-01,
     2     7.96666477413626740D-01,     9.60289856497536232D-01/
      DATA W1, W2, W3, W4/
     1     3.62683783378361983D-01,     3.13706645877887287D-01,
     2     2.22381034453374471D-01,     1.01228536290376259D-01/
      DATA ICALL  /  0  /
      DATA SQ2/1.41421356D0/
      DATA NLMN/1/,KMX/5000/,KML/6/
      G8(X,H)=H*((W1*(FUN(X-X1*H)*DBVALU(XT,BC,N,KK,ID,X-X1*H,INBV,WORK)
     1            + FUN(X+X1*H)*DBVALU(XT,BC,N,KK,ID,X+X1*H,INBV,WORK))
     2          +W2*(FUN(X-X2*H)*DBVALU(XT,BC,N,KK,ID,X-X2*H,INBV,WORK)+
     3              FUN(X+X2*H)*DBVALU(XT,BC,N,KK,ID,X+X2*H,INBV,WORK)))
     4         +(W3*(FUN(X-X3*H)*DBVALU(XT,BC,N,KK,ID,X-X3*H,INBV,WORK)+
     5               FUN(X+X3*H)*DBVALU(XT,BC,N,KK,ID,X+X3*H,INBV,WORK))
     6          +W4*(FUN(X-X4*H)*DBVALU(XT,BC,N,KK,ID,X-X4*H,INBV,WORK)+
     7             FUN(X+X4*H)*DBVALU(XT,BC,N,KK,ID,X+X4*H,INBV,WORK))))
C
C     INITIALIZE
C
C***FIRST EXECUTABLE STATEMENT  DBSGQ8
      IF (ICALL.NE.0)  CALL XERROR(  'DBSGQ8- DBSGQ8 CALLED RECURSIVELY.
     1  RECURSIVE CALLS ARE ILLEGAL IN FORTRAN.', 75, 7, 2)
      ICALL = 1
      K = I1MACH(14)
      ANIB = D1MACH(5)*DBLE(FLOAT(K))/0.30102000D0
      NBITS = INT(SNGL(ANIB))
      NLMX = MIN0((NBITS*5)/8,60)
      ANS = 0.0D0
      IERR = 1
      CE = 0.0D0
      IF (A.EQ.B) GO TO 140
      LMX = NLMX
      LMN = NLMN
      IF (B.EQ.0.0D0) GO TO 10
      IF (DSIGN(1.0D0,B)*A.LE.0.0D0) GO TO 10
      C = DABS(1.0D0-A/B)
      IF (C.GT.0.1D0) GO TO 10
      IF (C.LE.0.0D0) GO TO 140
      ANIB = 0.5D0 - DLOG(C)/0.69314718D0
      NIB = INT(SNGL(ANIB))
      LMX = MIN0(NLMX,NBITS-NIB-7)
      IF (LMX.LT.1) GO TO 130
      LMN = MIN0(LMN,LMX)
   10 TOL = DMAX1(DABS(ERR),2.0D0**(5-NBITS))/2.0D0
      IF (ERR.EQ.0.0D0) TOL = DSQRT(D1MACH(4))
      EPS = TOL
      HH(1) = (B-A)/4.0D0
      AA(1) = A
      LR(1) = 1
      L = 1
      EST = G8(AA(L)+2.0D0*HH(L),2.0D0*HH(L))
      K = 8
      AREA = DABS(EST)
      EF = 0.5D0
      MXL = 0
C
C     COMPUTE REFINED ESTIMATES, ESTIMATE THE ERROR, ETC.
C
   20 GL = G8(AA(L)+HH(L),HH(L))
      GR(L) = G8(AA(L)+3.0D0*HH(L),HH(L))
      K = K + 16
      AREA = AREA + (DABS(GL)+DABS(GR(L))-DABS(EST))
      GLR = GL + GR(L)
      EE = DABS(EST-GLR)*EF
      AE = DMAX1(EPS*AREA,TOL*DABS(GLR))
      IF (EE-AE) 40, 40, 50
   30 MXL = 1
   40 CE = CE + (EST-GLR)
      IF (LR(L)) 60, 60, 80
C
C     CONSIDER THE LEFT HALF OF THIS LEVEL
C
   50 IF (K.GT.KMX) LMX = KML
      IF (L.GE.LMX) GO TO 30
      L = L + 1
      EPS = EPS*0.5D0
      EF = EF/SQ2
      HH(L) = HH(L-1)*0.5D0
      LR(L) = -1
      AA(L) = AA(L-1)
      EST = GL
      GO TO 20
C
C     PROCEED TO RIGHT HALF AT THIS LEVEL
C
   60 VL(L) = GLR
   70 EST = GR(L-1)
      LR(L) = 1
      AA(L) = AA(L) + 4.0D0*HH(L)
      GO TO 20
C
C     RETURN ONE LEVEL
C
   80 VR = GLR
   90 IF (L.LE.1) GO TO 120
      L = L - 1
      EPS = EPS*2.0D0
      EF = EF*SQ2
      IF (LR(L)) 100, 100, 110
  100 VL(L) = VL(L+1) + VR
      GO TO 70
  110 VR = VL(L+1) + VR
      GO TO 90
C
C      EXIT
C
  120 ANS = VR
      IF ((MXL.EQ.0) .OR. (DABS(CE).LE.2.0D0*TOL*AREA)) GO TO 140
      IERR = 2
      CALL XERROR( 'DBSGQ8- ANS IS PROBABLY INSUFFICIENTLY ACCURATE.',
     1 48, 3, 1)
      GO TO 140
  130 IERR = -1
      CALL XERROR( 'DBSGQ8- THE FOLLOWING TEMPORARY DIAGNOSTIC WILL APPE
     1AR ONLY ONCE.  A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL INTEGR
     2ATION.  ANS IS SET TO ZERO, AND IERR=-1.', 158, 1, -1)
  140 ICALL = 0
      IF (ERR.LT.0.0D0) ERR = CE
      RETURN
      END
