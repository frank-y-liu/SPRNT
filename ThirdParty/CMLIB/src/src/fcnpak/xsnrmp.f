      SUBROUTINE XSNRMP(NU,MU1,MU2,SARG,MODE,SPN,IPN,ISIG)
C***BEGIN PROLOGUE  XSNRMP
C***DATE WRITTEN   820712   (YYMMDD)
C***REVISION DATE  871110   (YYMMDD)
C***CATEGORY NO.  C3a2,C9
C***KEYWORDS  LEGENDRE FUNCTIONS
C***AUTHOR  LOZIER, DANIEL W. (NATIONAL BUREAU OF STANDARDS)
C           SMITH, JOHN M. (NBS AND GEORGE MASON UNIVERSITY)
C***PURPOSE  TO COMPUTE THE NORMALIZED LEGENDRE POLYNOMIAL
C***DESCRIPTION
C
C        SUBROUTINE TO CALCULATE NORMALIZED LEGENDRE POLYNOMIALS
C        (XDNRMP is double-precision version)
C        XSNRMP calculates normalized Legendre polynomials of varying
C        order and fixed argument and degree. The order MU and degree
C        NU are nonegative integers and the argument is real. Because
C        the algorithm requires the use of numbers outside the normal
C        machine range, this subroutine employs a special arithmetic
C        called extended-range arithmetic. See J.M. Smith, F.W.J. Olver,
C        and D.W. Lozier, Extended-Range Arithmetic and Normalized
C        Legendre Polynomials, ACM Transactions on Mathematical Soft-
C        ware, 93-105, March 1981, for a complete description of the
C        algorithm and special arithmetic. Also see program comments
C        in XSSET.
C
C        The normalized Legendre polynomials are multiples of the
C        associated Legendre polynomials of the first kind where the
C        normalizing coefficients are chosen so as to make the integral
C        from -1 to 1 of the square of each function equal to 1. See
C        E. Jahnke, F. Emde and F. Losch, Tables of Higher Functions,
C        McGraw-Hill, New York, 1960, p. 121.
C
C        The input values to XSNRMP are NU, MU1, MU2, SARG, and MODE.
C        These must satisfy
C          1. NU .GE. 0 specifies the degree of the normalized Legendre
C             polynomial that is wanted.
C          2. MU1 .GE. 0 specifies the lowest-order normalized Legendre
C             polynomial that is wanted.
C          3. MU2 .GE. MU1 specifies the highest-order normalized Leg-
C             endre polynomial that is wanted.
C         4a. MODE = 1 and -1.0 .LE. SARG .LE. 1.0 specifies that
C             Normalized Legendre(NU, MU, SARG) is wanted for MU = MU1,
C             MU1 + 1, ..., MU2.
C         4b. MODE = 2 and -3.14159... .LT. SARG .LT. 3.14159... spec-
C             ifies that Normalized Legendre(NU, MU, COS(SARG)) is want-
C             ed for MU = MU1, MU1 + 1, ..., MU2.
C
C        The output of XSNRMP consists of the two vectors SPN and IPN
C        and the error estimate ISIG. The computed values are stored as
C        extended-range numbers such that
C             (SPN(1),IPN(1))=NORMALIZED LEGENDRE(NU,MU1,X)
C             (SPN(2),IPN(2))=NORMALIZED LEGENDRE(NU,MU1+1,X)
C                .
C                .
C             (SPN(K),IPN(K))=NORMALIZED LEGENDRE(NU,MU2,X)
C        where K = MU2 - MU1 + 1 and X = SARG or COS(SARG) according
C        to whether MODE = 1 or 2. Finally, ISIG is an estimate of the
C        number of decimal digits lost through rounding errors in the
C        computation. For example if SARG is accurate to 12 significant
C        decimals, then the computed function values are accurate to
C        12 - ISIG significant decimals (except in neighborhoods of
C        zeros).
C
C        The interpretation of (SPN(I),IPN(I)) is SPN(I)*(IR**IPN(I))
C        where IR is the internal radix of the computer arithmetic. When
C        IPN(I) = 0 the value of the normalized Legendre polynomial is
C        contained entirely in SPN(I) and subsequent single-precision
C        computations can be performed without further consideration of
C        extended-range arithmetic. However, if IPN(I) .NE. 0 the corre-
C        sponding value of the normalized Legendre polynomial cannot be
C        represented in single-precision because of overflow or under-
C        flow. THE USER MUST TEST IPN(I) IN HIS/HER PROGRAM. In the even
C        that IPN(I) is nonzero, the user should try using double pre-
C        cision if it has a wider exponent range. If double precision
C        fails, the user could rewrite his/her program to use extended-
C        range arithmetic.
C
C        The interpretation of (SPN(I),IPN(I)) can be changed to
C        SPN(I)*(10**IPN(I)) by calling the extended-range subroutine
C        XSCON. This should be done before printing the computed values.
C        As an example of usage, the Fortran coding
C              J = K
C              DO 20 I = 1, K
C              CALL XSCON(SPN(I), IPN(I))
C              PRINT 10, SPN(I), IPN(I)
C           10 FORMAT(1X, E30.8 , I15)
C              IF ((IPN(I) .EQ. 0) .OR. (J .LT. K)) GO TO 20
C              J = I - 1
C           20 CONTINUE
C        will print all computed values and determine the largest J
C        such that IPN(1) = IPN(2) = ... = IPN(J) = 0. Because of the
C        change of representation caused by calling XSCON, (SPN(I),
C        IPN(I)) for I = J+1, J+2, ... cannot be used in subsequent
C        extended-range computations.
C
C***REFERENCES  (SEE DESCRIPTION ABOVE)
C***ROUTINES CALLED  XERROR, XSADD, XSADJ, XSRED, XSSET
C***END PROLOGUE  XSNRMP
      INTEGER NU, MU1, MU2, MODE, IPN, ISIG
      REAL SARG, SPN
      DIMENSION SPN(*), IPN(*)
      REAL C1,C2,P,P1,P2,P3,S,SX,T,TX,X,RK
C CALL XSSET TO INITIALIZE EXTENDED-RANGE ARITHMETIC (SEE XSSET
C LISTING FOR DETAILS)
C***FIRST EXECUTABLE STATEMENT  XSNRMP
      CALL XSSET (0, 0, 0.0, 0)
C
C        TEST FOR PROPER INPUT VALUES.
C
      IF (NU.LT.0) GO TO 110
      IF (MU1.LT.0) GO TO 110
      IF (MU1.GT.MU2) GO TO 110
      IF (NU.EQ.0) GO TO 90
      IF (MODE.LT.1 .OR. MODE.GT.2) GO TO 110
      GO TO (10, 20), MODE
   10 IF (ABS(SARG).GT.1.0) GO TO 120
      IF (ABS(SARG).EQ.1.0) GO TO 90
      X = SARG
      SX = SQRT((1.0+ABS(X))*((0.5-ABS(X))+0.5))
      TX = X/SX
      ISIG = ALOG10(2.0*(FLOAT(NU))*(5.0+TX**2))
      GO TO 30
   20 IF (ABS(SARG).GT.4.0*ATAN(1.0)) GO TO 120
      IF (SARG.EQ.0.0) GO TO 90
      X = COS(SARG)
      SX = ABS(SIN(SARG))
      TX = X/SX
      ISIG = ALOG10(2.0*FLOAT(NU)*(5.0+ABS(SARG*TX)))
C
C        BEGIN CALCULATION
C
   30 MU = MU2
      I = MU2 - MU1 + 1
C
C        IF MU.GT.NU, NORMALIZED LEGENDRE(NU,MU,X)=0.
C
   40 IF (MU.LE.NU) GO TO 50
      SPN(I) = 0.0
      IPN(I) = 0
      I = I - 1
      MU = MU - 1
      IF (I .GT. 0) GO TO 40
      ISIG = 0
      GO TO 160
   50 MU = NU
C
C        P1 = 0. = NORMALIZED LEGENDRE(NU,NU+1,X)
C
      P1 = 0.0
      IP1 = 0
C
C        CALCULATE P2 = NORMALIZED LEGENDRE(NU,NU,X)
C
      P2 = 1.0
      IP2 = 0
      P3 = 0.5
      RK = 2.0
      DO 60 J=1,NU
        P3 = ((RK+1.0)/RK)*P3
        P2 = P2*SX
        CALL XSADJ(P2, IP2)
        RK = RK + 2.0
   60 CONTINUE
      P2 = P2*SQRT(P3)
      CALL XSADJ(P2, IP2)
      S = 2.0*TX
      T = 1.0/FLOAT(NU)
      IF (MU2.LT.NU) GO TO 70
      SPN(I) = P2
      IPN(I) = IP2
      I = I - 1
      IF (I .EQ. 0) GO TO 140
C
C        RECURRENCE PROCESS
C
   70 P = FLOAT(MU)*T
      C1 = 1.0/SQRT((1.0-P+T)*(1.0+P))
      C2 = S*P*C1*P2
      C1 = -SQRT((1.0+P+T)*(1.0-P))*C1*P1
      CALL XSADD(C2, IP2, C1, IP1, P, IP)
      MU = MU - 1
      IF (MU.GT.MU2) GO TO 80
C
C        STORE IN ARRAY SPN FOR RETURN TO CALLING ROUTINE.
C
      SPN(I) = P
      IPN(I) = IP
      I = I - 1
      IF (I .EQ. 0) GO TO 140
   80 P1 = P2
      IP1 = IP2
      P2 = P
      IP2 = IP
      IF (MU.LE.MU1) GO TO 140
      GO TO 70
C
C        SPECIAL CASE WHEN X=-1 OR +1, OR NU=0.
C
   90 K = MU2 - MU1 + 1
      DO 100 I=1,K
        SPN(I) = 0.0
        IPN(I) = 0
  100 CONTINUE
      ISIG = 0
      IF (MU1.GT.0) GO TO 160
      ISIG = 1
      SPN(1) = SQRT(FLOAT(NU)+0.5)
      IPN(1) = 0
      IF (MOD(NU,2).EQ.0) GO TO 160
      IF (MODE.EQ.1 .AND. SARG.EQ.1.0) GO TO 160
      IF (MODE.EQ.2) GO TO 160
      SPN(1) = -SPN(1)
      GO TO 160
C
C          ERROR PRINTOUTS AND TERMINATION.
C
  110 CALL XERROR('Err in XSNRMP...NU,MU1,MU2 or MODE not valid',44,1,1)
      GO TO 130
  120 CALL XERROR('Err in XSNRMP...SARG out of range',33,2,1)
  130 RETURN
C
C        RETURN TO CALLING PROGRAM
C
  140 K = MU2 - MU1 + 1
      DO 150 I=1,K
        CALL XSRED(SPN(I),IPN(I))
  150 CONTINUE
  160 RETURN
      END
