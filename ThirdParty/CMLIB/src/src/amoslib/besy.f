      SUBROUTINE BESY(X,FNU,N,Y)
C***BEGIN PROLOGUE  BESY
C***DATE WRITTEN   800501   (YYMMDD)
C***REVISION DATE  870625   (YYMMDD)
C***CATEGORY NO.  C10A3
C***KEYWORDS  LIBRARY=SLATEC,TYPE=SINGLE PRECISION(BESY-S DBESY-D),
C             BESSEL FUNCTION,SPECIAL FUNCTIONS,Y BESSEL FUNCTION
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  BESY implements forward recursion on the three term recur-
C            sion relation for a sequence of non-negative order Bessel
C            functions Y/SUB(FNU+I-1)/(X), I=1,N for real X.GT.0.0E0 and
C            non-negative orders FNU.
C***DESCRIPTION
C
C     Written by D. E. Amos, May, 1980
C
C     References
C         SAND-80-1498
C
C         On the Numerical Evaluation of the Ordinary Bessel Function
C         of the Second Kind by N. M. Temme, J. Comp. Physics, 21,
C         1976, pp. 343-350.
C
C         On the Numerical Evaluation of the Modified Bessel Function
C         of the Third Kind by N. M. Temme, J. Comp. Physics, 19, 1975,
C         pp. 324-337.
C
C         Tables of Bessel Functions of Moderate or Large Orders,
C         NPL Mathematical Tables, Vol. 6, by F. W. J. Olver, Her
C         Majesty's Stationery Office, London, 1962.
C
C     Abstract
C         BESY implements forward recursion on the three term
C         recursion relation for a sequence of non-negative order Bessel
C         functions Y/sub(FNU+I-1)/(X), I=1,N for real X .GT. 0.0E0 and
C         non-negative orders FNU.  If FNU .LT. NULIM, orders FNU and
C         FNU+1 are obtained from BESYNU which computes by a power
C         series for X .LE. 2, the K Bessel function of an imaginary
C         argument for 2 .LT. X .LE. 20 and the asymptotic expansion for
C         X .GT. 20.
C
C         If FNU .GE. NULIM, the uniform asymptotic expansion is coded
C         in ASYJY for orders FNU and FNU+1 to start the recursion.
C         NULIM is 70 or 100 depending on whether N=1 or N .GE. 2.  An
C         overflow test is made on the leading term of the asymptotic
C         expansion before any extensive computation is done.
C
C         BESY calls BESYNU,GAMMA,ASYJY,SIH,COSH,BESY0,BESY1
C                    I1MACH,R1MACH,XERROR
C
C     Description of Arguments
C
C         Input
C           X      - X .GT. 0.0E0
C           FNU    - order of the initial Y function, FNU .GE. 0.0E0
C           N      - number of members in the sequence, N .GE. 1
C
C         Output
C           Y      - a vector whose first N components contain values
C                    for the sequence Y(I)=Y/sub(FNU+I-1)/(X), I=1,N.
C
C     Error Conditions
C         Improper input arguments - a fatal error
C         Overflow - a fatal error
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ASYJY,BESY0,BESY1,BESYNU,I1MACH,R1MACH,XERROR,
C                    YAIRY
C***END PROLOGUE  BESY
C
      EXTERNAL YAIRY
      INTEGER I, IFLW, J, N, NB, ND, NN, NUD, NULIM
      INTEGER I1MACH
      REAL       AZN,CN,DNU,ELIM,FLGJY,FN,FNU,RAN,S,S1,S2,TM,TRX,
     1           W,WK,W2N,X,XLIM,XXN,Y
      REAL BESY0, BESY1, R1MACH
      DIMENSION W(2), NULIM(2), Y(*), WK(7)
      SAVE NULIM
      DATA NULIM(1),NULIM(2) / 70 , 100 /
C***FIRST EXECUTABLE STATEMENT  BESY
      NN = -I1MACH(12)
      ELIM = 2.303E0*(FLOAT(NN)*R1MACH(5)-3.0E0)
      XLIM = R1MACH(1)*1.0E+3
      IF (FNU.LT.0.0E0) GO TO 140
      IF (X.LE.0.0E0) GO TO 150
      IF (X.LT.XLIM) GO TO 170
      IF (N.LT.1) GO TO 160
C
C     ND IS A DUMMY VARIABLE FOR N
C
      ND = N
      NUD = INT(FNU)
      DNU = FNU - FLOAT(NUD)
      NN = MIN0(2,ND)
      FN = FNU + FLOAT(N-1)
      IF (FN.LT.2.0E0) GO TO 100
C
C     OVERFLOW TEST  (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
C     FOR THE LAST ORDER, FNU+N-1.GE.NULIM
C
      XXN = X/FN
      W2N = 1.0E0-XXN*XXN
      IF(W2N.LE.0.0E0) GO TO 10
      RAN = SQRT(W2N)
      AZN = ALOG((1.0E0+RAN)/XXN) - RAN
      CN = FN*AZN
      IF(CN.GT.ELIM) GO TO 170
   10 CONTINUE
      IF (NUD.LT.NULIM(NN)) GO TO 20
C
C     ASYMPTOTIC EXPANSION FOR ORDERS FNU AND FNU+1.GE.NULIM
C
      FLGJY = -1.0E0
      CALL ASYJY(YAIRY,X,FNU,FLGJY,NN,Y,WK,IFLW)
      IF(IFLW.NE.0) GO TO 170
      IF (NN.EQ.1) RETURN
      TRX = 2.0E0/X
      TM = (FNU+FNU+2.0E0)/X
      GO TO 80
C
   20 CONTINUE
      IF (DNU.NE.0.0E0) GO TO 30
      S1 = BESY0(X)
      IF (NUD.EQ.0 .AND. ND.EQ.1) GO TO 70
      S2 = BESY1(X)
      GO TO 40
   30 CONTINUE
      NB = 2
      IF (NUD.EQ.0 .AND. ND.EQ.1) NB = 1
      CALL BESYNU(X, DNU, NB, W)
      S1 = W(1)
      IF (NB.EQ.1) GO TO 70
      S2 = W(2)
   40 CONTINUE
      TRX = 2.0E0/X
      TM = (DNU+DNU+2.0E0)/X
C     FORWARD RECUR FROM DNU TO FNU+1 TO GET Y(1) AND Y(2)
      IF (ND.EQ.1) NUD = NUD - 1
      IF (NUD.GT.0) GO TO 50
      IF (ND.GT.1) GO TO 70
      S1 = S2
      GO TO 70
   50 CONTINUE
      DO 60 I=1,NUD
        S = S2
        S2 = TM*S2 - S1
        S1 = S
        TM = TM + TRX
   60 CONTINUE
      IF (ND.EQ.1) S1 = S2
   70 CONTINUE
      Y(1) = S1
      IF (ND.EQ.1) RETURN
      Y(2) = S2
   80 CONTINUE
      IF (ND.EQ.2) RETURN
C     FORWARD RECUR FROM FNU+2 TO FNU+N-1
      DO 90 I=3,ND
        Y(I) = TM*Y(I-1) - Y(I-2)
        TM = TM + TRX
   90 CONTINUE
      RETURN
C
  100 CONTINUE
C     OVERFLOW TEST
      IF (FN.LE.1.0E0) GO TO 110
      IF (-FN*(ALOG(X)-0.693E0).GT.ELIM) GO TO 170
  110 CONTINUE
      IF (DNU.EQ.0.0E0) GO TO 120
      CALL BESYNU(X, FNU, ND, Y)
      RETURN
  120 CONTINUE
      J = NUD
      IF (J.EQ.1) GO TO 130
      J = J + 1
      Y(J) = BESY0(X)
      IF (ND.EQ.1) RETURN
      J = J + 1
  130 CONTINUE
      Y(J) = BESY1(X)
      IF (ND.EQ.1) RETURN
      TRX = 2.0E0/X
      TM = TRX
      GO TO 80
C
C
C
  140 CONTINUE
      CALL XERROR( 'IN BESY, ORDER, FNU, LESS THAN ZERO', 35, 2, 1)
      RETURN
  150 CONTINUE
      CALL XERROR( 'IN BESY, X LESS THAN OR EQUAL TO ZERO', 37, 2, 1)
      RETURN
  160 CONTINUE
      CALL XERROR( 'IN BESY, N LESS THAN ONE', 24, 2, 1)
      RETURN
  170 CONTINUE
      CALL XERROR( 'IN BESY, OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL
     1', 52, 6,   1)
      RETURN
      END
