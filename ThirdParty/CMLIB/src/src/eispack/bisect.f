      SUBROUTINE BISECT(N,EPS1,D,E,E2,LB,UB,MM,M,W,IND,IERR,RV4,RV5)
C***BEGIN PROLOGUE  BISECT
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4A5,D4C2A
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Compute eigenvalues of symmetric tridiagonal matrix
C            given interval using Sturm sequencing.
C***DESCRIPTION
C
C     This subroutine is a translation of the bisection technique
C     in the ALGOL procedure TRISTURM by Peters and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     This subroutine finds those eigenvalues of a TRIDIAGONAL
C     SYMMETRIC matrix which lie in a specified interval,
C     using bisection.
C
C     On INPUT
C
C        N is the order of the matrix.
C
C        EPS1 is an absolute error tolerance for the computed
C          eigenvalues.  If the input EPS1 is non-positive,
C          it is reset for each submatrix to a default value,
C          namely, minus the product of the relative machine
C          precision and the 1-norm of the submatrix.
C
C        D contains the diagonal elements of the input matrix.
C
C        E contains the subdiagonal elements of the input matrix
C          in its last N-1 positions.  E(1) is arbitrary.
C
C        E2 contains the squares of the corresponding elements of E.
C          E2(1) is arbitrary.
C
C        LB and UB define the interval to be searched for eigenvalues.
C          If LB is not less than UB, no eigenvalues will be found.
C
C        MM should be set to an upper bound for the number of
C          eigenvalues in the interval.  WARNING. If more than
C          MM eigenvalues are determined to lie in the interval,
C          an error return is made with no eigenvalues found.
C
C     On OUTPUT
C
C        EPS1 is unaltered unless it has been reset to its
C          (last) default value.
C
C        D and E are unaltered.
C
C        Elements of E2, corresponding to elements of E regarded
C          as negligible, have been replaced by zero causing the
C          matrix to split into a direct sum of submatrices.
C          E2(1) is also set to zero.
C
C        M is the number of eigenvalues determined to lie in (LB,UB).
C
C        W contains the M eigenvalues in ascending order.
C
C        IND contains in its first M positions the submatrix indices
C          associated with the corresponding eigenvalues in W --
C          1 for eigenvalues belonging to the first submatrix from
C          the top, 2 for those belonging to the second submatrix, etc.
C
C        IERR is set to
C          Zero       for normal return,
C          3*N+1      if M exceeds MM.
C
C        RV4 and RV5 are temporary storage arrays.
C
C     The ALGOL procedure STURMCNT contained in TRISTURM
C     appears in BISECT in-line.
C
C     Note that subroutine TQL1 or IMTQL1 is generally faster than
C     BISECT, if more than N/4 eigenvalues are to be found.
C
C     Questions and comments should be directed to B. S. Garbow,
C     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  BISECT
C
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,MM,M1,M2,TAG,IERR,ISTURM
      REAL D(N),E(N),E2(N),W(MM),RV4(N),RV5(N)
      REAL U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,MACHEP,S1,S2
      INTEGER IND(MM)
C
      DATA MACHEP/1.0E0/
C***FIRST EXECUTABLE STATEMENT  BISECT
      IF (MACHEP .NE. 1.0E0) GO TO 10
C
C   --- This code fails to compute MACHEP correctly on IBM machines. ---
C   --- Replaced by call to R1MACH on 15 Jun 94 by Ron Boisvert.     ---
C
C   05 MACHEP = 0.5E0*MACHEP
C      IF (1.0E0 + MACHEP .GT. 1.0E0) GO TO 05
C      MACHEP = 2.0E0*MACHEP
C
      MACHEP = R1MACH(4)
C
   10 IERR = 0
      TAG = 0
      T1 = LB
      T2 = UB
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES ..........
      DO 40 I = 1, N
         IF (I .EQ. 1) GO TO 20
         S1 = ABS(D(I)) + ABS(D(I-1))
         S2 = S1 + ABS(E(I))
         IF (S2 .GT. S1) GO TO 40
   20    E2(I) = 0.0E0
   40 CONTINUE
C     .......... DETERMINE THE NUMBER OF EIGENVALUES
C                IN THE INTERVAL ..........
      P = 1
      Q = N
      X1 = UB
      ISTURM = 1
      GO TO 320
   60 M = S
      X1 = LB
      ISTURM = 2
      GO TO 320
   80 M = M - S
      IF (M .GT. MM) GO TO 980
      Q = 0
      R = 0
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
C                INTERVAL BY THE GERSCHGORIN BOUNDS ..........
  100 IF (R .EQ. M) GO TO 1001
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0E0
C
      DO 120 Q = P, N
         X1 = U
         U = 0.0E0
         V = 0.0E0
         IF (Q .EQ. N) GO TO 110
         U = ABS(E(Q+1))
         V = E2(Q+1)
  110    XU = AMIN1(D(Q)-(X1+U),XU)
         X0 = AMAX1(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0E0) GO TO 140
  120 CONTINUE
C
  140 X1 = AMAX1(ABS(XU),ABS(X0)) * MACHEP
      IF (EPS1 .LE. 0.0E0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
C     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * FLOAT(Q-P+1)
      LB = AMAX1(T1,XU-X1)
      UB = AMIN1(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
C     .......... FIND ROOTS BY BISECTION ..........
      X0 = UB
      ISTURM = 5
C
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
C     .......... LOOP FOR K-TH EIGENVALUE
C                FOR K=M2 STEP -1 UNTIL M1 DO --
C                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) ..........
      K = M2
  250    XU = LB
C     .......... FOR I=K STEP -1 UNTIL M1 DO -- ..........
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
C
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
C     .......... NEXT BISECTION STEP ..........
  300    X1 = (XU + X0) * 0.5E0
         S1 = 2.0E0*(ABS(XU) + ABS(X0) + ABS(EPS1))
         S2 = S1 + ABS(X0 - XU)
         IF (S2 .EQ. S1) GO TO 420
C     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
  320    S = P - 1
         U = 1.0E0
C
         DO 340 I = P, Q
            IF (U .NE. 0.0E0) GO TO 325
            V = ABS(E(I)) / MACHEP
            IF (E2(I) .EQ. 0.0E0) V = 0.0E0
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0E0) S = S + 1
  340    CONTINUE
C
         GO TO (60,80,200,220,360), ISTURM
C     .......... REFINE INTERVALS ..........
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
C     .......... K-TH EIGENVALUE FOUND ..........
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
C     .......... ORDER EIGENVALUES TAGGED WITH THEIR
C                SUBMATRIX ASSOCIATIONS ..........
  900 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
C
      DO 920 L = 1, R
         IF (J .GT. S) GO TO 910
         IF (K .GT. M2) GO TO 940
         IF (RV5(K) .GE. W(L)) GO TO 915
C
         DO 905 II = J, S
            I = L + S - II
            W(I+1) = W(I)
            IND(I+1) = IND(I)
  905    CONTINUE
C
  910    W(L) = RV5(K)
         IND(L) = TAG
         K = K + 1
         GO TO 920
  915    J = J + 1
  920 CONTINUE
C
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
C     .......... SET ERROR -- UNDERESTIMATE OF NUMBER OF
C                EIGENVALUES IN INTERVAL ..........
  980 IERR = 3 * N + 1
 1001 LB = T1
      UB = T2
      RETURN
      END
