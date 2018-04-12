      SUBROUTINE RATQR(N,EPS1,D,E,E2,M,W,IND,BD,TYPE,IDEF,IERR)
C***BEGIN PROLOGUE  RATQR
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4A5,D4C2A
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Computes largest or smallest eigenvalues of symmetric
C            tridiagonal matrix using rational QR method with Newton
C            correction.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure RATQR,
C     NUM. MATH. 11, 264-272(1968) by REINSCH and BAUER.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 257-265(1971).
C
C     This subroutine finds the algebraically smallest or largest
C     eigenvalues of a SYMMETRIC TRIDIAGONAL matrix by the
C     rational QR method with Newton corrections.
C
C     On Input
C
C        N is the order of the matrix.
C
C        EPS1 is a theoretical absolute error tolerance for the
C          computed eigenvalues.  If the input EPS1 is non-positive,
C          or indeed smaller than its default value, it is reset
C          at each iteration to the respective default value,
C          namely, the product of the relative machine precision
C          and the magnitude of the current eigenvalue iterate.
C          The theoretical absolute error in the K-th eigenvalue
C          is usually not greater than K times EPS1.
C
C        D contains the diagonal elements of the input matrix.
C
C        E contains the subdiagonal elements of the input matrix
C          in its last N-1 positions.  E(1) is arbitrary.
C
C        E2 contains the squares of the corresponding elements of E.
C          E2(1) is arbitrary.
C
C        M is the number of eigenvalues to be found.
C
C        IDEF should be set to 1 if the input matrix is known to be
C          positive definite, to -1 if the input matrix is known to
C          be negative definite, and to 0 otherwise.
C
C        TYPE should be set to .TRUE. if the smallest eigenvalues
C          are to be found, and to .FALSE. If the largest eigenvalues
C          are to be found.
C
C     On Output
C
C        EPS1 is unaltered unless it has been reset to its
C          (last) default value.
C
C        D and E are unaltered (unless W overwrites D).
C
C        ELEMENTS of E2, corresponding to elements of E regarded
C          as negligible, have been replaced by zero causing the
C          matrix to split into a direct sum of submatrices.
C          E2(1) is set to 0.0e0 if the smallest eigenvalues have been
C          found, and to 2.0e0 if the largest eigenvalues have been
C          found.  E2 is otherwise unaltered (unless overwritten by BD).
C
C        W contains the M algebraically smallest eigenvalues in
C          ascending order, or the M largest eigenvalues in
C          descending order.  If an error exit is made because of
C          an incorrect specification of IDEF, no eigenvalues
C          are found.  If the Newton iterates for a particular
C          eigenvalue are not monotone, the best estimate obtained
C          is returned and IERR is set.  W may coincide with D.
C
C        IND contains in its first M positions the submatrix indices
C          associated with the corresponding eigenvalues in W --
C          1 for eigenvalues belonging to the first submatrix from
C          the top, 2 for those belonging to the second submatrix, etc.
C
C        BD contains refined bounds for the theoretical errors of the
C          corresponding eigenvalues in W.  These bounds are usually
C          within the tolerance specified by EPS1.  BD may coincide
C          with E2.
C
C        IERR is set to
C          Zero       for normal return,
C          6*N+1      if  IDEF  is set to 1 and  type  to .TRUE.
C                     when the matrix is NOT positive definite, or
C                     if  IDEF  is set to -1 and  type  to .FALSE.
C                     when the matrix is NOT negative definite,
C          5*N+K      if successive iterates to the K-th eigenvalue
C                     are NOT monotone increasing, where K refers
C                     to the last such occurrence.
C
C     Note that subroutine TRIDIB is generally faster and more
C     accurate than RATQR if the eigenvalues are clustered.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  RATQR
C
      INTEGER I,J,K,M,N,II,JJ,K1,IDEF,IERR,JDEF
      REAL D(N),E(N),E2(N),W(N),BD(N)
      REAL F,P,Q,R,S,EP,QP,ERR,TOT,EPS1,DELTA,MACHEP
      INTEGER IND(N)
      LOGICAL TYPE
C
      DATA MACHEP/1.0E0/
C***FIRST EXECUTABLE STATEMENT  RATQR
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
      JDEF = IDEF
C     .......... COPY D ARRAY INTO W ..........
      DO 20 I = 1, N
   20 W(I) = D(I)
C
      IF (TYPE) GO TO 40
      J = 1
      GO TO 400
   40 ERR = 0.0E0
      S = 0.0E0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DEFINE
C                INITIAL SHIFT FROM LOWER GERSCHGORIN BOUND.
C                COPY E2 ARRAY INTO BD ..........
      TOT = W(1)
      Q = 0.0E0
      J = 0
C
      DO 100 I = 1, N
         P = Q
         IF (I .EQ. 1) GO TO 60
         IF (P .GT. MACHEP * (ABS(D(I)) + ABS(D(I-1)))) GO TO 80
   60    E2(I) = 0.0E0
   80    BD(I) = E2(I)
C     .......... COUNT ALSO IF ELEMENT OF E2 HAS UNDERFLOWED ..........
         IF (E2(I) .EQ. 0.0E0) J = J + 1
         IND(I) = J
         Q = 0.0E0
         IF (I .NE. N) Q = ABS(E(I+1))
         TOT = AMIN1(W(I)-P-Q,TOT)
  100 CONTINUE
C
      IF (JDEF .EQ. 1 .AND. TOT .LT. 0.0E0) GO TO 140
C
      DO 110 I = 1, N
  110 W(I) = W(I) - TOT
C
      GO TO 160
  140 TOT = 0.0E0
C
  160 DO 360 K = 1, M
C     .......... NEXT QR TRANSFORMATION ..........
  180    TOT = TOT + S
         DELTA = W(N) - S
         I = N
         F = ABS(MACHEP*TOT)
         IF (EPS1 .LT. F) EPS1 = F
         IF (DELTA .GT. EPS1) GO TO 190
         IF (DELTA .LT. (-EPS1)) GO TO 1000
         GO TO 300
C     .......... REPLACE SMALL SUB-DIAGONAL SQUARES BY ZERO
C                TO REDUCE THE INCIDENCE OF UNDERFLOWS ..........
  190    IF (K .EQ. N) GO TO 210
         K1 = K + 1
         DO 200 J = K1, N
            IF (BD(J) .LE. (MACHEP*(W(J)+W(J-1))) ** 2) BD(J) = 0.0E0
  200    CONTINUE
C
  210    F = BD(N) / DELTA
         QP = DELTA + F
         P = 1.0E0
         IF (K .EQ. N) GO TO 260
         K1 = N - K
C     .......... FOR I=N-1 STEP -1 UNTIL K DO -- ..........
         DO 240 II = 1, K1
            I = N - II
            Q = W(I) - S - F
            R = Q / QP
            P = P * R + 1.0E0
            EP = F * R
            W(I+1) = QP + EP
            DELTA = Q - EP
            IF (DELTA .GT. EPS1) GO TO 220
            IF (DELTA .LT. (-EPS1)) GO TO 1000
            GO TO 300
  220       F = BD(I) / Q
            QP = DELTA + F
            BD(I+1) = QP * EP
  240    CONTINUE
C
  260    W(K) = QP
         S = QP / P
         IF (TOT + S .GT. TOT) GO TO 180
C     .......... SET ERROR -- IRREGULAR END OF ITERATION.
C                DEFLATE MINIMUM DIAGONAL ELEMENT ..........
         IERR = 5 * N + K
         S = 0.0E0
         DELTA = QP
C
         DO 280 J = K, N
            IF (W(J) .GT. DELTA) GO TO 280
            I = J
            DELTA = W(J)
  280    CONTINUE
C     .......... CONVERGENCE ..........
  300    IF (I .LT. N) BD(I+1) = BD(I) * F / QP
         II = IND(I)
         IF (I .EQ. K) GO TO 340
         K1 = I - K
C     .......... FOR J=I-1 STEP -1 UNTIL K DO -- ..........
         DO 320 JJ = 1, K1
            J = I - JJ
            W(J+1) = W(J) - S
            BD(J+1) = BD(J)
            IND(J+1) = IND(J)
  320    CONTINUE
C
  340    W(K) = TOT
         ERR = ERR + ABS(DELTA)
         BD(K) = ERR
         IND(K) = II
  360 CONTINUE
C
      IF (TYPE) GO TO 1001
      F = BD(1)
      E2(1) = 2.0E0
      BD(1) = F
      J = 2
C     .......... NEGATE ELEMENTS OF W FOR LARGEST VALUES ..........
  400 DO 500 I = 1, N
  500 W(I) = -W(I)
C
      JDEF = -JDEF
      GO TO (40,1001), J
C     .......... SET ERROR -- IDEF SPECIFIED INCORRECTLY ..........
 1000 IERR = 6 * N + 1
 1001 RETURN
      END
