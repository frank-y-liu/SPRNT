      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z,IERR,RV1,RV2,RV3,RV4,RV6)
C***BEGIN PROLOGUE  TINVIT
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C3
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Eigenvectors of symmetric tridiagonal matrix corresponding
C            to some specified eigenvalues, using inverse iteration
C***DESCRIPTION
C
C     This subroutine is a translation of the inverse iteration tech-
C     nique in the ALGOL procedure TRISTURM by Peters and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     This subroutine finds those eigenvectors of a TRIDIAGONAL
C     SYMMETRIC matrix corresponding to specified eigenvalues,
C     using inverse iteration.
C
C     On Input
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        D contains the diagonal elements of the input matrix.
C
C        E contains the subdiagonal elements of the input matrix
C          in its last N-1 positions.  E(1) is arbitrary.
C
C        E2 contains the squares of the corresponding elements of E,
C          with zeros corresponding to negligible elements of E.
C          E(I) is considered negligible if it is not larger than
C          the product of the relative machine precision and the sum
C          of the magnitudes of D(I) and D(I-1).  E2(1) must contain
C          0.0e0 if the eigenvalues are in ascending order, or 2.0e0
C          if the eigenvalues are in descending order.  If  BISECT,
C          TRIDIB, or  IMTQLV  has been used to find the eigenvalues,
C          their output E2 array is exactly what is expected here.
C
C        M is the number of specified eigenvalues.
C
C        W CONTAINS the M eigenvalues in ascending or descending order.
C
C        IND contains in its first M positions the submatrix indices
C          associated with the corresponding eigenvalues in W --
C          1 for eigenvalues belonging to the first submatrix from
C          the top, 2 for those belonging to the second submatrix, etc.
C
C     On Output
C
C       ** All input arrays are unaltered.**
C
C        Z contains the associated set of orthonormal eigenvectors.
C          any vector which fails to converge is set to zero.
C
C        IERR is set to
C          Zero       for normal return,
C          -R         if the eigenvector corresponding to the R-th
C                     eigenvalue fails to converge in 5 iterations.
C
C        RV1, RV2, RV3, RV4, and RV6 are temporary storage arrays.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  TINVIT
C
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP
      INTEGER IND(M)
      REAL D(N),E(N),E2(N),W(M),Z(NM,M)
      REAL RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      REAL U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER
C
C***FIRST EXECUTABLE STATEMENT  TINVIT
      IERR = 0
      IF (M .EQ. 0) GO TO 1001
      TAG = 0
      ORDER = 1.0E0 - E2(1)
      Q = 0
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX ..........
  100 P = Q + 1
C
      DO 120 Q = P, N
         IF (Q .EQ. N) GO TO 140
         IF (E2(Q+1) .EQ. 0.0E0) GO TO 140
  120 CONTINUE
C     .......... FIND VECTORS BY INVERSE ITERATION ..........
  140 TAG = TAG + 1
      S = 0
C
      DO 920 R = 1, M
         IF (IND(R) .NE. TAG) GO TO 920
         ITS = 1
         X1 = W(R)
         IF (S .NE. 0) GO TO 510
C     .......... CHECK FOR ISOLATED ROOT ..........
         XU = 1.0E0
         IF (P .NE. Q) GO TO 490
         RV6(P) = 1.0E0
         GO TO 870
  490    NORM = ABS(D(P))
         IP = P + 1
C
         DO 500 I = IP, Q
  500    NORM = AMAX1(NORM, ABS(D(I)) + ABS(E(I)))
C     .......... EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
         EPS2 = 1.0E-3 * NORM
         EPS3 = NORM
  502    EPS3 = 0.5E0*EPS3
         IF (NORM + EPS3 .GT. NORM) GO TO 502
         UK = SQRT(FLOAT(Q-P+5))
         EPS3 = UK * EPS3
         EPS4 = UK * EPS3
         UK = EPS4 / UK
         S = P
  505    GROUP = 0
         GO TO 520
C     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
  510    IF (ABS(X1-X0) .GE. EPS2) GO TO 505
         GROUP = GROUP + 1
         IF (ORDER * (X1 - X0) .LE. 0.0E0) X1 = X0 + ORDER * EPS3
C     .......... ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR ..........
  520    V = 0.0E0
C
         DO 580 I = P, Q
            RV6(I) = UK
            IF (I .EQ. P) GO TO 560
            IF (ABS(E(I)) .LT. ABS(U)) GO TO 540
C     .......... WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY ..........
            XU = U / E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = 0.0E0
            IF (I .NE. Q) RV3(I-1) = E(I+1)
            U = V - XU * RV2(I-1)
            V = -XU * RV3(I-1)
            GO TO 580
  540       XU = E(I) / U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = 0.0E0
  560       U = D(I) - X1 - XU * V
            IF (I .NE. Q) V = E(I+1)
  580    CONTINUE
C
         IF (U .EQ. 0.0E0) U = EPS3
         RV1(Q) = U
         RV2(Q) = 0.0E0
         RV3(Q) = 0.0E0
C     .......... BACK SUBSTITUTION
C                FOR I=Q STEP -1 UNTIL P DO -- ..........
  600    DO 620 II = P, Q
            I = P + Q - II
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
            V = U
            U = RV6(I)
  620    CONTINUE
C     .......... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP ..........
         IF (GROUP .EQ. 0) GO TO 700
         J = R
C
         DO 680 JJ = 1, GROUP
  630       J = J - 1
            IF (IND(J) .NE. TAG) GO TO 630
            XU = 0.0E0
C
            DO 640 I = P, Q
  640       XU = XU + RV6(I) * Z(I,J)
C
            DO 660 I = P, Q
  660       RV6(I) = RV6(I) - XU * Z(I,J)
C
  680    CONTINUE
C
  700    NORM = 0.0E0
C
         DO 720 I = P, Q
  720    NORM = NORM + ABS(RV6(I))
C
         IF (NORM .GE. 1.0E0) GO TO 840
C     .......... FORWARD SUBSTITUTION ..........
         IF (ITS .EQ. 5) GO TO 830
         IF (NORM .NE. 0.0E0) GO TO 740
         RV6(S) = EPS4
         S = S + 1
         IF (S .GT. Q) S = P
         GO TO 780
  740    XU = EPS4 / NORM
C
         DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
C     .......... ELIMINATION OPERATIONS ON NEXT VECTOR
C                ITERATE ..........
  780    DO 820 I = IP, Q
            U = RV6(I)
C     .......... IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                WAS PERFORMED EARLIER IN THE
C                TRIANGULARIZATION PROCESS ..........
            IF (RV1(I-1) .NE. E(I)) GO TO 800
            U = RV6(I-1)
            RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
C
         ITS = ITS + 1
         GO TO 600
C     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
  830    IERR = -R
         XU = 0.0E0
         GO TO 870
C     .......... NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER ..........
  840    U = 0.0E0
C
         DO 860 I = P, Q
  860    U = U + RV6(I)**2
C
         XU = 1.0E0 / SQRT(U)
C
  870    DO 880 I = 1, N
  880    Z(I,R) = 0.0E0
C
         DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
C
         X0 = X1
  920 CONTINUE
C
      IF (Q .LT. N) GO TO 100
 1001 RETURN
      END
