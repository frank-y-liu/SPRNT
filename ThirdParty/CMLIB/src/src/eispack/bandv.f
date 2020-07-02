      SUBROUTINE BANDV(NM,N,MBW,A,E21,M,W,Z,IERR,NV,RV,RV6)
C***BEGIN PROLOGUE  BANDV
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C3
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Forms eigenvectors of real symmetric band matrix
C            associated with a set of ordered approximate eigenvalue
C            by inverse iteration.
C***DESCRIPTION
C
C     This subroutine finds those eigenvectors of a REAL SYMMETRIC
C     BAND matrix corresponding to specified eigenvalues, using inverse
C     iteration.  The subroutine may also be used to solve systems
C     of linear equations with a symmetric or non-symmetric band
C     coefficient matrix.
C
C     On INPUT
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        MBW is the number of columns of the array A used to store the
C          band matrix.  If the matrix is symmetric, MBW is its (half)
C          band width, denoted MB and defined as the number of adjacent
C          diagonals, including the principal diagonal, required to
C          specify the non-zero portion of the lower triangle of the
C          matrix.  If the subroutine is being used to solve systems
C          of linear equations and the coefficient matrix is not
C          symmetric, it must however have the same number of adjacent
C          diagonals above the main diagonal as below, and in this
C          case, MBW=2*MB-1.
C
C        A contains the lower triangle of the symmetric band input
C          matrix stored as an N by MB array.  Its lowest subdiagonal
C          is stored in the last N+1-MB positions of the first column,
C          its next subdiagonal in the last N+2-MB positions of the
C          second column, further subdiagonals similarly, and finally
C          its principal diagonal in the N positions of column MB.
C          If the subroutine is being used to solve systems of linear
C          equations and the coefficient matrix is not symmetric, A is
C          N by 2*MB-1 instead with lower triangle as above and with
C          its first superdiagonal stored in the first N-1 positions of
C          column MB+1, its second superdiagonal in the first N-2
C          positions of column MB+2, further superdiagonals similarly,
C          and finally its highest superdiagonal in the first N+1-MB
C          positions of the last column.
C          Contents of storages not part of the matrix are arbitrary.
C
C        E21 specifies the ordering of the eigenvalues and contains
C            0.0E0 if the eigenvalues are in ascending order, or
C            2.0E0 if the eigenvalues are in descending order.
C          If the subroutine is being used to solve systems of linear
C          equations, E21 should be set to 1.0E0 if the coefficient
C          matrix is symmetric and to -1.0E0 if not.
C
C        M is the number of specified eigenvalues or the number of
C          systems of linear equations.
C
C        W contains the M eigenvalues in ascending or descending order.
C          If the subroutine is being used to solve systems of linear
C          equations (A-W(R)*I)*X(R)=B(R), where I is the identity
C          matrix, W(R) should be set accordingly, for R=1,2,...,M.
C
C        Z contains the constant matrix columns (B(R),R=1,2,...,M), if
C          the subroutine is used to solve systems of linear equations.
C
C        NV must be set to the dimension of the array parameter RV
C          as declared in the calling program dimension statement.
C
C     On OUTPUT
C
C        A and W are unaltered.
C
C        Z contains the associated set of orthogonal eigenvectors.
C          Any vector which fails to converge is set to zero.  If the
C          subroutine is used to solve systems of linear equations,
C          Z contains the solution matrix columns (X(R),R=1,2,...,M).
C
C        IERR is set to
C          Zero       for normal return,
C          -R         if the eigenvector corresponding to the R-th
C                     eigenvalue fails to converge, or if the R-th
C                     system of linear equations is nearly singular.
C
C        RV and RV6 are temporary storage arrays.  Note that RV is
C          of dimension at least N*(2*MB-1).  If the subroutine
C          is being used to solve systems of linear equations, the
C          determinant (up to sign) of A-W(M)*I is available, upon
C          return, as the product of the first N elements of RV.
C
C     Questions and comments should be directed to B. S. Garbow,
C     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  BANDV
C
      INTEGER I,J,K,M,N,R,II,IJ,JJ,KJ,MB,M1,NM,NV,IJ1,ITS,KJ1,MBW,M21
      INTEGER IERR,MAXJ,MAXK,GROUP
      REAL A(NM,MBW),W(M),Z(NM,M),RV(NV),RV6(N)
      REAL U,V,UK,XU,X0,X1,E21,EPS2,EPS3,EPS4,NORM,ORDER,S
C
C***FIRST EXECUTABLE STATEMENT  BANDV
      IERR = 0
      IF (M .EQ. 0) GO TO 1001
      MB = MBW
      IF (E21 .LT. 0.0E0) MB = (MBW + 1) / 2
      M1 = MB - 1
      M21 = M1 + MB
      ORDER = 1.0E0 - ABS(E21)
C     .......... FIND VECTORS BY INVERSE ITERATION ..........
      DO 920 R = 1, M
         ITS = 1
         X1 = W(R)
         IF (R .NE. 1) GO TO 100
C     .......... COMPUTE NORM OF MATRIX ..........
         NORM = 0.0E0
C
         DO 60 J = 1, MB
            JJ = MB + 1 - J
            KJ = JJ + M1
            IJ = 1
            S = 0.0E0
C
            DO 40 I = JJ, N
               S = S + ABS(A(I,J))
               IF (E21 .GE. 0.0E0) GO TO 40
               S = S + ABS(A(IJ,KJ))
               IJ = IJ + 1
   40       CONTINUE
C
            NORM = AMAX1(NORM,S)
   60    CONTINUE
C
         IF (E21 .LT. 0.0E0) NORM = 0.5E0 * NORM
C     .......... EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
         IF (NORM .EQ. 0.0E0) NORM = 1.0E0
         EPS2 = 1.0E-3 * NORM * ABS(ORDER)
         EPS3 = NORM
   70    EPS3 = 0.5E0*EPS3
         IF (NORM + EPS3 .GT. NORM) GO TO 70
         UK = SQRT(FLOAT(N))
         EPS3 = UK * EPS3
         EPS4 = UK * EPS3
   80    GROUP = 0
         GO TO 120
C     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
  100    IF (ABS(X1-X0) .GE. EPS2) GO TO 80
         GROUP = GROUP + 1
         IF (ORDER * (X1 - X0) .LE. 0.0E0) X1 = X0 + ORDER * EPS3
C     .......... EXPAND MATRIX, SUBTRACT EIGENVALUE,
C                AND INITIALIZE VECTOR ..........
  120    DO 200 I = 1, N
            IJ = I + MIN0(0,I-M1) * N
            KJ = IJ + MB * N
            IJ1 = KJ + M1 * N
            IF (M1 .EQ. 0) GO TO 180
C
            DO 150 J = 1, M1
               IF (IJ .GT. M1) GO TO 125
               IF (IJ .GT. 0) GO TO 130
               RV(IJ1) = 0.0E0
               IJ1 = IJ1 + N
               GO TO 130
  125          RV(IJ) = A(I,J)
  130          IJ = IJ + N
               II = I + J
               IF (II .GT. N) GO TO 150
               JJ = MB - J
               IF (E21 .GE. 0.0E0) GO TO 140
               II = I
               JJ = MB + J
  140          RV(KJ) = A(II,JJ)
               KJ = KJ + N
  150       CONTINUE
C
  180       RV(IJ) = A(I,MB) - X1
            RV6(I) = EPS4
            IF (ORDER .EQ. 0.0E0) RV6(I) = Z(I,R)
  200    CONTINUE
C
         IF (M1 .EQ. 0) GO TO 600
C     .......... ELIMINATION WITH INTERCHANGES ..........
         DO 580 I = 1, N
            II = I + 1
            MAXK = MIN0(I+M1-1,N)
            MAXJ = MIN0(N-I,M21-2) * N
C
            DO 360 K = I, MAXK
               KJ1 = K
               J = KJ1 + N
               JJ = J + MAXJ
C
               DO 340 KJ = J, JJ, N
                  RV(KJ1) = RV(KJ)
                  KJ1 = KJ
  340          CONTINUE
C
               RV(KJ1) = 0.0E0
  360       CONTINUE
C
            IF (I .EQ. N) GO TO 580
            U = 0.0E0
            MAXK = MIN0(I+M1,N)
            MAXJ = MIN0(N-II,M21-2) * N
C
            DO 450 J = I, MAXK
               IF (ABS(RV(J)) .LT. ABS(U)) GO TO 450
               U = RV(J)
               K = J
  450       CONTINUE
C
            J = I + N
            JJ = J + MAXJ
            IF (K .EQ. I) GO TO 520
            KJ = K
C
            DO 500 IJ = I, JJ, N
               V = RV(IJ)
               RV(IJ) = RV(KJ)
               RV(KJ) = V
               KJ = KJ + N
  500       CONTINUE
C
            IF (ORDER .NE. 0.0E0) GO TO 520
            V = RV6(I)
            RV6(I) = RV6(K)
            RV6(K) = V
  520       IF (U .EQ. 0.0E0) GO TO 580
C
            DO 560 K = II, MAXK
               V = RV(K) / U
               KJ = K
C
               DO 540 IJ = J, JJ, N
                  KJ = KJ + N
                  RV(KJ) = RV(KJ) - V * RV(IJ)
  540          CONTINUE
C
               IF (ORDER .EQ. 0.0E0) RV6(K) = RV6(K) - V * RV6(I)
  560       CONTINUE
C
  580    CONTINUE
C     .......... BACK SUBSTITUTION
C                FOR I=N STEP -1 UNTIL 1 DO -- ..........
  600    DO 630 II = 1, N
            I = N + 1 - II
            MAXJ = MIN0(II,M21)
            IF (MAXJ .EQ. 1) GO TO 620
            IJ1 = I
            J = IJ1 + N
            JJ = J + (MAXJ - 2) * N
C
            DO 610 IJ = J, JJ, N
               IJ1 = IJ1 + 1
               RV6(I) = RV6(I) - RV(IJ) * RV6(IJ1)
  610       CONTINUE
C
  620       V = RV(I)
            IF (ABS(V) .GE. EPS3) GO TO 625
C     .......... SET ERROR -- NEARLY SINGULAR LINEAR SYSTEM ..........
            IF (ORDER .EQ. 0.0E0) IERR = -R
            V = SIGN(EPS3,V)
  625       RV6(I) = RV6(I) / V
  630    CONTINUE
C
         XU = 1.0E0
         IF (ORDER .EQ. 0.0E0) GO TO 870
C     .......... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP ..........
         IF (GROUP .EQ. 0) GO TO 700
C
         DO 680 JJ = 1, GROUP
            J = R - GROUP - 1 + JJ
            XU = 0.0E0
C
            DO 640 I = 1, N
  640       XU = XU + RV6(I) * Z(I,J)
C
            DO 660 I = 1, N
  660       RV6(I) = RV6(I) - XU * Z(I,J)
C
  680    CONTINUE
C
  700    NORM = 0.0E0
C
         DO 720 I = 1, N
  720    NORM = NORM + ABS(RV6(I))
C
         IF (NORM .GE. 0.1E0) GO TO 840
C     .......... IN-LINE PROCEDURE FOR CHOOSING
C                A NEW STARTING VECTOR ..........
         IF (ITS .GE. N) GO TO 830
         ITS = ITS + 1
         XU = EPS4 / (UK + 1.0E0)
         RV6(1) = EPS4
C
         DO 760 I = 2, N
  760    RV6(I) = XU
C
         RV6(ITS) = RV6(ITS) - EPS4 * UK
         GO TO 600
C     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
  830    IERR = -R
         XU = 0.0E0
         GO TO 870
C     .......... NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER ..........
  840    U = 0.0E0
C
         DO 860 I = 1, N
  860    U = U + RV6(I)**2
C
         XU = 1.0E0 / SQRT(U)
C
  870    DO 900 I = 1, N
  900    Z(I,R) = RV6(I) * XU
C
         X0 = X1
  920 CONTINUE
C
 1001 RETURN
      END
