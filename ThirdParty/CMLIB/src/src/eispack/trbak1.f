      SUBROUTINE TRBAK1(NM,N,A,E,M,Z)
C***BEGIN PROLOGUE  TRBAK1
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C4
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Forms the eigenvectors of real symmetric matrix from
C            eigenvectors of symmetric tridiagonal matrix formed
C            TRED1.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TRBAK1,
C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     This subroutine forms the eigenvectors of a REAL SYMMETRIC
C     matrix by back transforming those of the corresponding
C     symmetric tridiagonal matrix determined by  TRED1.
C
C     On Input
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        A contains information about the orthogonal trans-
C          formations used in the reduction by  TRED1
C          in its strict lower triangle.
C
C        E contains the subdiagonal elements of the tridiagonal
C          matrix in its last N-1 positions.  E(1) is arbitrary.
C
C        M is the number of eigenvectors to be back transformed.
C
C        Z contains the eigenvectors to be back transformed
C          in its first M columns.
C
C     On Output
C
C        Z contains the transformed eigenvectors
C          in its first M columns.
C
C     Note that TRBAK1 preserves vector Euclidean norms.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  TRBAK1
C
      INTEGER I,J,K,L,M,N,NM
      REAL A(NM,N),E(N),Z(NM,M)
      REAL S
C
C***FIRST EXECUTABLE STATEMENT  TRBAK1
      IF (M .EQ. 0) GO TO 200
      IF (N .EQ. 1) GO TO 200
C
      DO 140 I = 2, N
         L = I - 1
         IF (E(I) .EQ. 0.0E0) GO TO 140
C
         DO 130 J = 1, M
            S = 0.0E0
C
            DO 110 K = 1, L
  110       S = S + A(I,K) * Z(K,J)
C     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN TRED1.
C                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            S = (S / A(I,L)) / E(I)
C
            DO 120 K = 1, L
  120       Z(K,J) = Z(K,J) + S * A(I,K)
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
