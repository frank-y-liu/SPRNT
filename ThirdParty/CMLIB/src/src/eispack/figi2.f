      SUBROUTINE FIGI2(NM,N,T,D,E,Z,IERR)
C***BEGIN PROLOGUE  FIGI2
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C1C
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Transforms certain real non-symmetric tridiagonal matrix
C            to symmetric tridiagonal matrix.
C***DESCRIPTION
C
C     Given a NONSYMMETRIC TRIDIAGONAL matrix such that the products
C     of corresponding pairs of off-diagonal elements are all
C     non-negative, and zero only when both factors are zero, this
C     subroutine reduces it to a SYMMETRIC TRIDIAGONAL matrix
C     using and accumulating diagonal similarity transformations.
C
C     On INPUT
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        T contains the input matrix.  Its subdiagonal is
C          stored in the last N-1 positions of the first column,
C          its diagonal in the N positions of the second column,
C          and its superdiagonal in the first N-1 positions of
C          the third column.  T(1,1) and T(N,3) are arbitrary.
C
C     On OUTPUT
C
C        T is unaltered.
C
C        D contains the diagonal elements of the symmetric matrix.
C
C        E contains the subdiagonal elements of the symmetric
C          matrix in its last N-1 positions.  E(1) is not set.
C
C        Z contains the transformation matrix produced in
C          the reduction.
C
C        IERR is set to
C          Zero       for normal return,
C          N+I        if T(I,1)*T(I-1,3) is negative,
C          2*N+I      if T(I,1)*T(I-1,3) is zero with
C                     one factor non-zero.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  FIGI2
C
      INTEGER I,J,N,NM,IERR
      REAL T(NM,3),D(N),E(N),Z(NM,N)
      REAL H
C
C***FIRST EXECUTABLE STATEMENT  FIGI2
      IERR = 0
C
      DO 100 I = 1, N
C
         DO 50 J = 1, N
   50    Z(I,J) = 0.0E0
C
         IF (I .EQ. 1) GO TO 70
         H = T(I,1) * T(I-1,3)
         IF (H) 900, 60, 80
   60    IF (T(I,1) .NE. 0.0E0 .OR. T(I-1,3) .NE. 0.0E0) GO TO 1000
         E(I) = 0.0E0
   70    Z(I,I) = 1.0E0
         GO TO 90
   80    E(I) = SQRT(H)
         Z(I,I) = Z(I-1,I-1) * E(I) / T(I-1,3)
   90    D(I) = T(I,2)
  100 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL
C                ELEMENTS IS NEGATIVE ..........
  900 IERR = N + I
      GO TO 1001
C     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL
C                ELEMENTS IS ZERO WITH ONE MEMBER NON-ZERO ..........
 1000 IERR = 2 * N + I
 1001 RETURN
      END
