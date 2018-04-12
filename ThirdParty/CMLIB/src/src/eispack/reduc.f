      SUBROUTINE REDUC(NM,N,A,B,DL,IERR)
C***BEGIN PROLOGUE  REDUC
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C1C
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Reduces generalized symmetric eigenproblem
C            A*X=(LAMBDA)*B*X, to standard symmetric eigenproblem
C            using Cholesky factorization.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure REDUC1,
C     NUM. MATH. 11, 99-110(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971).
C
C     This subroutine reduces the generalized SYMMETRIC eigenproblem
C     Ax=(LAMBDA)Bx, where B is POSITIVE DEFINITE, to the standard
C     symmetric eigenproblem using the Cholesky factorization of B.
C
C     On Input
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrices A and B.  If the Cholesky
C          factor L of B is already available, N should be prefixed
C          with a minus sign.
C
C        A and B contain the real symmetric input matrices.  Only the
C          full upper triangles of the matrices need be supplied.  If
C          N is negative, the strict lower triangle of B contains,
C          instead, the strict lower triangle of its Cholesky factor L.
C
C        DL contains, if N is negative, the diagonal elements of L.
C
C     On Output
C
C        A contains in its full lower triangle the full lower triangle
C          of the symmetric matrix derived from the reduction to the
C          standard form.  The strict upper triangle of A is unaltered.
C
C        B contains in its strict lower triangle the strict lower
C          triangle of its Cholesky factor L.  The full upper
C          triangle of B is unaltered.
C
C        DL contains the diagonal elements of L.
C
C        IERR is set to
C          Zero       for normal return,
C          7*N+1      if B is not positive definite.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  REDUC
C
      INTEGER I,J,K,N,I1,J1,NM,NN,IERR
      REAL A(NM,N),B(NM,N),DL(N)
      REAL X,Y
C
C***FIRST EXECUTABLE STATEMENT  REDUC
      IERR = 0
      NN = IABS(N)
      IF (N .LT. 0) GO TO 100
C     .......... FORM L IN THE ARRAYS B AND DL ..........
      DO 80 I = 1, N
         I1 = I - 1
C
         DO 80 J = I, N
            X = B(I,J)
            IF (I .EQ. 1) GO TO 40
C
            DO 20 K = 1, I1
   20       X = X - B(I,K) * B(J,K)
C
   40       IF (J .NE. I) GO TO 60
            IF (X .LE. 0.0E0) GO TO 1000
            Y = SQRT(X)
            DL(I) = Y
            GO TO 80
   60       B(J,I) = X / Y
   80 CONTINUE
C     .......... FORM THE TRANSPOSE OF THE UPPER TRIANGLE OF INV(L)*A
C                IN THE LOWER TRIANGLE OF THE ARRAY A ..........
  100 DO 200 I = 1, NN
         I1 = I - 1
         Y = DL(I)
C
         DO 200 J = I, NN
            X = A(I,J)
            IF (I .EQ. 1) GO TO 180
C
            DO 160 K = 1, I1
  160       X = X - B(I,K) * A(J,K)
C
  180       A(J,I) = X / Y
  200 CONTINUE
C     .......... PRE-MULTIPLY BY INV(L) AND OVERWRITE ..........
      DO 300 J = 1, NN
         J1 = J - 1
C
         DO 300 I = J, NN
            X = A(I,J)
            IF (I .EQ. J) GO TO 240
            I1 = I - 1
C
            DO 220 K = J, I1
  220       X = X - A(K,J) * B(I,K)
C
  240       IF (J .EQ. 1) GO TO 280
C
            DO 260 K = 1, J1
  260       X = X - A(J,K) * B(I,K)
C
  280       A(I,J) = X / DL(I)
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- B IS NOT POSITIVE DEFINITE ..........
 1000 IERR = 7 * N + 1
 1001 RETURN
      END
