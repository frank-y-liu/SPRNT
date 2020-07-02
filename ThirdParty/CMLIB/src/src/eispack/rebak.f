      SUBROUTINE REBAK(NM,N,B,DL,M,Z)
C***BEGIN PROLOGUE  REBAK
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C4
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Forms eigenvectors of generalized symmetric eigensystem
C            from eigenvectors of derived matrix output from REDUC
C            REDUC2.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure REBAKA,
C     NUM. MATH. 11, 99-110(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971).
C
C     This subroutine forms the eigenvectors of a generalized
C     SYMMETRIC eigensystem by back transforming those of the
C     derived symmetric matrix determined by  REDUC.
C
C     On Input
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix system.
C
C        B contains information about the similarity transformation
C          (Cholesky decomposition) used in the reduction by  REDUC
C          in its strict lower triangle.
C
C        DL contains further information about the transformation.
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
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  REBAK
C
      INTEGER I,J,K,M,N,I1,II,NM
      REAL B(NM,N),DL(N),Z(NM,M)
      REAL X
C
C***FIRST EXECUTABLE STATEMENT  REBAK
      IF (M .EQ. 0) GO TO 200
C
      DO 100 J = 1, M
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
         DO 100 II = 1, N
            I = N + 1 - II
            I1 = I + 1
            X = Z(I,J)
            IF (I .EQ. N) GO TO 80
C
            DO 60 K = I1, N
   60       X = X - B(K,I) * Z(K,J)
C
   80       Z(I,J) = X / DL(I)
  100 CONTINUE
C
  200 RETURN
      END
