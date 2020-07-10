      SUBROUTINE ELMHES(NM,N,LOW,IGH,A,INT)
C***BEGIN PROLOGUE  ELMHES
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C1B2
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Reduces real general matrix to upper Hessenberg form
C            stabilized elementary similarity transformations.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure ELMHES,
C     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     Given a REAL GENERAL matrix, this subroutine
C     reduces a submatrix situated in rows and columns
C     LOW through IGH to upper Hessenberg form by
C     stabilized elementary similarity transformations.
C
C     On INPUT
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        LOW and IGH are integers determined by the balancing
C          subroutine  BALANC.  If  BALANC  has not been used,
C          set LOW=1, IGH=N.
C
C        A contains the input matrix.
C
C     On OUTPUT
C
C        A contains the Hessenberg matrix.  The multipliers
C          which were used in the reduction are stored in the
C          remaining triangle under the Hessenberg matrix.
C
C        INT contains information on the rows and columns
C          interchanged in the reduction.
C          Only elements LOW through IGH are used.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ELMHES
C
      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
      REAL A(NM,N)
      REAL X,Y
      INTEGER INT(IGH)
C
C***FIRST EXECUTABLE STATEMENT  ELMHES
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C
      DO 180 M = KP1, LA
         MM1 = M - 1
         X = 0.0E0
         I = M
C
         DO 100 J = M, IGH
            IF (ABS(A(J,MM1)) .LE. ABS(X)) GO TO 100
            X = A(J,MM1)
            I = J
  100    CONTINUE
C
         INT(M) = I
         IF (I .EQ. M) GO TO 130
C    .......... INTERCHANGE ROWS AND COLUMNS OF A ..........
         DO 110 J = MM1, N
            Y = A(I,J)
            A(I,J) = A(M,J)
            A(M,J) = Y
  110    CONTINUE
C
         DO 120 J = 1, IGH
            Y = A(J,I)
            A(J,I) = A(J,M)
            A(J,M) = Y
  120    CONTINUE
C    .......... END INTERCHANGE ..........
  130    IF (X .EQ. 0.0E0) GO TO 180
         MP1 = M + 1
C
         DO 160 I = MP1, IGH
            Y = A(I,MM1)
            IF (Y .EQ. 0.0E0) GO TO 160
            Y = Y / X
            A(I,MM1) = Y
C
            DO 140 J = M, N
  140       A(I,J) = A(I,J) - Y * A(M,J)
C
            DO 150 J = 1, IGH
  150       A(J,M) = A(J,M) + Y * A(J,I)
C
  160    CONTINUE
C
  180 CONTINUE
C
  200 RETURN
      END
