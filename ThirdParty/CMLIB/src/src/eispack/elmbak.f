      SUBROUTINE ELMBAK(NM,LOW,IGH,A,INT,M,Z)
C***BEGIN PROLOGUE  ELMBAK
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C4
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Forms eigenvectors of real general matrix from
C            eigenvectors of upper Hessenberg matrix output from
C            ELMHES.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure ELMBAK,
C     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     This subroutine forms the eigenvectors of a REAL GENERAL
C     matrix by back transforming those of the corresponding
C     upper Hessenberg matrix determined by  ELMHES.
C
C     On INPUT
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        LOW and IGH are integers determined by the balancing
C          subroutine  BALANC.  If  BALANC  has not been used,
C          set LOW=1 and IGH equal to the order of the matrix.
C
C        A contains the multipliers which were used in the
C          reduction by  ELMHES  in its lower triangle
C          below the subdiagonal.
C
C        INT contains information on the rows and columns
C          interchanged in the reduction by  ELMHES.
C          Only elements LOW through IGH are used.
C
C        M is the number of columns of Z to be back transformed.
C
C        Z contains the real and imaginary parts of the eigen-
C          vectors to be back transformed in its first M columns.
C
C     On OUTPUT
C
C        Z contains the real and imaginary parts of the
C          transformed eigenvectors in its first M columns.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ELMBAK
C
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL A(NM,IGH),Z(NM,M)
      REAL X
      INTEGER INT(IGH)
C
C***FIRST EXECUTABLE STATEMENT  ELMBAK
      IF (M .EQ. 0) GO TO 200
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = KP1, LA
         MP = LOW + IGH - MM
         MP1 = MP + 1
C
         DO 110 I = MP1, IGH
            X = A(I,MP-1)
            IF (X .EQ. 0.0E0) GO TO 110
C
            DO 100 J = 1, M
  100       Z(I,J) = Z(I,J) + X * Z(MP,J)
C
  110    CONTINUE
C
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
C
         DO 130 J = 1, M
            X = Z(I,J)
            Z(I,J) = Z(MP,J)
            Z(MP,J) = X
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
