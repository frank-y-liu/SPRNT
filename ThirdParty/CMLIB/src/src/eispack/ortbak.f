      SUBROUTINE ORTBAK(NM,LOW,IGH,A,ORT,M,Z)
C***BEGIN PROLOGUE  ORTBAK
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C4
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Forms eigenvectors of general real matrix from
C            eigenvectors of upper Hesenberg matrix output from ORTHES
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure ORTBAK,
C     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     This subroutine forms the eigenvectors of a REAL GENERAL
C     matrix by back transforming those of the corresponding
C     upper Hessenberg matrix determined by  ORTHES.
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
C        A contains information about the orthogonal trans-
C          formations used in the reduction by  ORTHES
C          in its strict lower triangle.
C
C        ORT contains further information about the trans-
C          formations used in the reduction by  ORTHES.
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
C        ORT has been altered.
C
C     NOTE that ORTBAK preserves vector Euclidean norms.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ORTBAK
C
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL A(NM,IGH),ORT(IGH),Z(NM,M)
      REAL G
C
C***FIRST EXECUTABLE STATEMENT  ORTBAK
      IF (M .EQ. 0) GO TO 200
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = KP1, LA
         MP = LOW + IGH - MM
         IF (A(MP,MP-1) .EQ. 0.0E0) GO TO 140
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
  100    ORT(I) = A(I,MP-1)
C
         DO 130 J = 1, M
            G = 0.0E0
C
            DO 110 I = MP, IGH
  110       G = G + ORT(I) * Z(I,J)
C     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN ORTHES.
C                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            G = (G / ORT(MP)) / A(MP,MP-1)
C
            DO 120 I = MP, IGH
  120       Z(I,J) = Z(I,J) + G * ORT(I)
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
