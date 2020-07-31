      SUBROUTINE ELTRAN(NM,N,LOW,IGH,A,INT,Z)
C***BEGIN PROLOGUE  ELTRAN
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C4
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Accumulates the stabilized elementary similarity
C            transformations used in the reduction of a real general
C            matrix to upper Hessenberg form by ELMHES.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure ELMTRANS,
C     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     This subroutine accumulates the stabilized elementary
C     similarity transformations used in the reduction of a
C     REAL GENERAL matrix to upper Hessenberg form by  ELMHES.
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
C        A contains the multipliers which were used in the
C          reduction by  ELMHES  in its lower triangle
C          below the subdiagonal.
C
C        INT contains information on the rows and columns
C          interchanged in the reduction by  ELMHES.
C          only elements LOW through IGH are used.
C
C     On OUTPUT
C
C        Z contains the transformation matrix produced in the
C          reduction by  ELMHES.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ELTRAN
C
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      REAL A(NM,IGH),Z(NM,N)
      INTEGER INT(IGH)
C
C***FIRST EXECUTABLE STATEMENT  ELTRAN
      DO 80 I = 1, N
C
         DO 60 J = 1, N
   60    Z(I,J) = 0.0E0
C
         Z(I,I) = 1.0E0
   80 CONTINUE
C
      KL = IGH - LOW - 1
      IF (KL .LT. 1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = 1, KL
         MP = IGH - MM
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
  100    Z(I,MP) = A(I,MP-1)
C
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
C
         DO 130 J = MP, IGH
            Z(MP,J) = Z(I,J)
            Z(I,J) = 0.0E0
  130    CONTINUE
C
         Z(I,MP) = 1.0E0
  140 CONTINUE
C
  200 RETURN
      END
