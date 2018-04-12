      SUBROUTINE COMBAK(NM,LOW,IGH,AR,AI,INT,M,ZR,ZI)
C***BEGIN PROLOGUE  COMBAK
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C4
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Forms eigenvectors of complex general matrix from
C            eigenvectors of upper Hessenberg matrix output from
C            COMHES.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure COMBAK,
C     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     This subroutine forms the eigenvectors of a COMPLEX GENERAL
C     matrix by back transforming those of the corresponding
C     upper Hessenberg matrix determined by  COMHES.
C
C     On INPUT
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        LOW and IGH are integers determined by the balancing
C          subroutine  CBAL.  If  CBAL  has not been used,
C          set LOW=1 and IGH equal to the order of the matrix.
C
C        AR and AI contain the multipliers which were used in the
C          reduction by  COMHES  in their lower triangles
C          below the subdiagonal.
C
C        INT contains information on the rows and columns
C          interchanged in the reduction by  COMHES.
C          Only elements LOW through IGH are used.
C
C        M is the number of eigenvectors to be back transformed.
C
C        ZR and ZI contain the real and imaginary parts,
C          respectively, of the eigenvectors to be
C          back transformed in their first M columns.
C
C     On OUTPUT
C
C        ZR and ZI contain the real and imaginary parts,
C          respectively, of the transformed eigenvectors
C          in their first M columns.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  COMBAK
C
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL AR(NM,IGH),AI(NM,IGH),ZR(NM,M),ZI(NM,M)
      REAL XR,XI
      INTEGER INT(IGH)
C
C***FIRST EXECUTABLE STATEMENT  COMBAK
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
            XR = AR(I,MP-1)
            XI = AI(I,MP-1)
            IF (XR .EQ. 0.0E0 .AND. XI .EQ. 0.0E0) GO TO 110
C
            DO 100 J = 1, M
               ZR(I,J) = ZR(I,J) + XR * ZR(MP,J) - XI * ZI(MP,J)
               ZI(I,J) = ZI(I,J) + XR * ZI(MP,J) + XI * ZR(MP,J)
  100       CONTINUE
C
  110    CONTINUE
C
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
C
         DO 130 J = 1, M
            XR = ZR(I,J)
            ZR(I,J) = ZR(MP,J)
            ZR(MP,J) = XR
            XI = ZI(I,J)
            ZI(I,J) = ZI(MP,J)
            ZI(MP,J) = XI
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
