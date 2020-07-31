      SUBROUTINE CORTB(NM,LOW,IGH,AR,AI,ORTR,ORTI,M,ZR,ZI)
C***BEGIN PROLOGUE  CORTB
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C4
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Forms eigenvectors of complex general matrix from
C            eigenvectors of upper Hessenberg matrix output from
C             CORTH
C***DESCRIPTION
C
C     This subroutine is a translation of a complex analogue of
C     the ALGOL procedure ORTBAK, NUM. MATH. 12, 349-368(1968)
C     by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     This subroutine forms the eigenvectors of a COMPLEX GENERAL
C     matrix by back transforming those of the corresponding
C     upper Hessenberg matrix determined by  CORTH.
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
C        AR and AI contain information about the unitary
C          transformations used in the reduction by  CORTH
C          in their strict lower triangles.
C
C        ORTR and ORTI contain further information about the
C          transformations used in the reduction by  CORTH.
C          Only elements LOW through IGH are used.
C
C        M is the number of columns of ZR and ZI to be back transformed.
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
C        ORTR and ORTI have been altered.
C
C     Note that CORTB preserves vector Euclidean norms.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CORTB
C
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL AR(NM,IGH),AI(NM,IGH),ORTR(IGH),ORTI(IGH)
      REAL ZR(NM,M),ZI(NM,M)
      REAL H,GI,GR
C
C***FIRST EXECUTABLE STATEMENT  CORTB
      IF (M .EQ. 0) GO TO 200
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = KP1, LA
         MP = LOW + IGH - MM
         IF (AR(MP,MP-1) .EQ. 0.0E0 .AND. AI(MP,MP-1) .EQ. 0.0E0)
     1      GO TO 140
C     .......... H BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
         H = AR(MP,MP-1) * ORTR(MP) + AI(MP,MP-1) * ORTI(MP)
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
            ORTR(I) = AR(I,MP-1)
            ORTI(I) = AI(I,MP-1)
  100    CONTINUE
C
         DO 130 J = 1, M
            GR = 0.0E0
            GI = 0.0E0
C
            DO 110 I = MP, IGH
               GR = GR + ORTR(I) * ZR(I,J) + ORTI(I) * ZI(I,J)
               GI = GI + ORTR(I) * ZI(I,J) - ORTI(I) * ZR(I,J)
  110       CONTINUE
C
            GR = GR / H
            GI = GI / H
C
            DO 120 I = MP, IGH
               ZR(I,J) = ZR(I,J) + GR * ORTR(I) - GI * ORTI(I)
               ZI(I,J) = ZI(I,J) + GR * ORTI(I) + GI * ORTR(I)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
