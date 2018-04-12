      SUBROUTINE SGBDI(ABD,LDA,N,ML,MU,IPVT,DET)
C***BEGIN PROLOGUE  SGBDI
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D3A2
C***KEYWORDS  BANDED,DETERMINANT,FACTOR,INVERSE,LINEAR ALGEBRA,LINPACK,
C             MATRIX
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Computes the determinant of a BAND matrix
C            using the factors computed by SGBCO or SGBFA.
C            If the inverse is needed, use SGBSL  N  times.
C***DESCRIPTION
C
C     SGBDI computes the determinant of a band matrix
C     using the factors computed by SBGCO or SGBFA.
C     If the inverse is needed, use SGBSL  N  times.
C
C     On Entry
C
C        ABD     REAL(LDA, N)
C                the output from SBGCO or SGBFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C
C        IPVT    INTEGER(N)
C                the pivot vector from SBGCO or SGBFA.
C
C     On Return
C
C        DET     REAL(2)
C                determinant of original matrix.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. ABS(DET(1)) .LT. 10.0
C                or  DET(1) = 0.0 .
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     Fortran ABS
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SGBDI
      INTEGER LDA,N,ML,MU,IPVT(*)
      REAL ABD(LDA,*),DET(2)
C
      REAL TEN
      INTEGER I,M
C
C***FIRST EXECUTABLE STATEMENT  SGBDI
      M = ML + MU + 1
      DET(1) = 1.0E0
      DET(2) = 0.0E0
      TEN = 10.0E0
      DO 50 I = 1, N
         IF (IPVT(I) .NE. I) DET(1) = -DET(1)
         DET(1) = ABD(M,I)*DET(1)
C     ...EXIT
         IF (DET(1) .EQ. 0.0E0) GO TO 60
   10    IF (ABS(DET(1)) .GE. 1.0E0) GO TO 20
            DET(1) = TEN*DET(1)
            DET(2) = DET(2) - 1.0E0
         GO TO 10
   20    CONTINUE
   30    IF (ABS(DET(1)) .LT. TEN) GO TO 40
            DET(1) = DET(1)/TEN
            DET(2) = DET(2) + 1.0E0
         GO TO 30
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
