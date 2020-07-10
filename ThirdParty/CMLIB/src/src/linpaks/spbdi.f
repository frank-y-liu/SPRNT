      SUBROUTINE SPBDI(ABD,LDA,N,M,DET)
C***BEGIN PROLOGUE  SPBDI
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D3B2
C***KEYWORDS  BANDED,DETERMINANT,FACTOR,INVERSE,LINEAR ALGEBRA,LINPACK,
C             MATRIX,POSITIVE DEFINITE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Computes the determinant of a real SYMMETRIC POSITIVE DEFI-
C            NITE BAND matrix using the factors computed by SPBCO or
C            SPBFA.  If the inverse is needed, use SPBSL  N  times.
C***DESCRIPTION
C
C     SPBDI computes the determinant
C     of a real symmetric positive definite band matrix
C     using the factors computed by SPBCO or SPBFA.
C     If the inverse is needed, use SPBSL  N  times.
C
C     On Entry
C
C        ABD     REAL(LDA, N)
C                the output from SPBCO or SPBFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        M       INTEGER
C                the number of diagonals above the main diagonal.
C
C     On Return
C
C        DET     REAL(2)
C                determinant of original matrix in the form
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. DET(1) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SPBDI
      INTEGER LDA,N,M
      REAL ABD(LDA,*)
      REAL DET(2)
C
      REAL S
      INTEGER I
C
C     COMPUTE DETERMINANT
C
C***FIRST EXECUTABLE STATEMENT  SPBDI
      DET(1) = 1.0E0
      DET(2) = 0.0E0
      S = 10.0E0
      DO 50 I = 1, N
         DET(1) = ABD(M+1,I)**2*DET(1)
C     ...EXIT
         IF (DET(1) .EQ. 0.0E0) GO TO 60
   10    IF (DET(1) .GE. 1.0E0) GO TO 20
            DET(1) = S*DET(1)
            DET(2) = DET(2) - 1.0E0
         GO TO 10
   20    CONTINUE
   30    IF (DET(1) .LT. S) GO TO 40
            DET(1) = DET(1)/S
            DET(2) = DET(2) + 1.0E0
         GO TO 30
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
