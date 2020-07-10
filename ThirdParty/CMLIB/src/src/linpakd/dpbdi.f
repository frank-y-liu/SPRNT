      SUBROUTINE DPBDI(ABD,LDA,N,M,DET)
C***BEGIN PROLOGUE  DPBDI
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D3B2
C***KEYWORDS  BANDED,DETERMINANT,DOUBLE PRECISION,FACTOR,INVERSE,
C             LINEAR ALGEBRA,LINPACK,MATRIX,POSITIVE DEFINITE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Computes the determinant of a d.p. SYMMETRIC POSITIVE
C            DEFINITE BAND matrix using factors of DPBCO or DPBFA.
C            Use DPBSL N times for inverse.
C***DESCRIPTION
C
C     DPBDI computes the determinant
C     of a double precision symmetric positive definite band matrix
C     using the factors computed by DPBCO or DPBFA.
C     If the inverse is needed, use DPBSL  N  times.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                the output from DPBCO or DPBFA.
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
C        DET     DOUBLE PRECISION(2)
C                determinant of original matrix in the form
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. DET(1) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DPBDI
      INTEGER LDA,N,M
      DOUBLE PRECISION ABD(LDA,*)
      DOUBLE PRECISION DET(2)
C
      DOUBLE PRECISION S
      INTEGER I
C
C     COMPUTE DETERMINANT
C
C***FIRST EXECUTABLE STATEMENT  DPBDI
      DET(1) = 1.0D0
      DET(2) = 0.0D0
      S = 10.0D0
      DO 50 I = 1, N
         DET(1) = ABD(M+1,I)**2*DET(1)
C     ...EXIT
         IF (DET(1) .EQ. 0.0D0) GO TO 60
   10    IF (DET(1) .GE. 1.0D0) GO TO 20
            DET(1) = S*DET(1)
            DET(2) = DET(2) - 1.0D0
         GO TO 10
   20    CONTINUE
   30    IF (DET(1) .LT. S) GO TO 40
            DET(1) = DET(1)/S
            DET(2) = DET(2) + 1.0D0
         GO TO 30
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
