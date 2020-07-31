      SUBROUTINE DNBDI(ABE,LDA,N,ML,MU,IPVT,DET)
C***BEGIN PROLOGUE  DNBDI
C***DATE WRITTEN   800728   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D3A2
C***KEYWORDS  BAND,DOUBLE PRECISION,LINEAR EQUATIONS,NONSYMMETRIC
C***AUTHOR  VOORHEES, E., (LANL)
C***PURPOSE  Computes the determinant of a Band matrix using
C            the factors computed by DNBCO or DNBFA.
C***DESCRIPTION
C
C     DNBDI computes the determinant of a band matrix
C     using the factors computed by DNBCO or DNBFA.
C     If the inverse is needed, use DNBSL  N  times.
C
C     On Entry
C
C        ABE     DOUBLE PRECISION(LDA, NC)
C                the output from DNBCO or DNBFA.
C                NC must be .GE. 2*ML+MU+1 .
C
C        LDA     INTEGER
C                the leading dimension of the array  ABE .
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
C                the pivot vector from DNBCO or DNBFA.
C
C     On Return
C
C        DET     DOUBLE PRECISION(2)
C                determinant of original matrix.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. DABS(DET(1)) .LT. 10.0
C                or  DET(1) = 0.0 .
C
C     SLATEC.  This version dated 07/28/80 .
C     E. A. Voorhees, Los Alamos Scientific Laboratory
C
C     Subroutines and Functions
C
C     Fortran DABS
C***REFERENCES  SUBROUTINE DNBDI WAS DEVELOPED BY GROUP C-3, LOS ALAMOS
C                 SCIENTIFIC LABORATORY, LOS ALAMOS, NM 87545.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DNBDI
      INTEGER LDA,N,ML,MU,IPVT(*)
      DOUBLE PRECISION ABE(LDA,*),DET(2)
C
      DOUBLE PRECISION TEN
      INTEGER I
C***FIRST EXECUTABLE STATEMENT  DNBDI
      DET(1) = 1.0D0
      DET(2) = 0.0D0
      TEN = 10.0D0
      DO 50 I = 1, N
         IF (IPVT(I) .NE. I) DET(1) = -DET(1)
         DET(1) = ABE(I,ML+1)*DET(1)
C     ...EXIT
         IF (DET(1) .EQ. 0.0D0) GO TO 60
   10    IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20
            DET(1) = TEN*DET(1)
            DET(2) = DET(2) - 1.0D0
         GO TO 10
   20    CONTINUE
   30    IF (DABS(DET(1)) .LT. TEN) GO TO 40
            DET(1) = DET(1)/TEN
            DET(2) = DET(2) + 1.0D0
         GO TO 30
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
