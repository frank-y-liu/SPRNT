      SUBROUTINE FIRST(DMCOV, N, COV, X, Y)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      COMPUTES LARGEST EIGENVECTOR OF SUBMATRIX OF COV OVER ROWS I
C         WHERE X(I).NE.0. FOR I=1,...,N
C
C   INPUT PARAMETERS
C   ----------------
C
C   DMCOV INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE FIRST DIMENSION OF THE MATRIX COV.  MUST BE AT LEAST N.
C
C   N     INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE NUMBER OF ROWS AND COLUMNS OF MATRIX COV.
C
C   COV   REAL MATRIX WHOSE FIRST DIMENSION MUST BE DMCOV AND WHOSE
C            SECOND DIMENSION MUST BE AT LEAST N. (CHANGED ON OUTPUT).
C         THE COVARIANCE MATRIX.
C
C   X     REAL VECTOR DIMENSIONED AT LEAST N (CHANGED ON OUTPUT).
C         ON INPUT, X(J) .NE. 0. MEANS ROW J WILL BE CONSIDERED IN THE
C            CALCULATION.  OTHERWISE, ROW J WILL NOT BE CONSIDERED.
C
C   Y     REAL VECTOR DIMENSIONED AT LEAST N.
C         WORK VECTOR.
C
C   OUTPUT PARAMETER
C   ----------------
C
C   X     REAL VECTOR DIMENSIONED AT LEAST N.
C         ON OUTPUT, X IS THE EIGENVECTOR.
C
C   REFERENCE
C   ---------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGE 328.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      INTEGER DMCOV
      DIMENSION COV(DMCOV,*), X(*), Y(*)
C
      TH=.0005
      ICNT=0
      DO 20 I=1,N
         IF(COV(I,I).LE.TH) THEN
            DO 10 J = 1 , N
               COV(I,J) = 0.
   10          COV(J,I) = 0.
         ENDIF
         IF (X(I).NE.0.) X(I)=SQRT(COV(I,I))
   20 CONTINUE
   30 SN=0.
      DO 40 I=1,N
   40    SN=SN+X(I)*X(I)
      SN=SQRT(SN)
      DO 50 I=1,N
   50    IF(SN.NE.0.) X(I)=X(I)/SN
      SYY=0.
      SXY=0.
      DO 70 I=1,N
         Y(I)=0.
         IF(X(I).NE.0.) THEN
            DO 60 J=1,N
   60          IF (X(J).NE.0.) Y(I)=Y(I)+COV(I,J)*X(J)
            SXY=SXY+X(I)*Y(I)
            SYY=SYY+Y(I)*Y(I)
         ENDIF
   70 CONTINUE
      IF(SXY.LT.0.) SXY=0.
      E=SQRT(SXY)
      SYY=SQRT(SYY)
      ERR=0.
      DO 80 I=1,N
         IF(SYY.NE.0.) Y(I)=Y(I)/SYY
         ERR=ERR+(X(I)-Y(I))**2
   80    X(I)=E*(1.2*Y(I)-0.2*X(I))
      ICNT=ICNT+1
      IF (ICNT.GT.20) RETURN
      IF (ERR.GT.TH*TH) GO TO 30
      RETURN
      END
