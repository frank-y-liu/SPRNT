      SUBROUTINE FIND(DMCOV, N, COV, Y, X, XY)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      FINDS BEST EIGENVECTOR FITTING TO ARRAY COV
C
C   INPUT PARAMETERS
C   ----------------
C
C   DMCOV INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE FIRST DIMENSION OF THE MATRIX COV.  MUST BE AT LEAST N.
C
C   N     INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE NUMBER OF VARIABLES.
C
C   COV   REAL MATRIX WHOSE FIRST DIMENSION MUST BE DMCOV AND WHOSE
C            SECOND DIMENSION MUST BE AT LEAST N. (CHANGED ON OUTPUT).
C         THE COVARIANCE MATRIX.
C
C   X     REAL VECTOR DIMENSIONED AT LEAST N.
C         WORK VECTOR.
C
C   XY    REAL VECTOR DIMENSIONED AT LEAST N.
C         WORK VECTOR.
C
C   OUTPUT PARAMETER
C   ----------------
C
C   Y     REAL VECTOR DIMENSIONED AT LEAST N.
C         THE BEST EIGENVECTOR.
C
C   REFERENCE
C   ---------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGE 329.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      INTEGER DMCOV
      DIMENSION COV(DMCOV,*), X(*), Y(*), XY(*)
C
      DO 10 I=1,N
   10    X(I)=1.
      CMAX=0.
C
C     CALCULATE LARGEST EIGENVECTOR
C
   20 CALL FIRST(DMCOV,N,COV,X,XY)
      XS=0.
      SS=0.
      DO 30 I=1,N
         IF (X(I).NE.0.) THEN
            XS=XS+1.
            SS=SS+X(I)*X(I)
         ENDIF
   30 CONTINUE
      XMIN=R1MACH(2)
C
C     CHOOSE IMIN AS IN STEP 3 OF ALGORITHM
C
      DO 50 J=1,N
         IF(X(J).NE.0.) THEN
            XX=0.
            DO 40 I=1,N
   40          XX=XX+X(I)*COV(I,J)
            XX=XX*XX/(SS*SS*COV(J,J))
            IF (XX.LE.XMIN) THEN
               XMIN=XX
               IMIN=J
            ENDIF
         ENDIF
   50 CONTINUE
      IF(XS.NE.0.) CC=SS/XS
      IF(CC.GT.CMAX) THEN
         CMAX=CC
         DO 60 I=1,N
   60       Y(I)=X(I)
      ENDIF
C
C     REMOVE IMIN AND COMPUTE PARTIAL CORRELATION MATRIX OF COV
C
      CALL REMOVE(DMCOV,N,COV,IMIN,X)
      X(IMIN)=0.
      IF (XS.GT.1.) GO TO 20
C
C     STOP IF ONLY ONE VARIABLE LEFT
C
      RETURN
      END
