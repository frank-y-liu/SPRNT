      SUBROUTINE PICK(K, B1, B2, BB1, BB2)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      CHANGES THE INTERVAL FOR A VARIABLE IN A BLOCK
C
C   INPUT PARAMETERS
C   -----------------
C
C   K     INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE NUMBER OF THE CHANGE.
C
C   B1    REAL SCALAR (UNCHANGED ON OUTPUT).
C         THE ORIGINAL MINIMUM OF THE BLOCK.
C
C   B2    REAL SCALAR (UNCHANGED ON OUTPUT).
C         THE ORIGINAL MAXIMUM OF THE BLOCK.
C
C   OUTPUT PARAMETERS
C   -----------------
C
C   BB1   REAL SCALAR.
C         THE NEW MINIMUM OF THE BLOCK.
C
C   BB2   REAL SCALAR.
C         THE NEW MINIMUM OF THE BLOCK.
C
C   REFERENCE
C   ---------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGE 53.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      BB=(B2-B1)/2.
      IF(K.EQ.1) THEN
         BB1=B1
         BB2=B2-BB
      ELSE IF(K.EQ.2) THEN
         BB1=B1+BB
         BB2=B2+BB
      ELSE IF(K.EQ.3) THEN
         BB1=B1-BB
         BB2=B2-BB
      ELSE IF(K.EQ.4) THEN
         BB1=B1+BB/2
         BB2=B2-BB/2
      ELSE IF(K.EQ.5) THEN
         BB1=B1+BB
         BB2=B2
      ENDIF
      RETURN
      END
