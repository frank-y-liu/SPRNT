      SUBROUTINE INVERT(MM, M, A, DET, WORK, IWORK, IERR, OUNIT)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      COMPUTES THE INVERSE AND DETERMINANT OF THE SYMMETRIC MATRIX
C      (E.G., A COVARIANCE MATRIX)
C
C   DESCRIPTION
C   -----------
C
C   1.  THE LINPACK SUBROUTINE SSIFA IS CALLED TO FACTOR THE MATRIX AND
C       THEN THE LINPACK SUBROUTINE SSIDI IS CALLED TO USE THE
C       FACTORIZATION TO FIND THE INVERSE AND DETERMINANT.  THE INPUT
C       MATRIX MUST BE SYMMETRIC AND IS OVERWRITTEN WITH ITS INVERSE ON
C       OUTPUT.
C
C   INPUT PARAMETERS
C   ----------------
C
C   MM    INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE FIRST DIMENSION OF THE MATRIX A.  MUST BE AT LEAST M.
C
C   M     INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE NUMBER OF ROWS AND COLUMNS IN THE MATRIX A.
C
C   A     REAL SYMMETRIC MATRIX WHOSE FIRST DIMENSION MUST BE MM AND
C            WHOSE SECOND DIMENSION MUST BE AT LEAST M (CHANGED ON
C            OUTPUT).
C         THE MATRIX OF DATA VALUES.
C
C         A(I,J) IS THE VALUE FOR THE J-TH VARIABLE FOR THE I-TH CASE.
C
C   WORK  REAL VECTOR DIMENSIONED AT LEAST N.
C         WORK VECTOR.
C
C   IWORK INTEGER VECTOR DIMENSIONED AT LEAST N.
C         WORK VECTOR.
C
C   OUNIT INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         UNIT NUMBER FOR ERROR MESSAGES.
C
C   OUTPUT PARAMETERS
C   -----------------
C
C   A     REAL MATRIX WHOSE FIRST DIMENSION MUST BE MM AND SECOND
C            DIMENSION MUST BE AT LEAST N.
C         THE INVERSE OF THE INPUT MATRIX.
C
C   DET   REAL VECTOR DIMENSIONED AT LEAST 2.
C         THE DETERMINANT OF THE MATRIX.
C
C         THE DETERMINANT IS  DET(1) ** DET(2).
C
C   IERR  INTEGER SCALAR.
C         ERROR FLAG.
C
C         IF IERR = 0, NO ERROR CONDITION WAS DETECTED.
C
C         IF IERR = K, THE K-TH PIVOT BLOCK IS SINGULAR.  THE INVERSE IS
C                      NOT COMPUTED.  ERROR CONDITION SET IN CMLIB
C                      ROUTINE SSIFA.
C
C   REFERENCES
C   ----------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGE 69.
C
C     NBS CORE MATH LIBRARY, VOLS. 1-4 (GAITHERSBURG: QA297.C69 IN NBS
C     LIBRARY, ADMIN E-120).
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      INTEGER OUNIT
      DIMENSION A(MM,*), IWORK(*), WORK(*), DET(*), INERT(3)
C
      CALL SSIFA(A,MM,M,IWORK,IERR)
      IF (IERR .NE. 0) THEN
         IF (OUNIT .GT. 0)
     *      WRITE(OUNIT,*) 'MATRIX TO BE INVERTED MAY BE SINGULAR '
         RETURN
      ENDIF
      JOB = 111
      CALL SSIDI(A,MM,M,IWORK,DET,INERT,WORK,JOB)
      DO 10 I = 1 , M
         DO 10 J = I , M
 10         A(J,I) = A(I,J)
      RETURN
      END
