      SUBROUTINE PLACE(I1, I2, J1, J2, DMF, F, D)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C     REPLACES THE FREQUENCIES WITHIN A BLOCK WITH A NULL CHARACTER
C
C   DESCRIPTION
C   -----------
C
C   1.  THE ROUTINE REPLACES ALL FREQUENCIES WITHIN THE INPUT BLOCK TO A
C       -D.
C
C   INPUT PARAMETERS
C   ----------------
C
C   I1, I2, J1, J2, DMF, F -- SEE SUBROUTINE DENSTY
C
C   D     REAL SCALAR (UNCHANGED ON OUTPUT).
C         THE NUMBER TO BE REPLACED AS THE NEW FREQUENCY.
C
C   REFERENCE
C   ---------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGE 51.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      INTEGER DMF
      DIMENSION F(DMF,*)
C
      DO 10 I=I1,I2
         DO 10 J=J1,J2
   10       IF(F(I,J).GE.0.) F(I,J)=-D
      RETURN
      END
