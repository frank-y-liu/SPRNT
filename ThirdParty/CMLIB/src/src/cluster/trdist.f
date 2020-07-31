      SUBROUTINE TRDIST(M, DMD, D, DMT1, DMT2, T)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      PRODUCES TRIADS FROM A DISTANCE MATRIX
C
C   DESCRIPTION
C   -----------
C
C   1.  (IJ)K IS A LEGAL TRIAD IF THE DISTANCE BETWEEN OBJECTS I AND J
C       IS SMALLER THAN THE DISTANCE BETWEEN OBJECTS I AND K AND THE
C       DISTANCE BETWEEN OBJECTS J AND K.
C
C   2.  THE ROUTINE RETURNS THE MATRIX T.  T(I,J,K)=1.  IF (IJ)K IS A
C       LEGAL TRIAD.  T(I,J,K)=0.  OTHERWISE.
C
C   INPUT PARAMETERS
C   ----------------
C
C   M     INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE NUMBER OF OBJECTS.
C
C   DMD   INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE LEADING DIMENSION OF THE MATRIX D.  MUST BE AT LEAST M.
C
C   D     REAL MATRIX WHOSE FIRST DIMENSION MUST BE DMD AND SECOND
C            DIMENSION MUST BE AT LEAST M (UNCHANGED ON OUTPUT).
C         D(I,J) IS THE DISTANCE FROM OBJECT I TO OBJECT J.
C
C   DMT1  INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE FIRST DIMENSION OF THE MATRIX T.  MUST BE AT LEAST M.
C
C   DMT2  INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE SECOND DIMENSION OF THE MATRIX T.  MUST BE AT LEAST M.
C
C   OUTPUT PARAMETER
C   ----------------
C
C   T     REAL MATRIX WHOSE FIRST DIMENSION MUST BE DMT1, WHOSE SECOND
C            DIMENSION MUST BE AT DMT2, AND WHOSE THIRD DIMENSION MUST
C            BE AT LEAST M
C         T(I,J,K) = 1. IF (IJ)K IS A LEGAL TRIAD, IE. IF THE DISTANCE
C         FROM OBJECT I TO OBJECT J IS SMALLER THAN THE DISTANCE FROM
C         OBJECT J TO OBJECT K AND THE DISTANCE FROM OBJECT I TO OBJECT
C         K.
C
C   REFERENCE
C   ---------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGE 190.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      INTEGER DMD, DMT1, DMT2
      DIMENSION D(DMD,*), T(DMT1,DMT2,*)
C
      DO 10 I=1,M
         DO 10 J=1,M
            DO 10 K=1,M
               T(I,J,K)=0.
               IF(D(I,J).LT.D(J,K).AND.D(I,J).LT.D(K,I)) T(I,J,K)=1.
   10 CONTINUE
      RETURN
      END