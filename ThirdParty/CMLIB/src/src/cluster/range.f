      SUBROUTINE RANGE(N, X, XMIN, XMAX)
C
C<><><><><><><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><>
C
C   DESCRIPTION
C   -----------
C
C     DETERMINES THE RANGE OF THE VECTOR X.  THE LARGEST VALUE IS
C     RETURNED IN XMAX AND THE SMALLEST VALUE IS RETURNED IN XMIN
C
C   INPUT PARAMETERS
C   ----------------
C
C   N     INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         NUMBER OF ELEMENTS IN VECTOR X.
C
C   X     REAL VECTOR DIMENSIONED AT LEAST N (UNCHANGED ON OUTPUT).
C         VECTOR OF ELEMENTS.
C
C   OUTPUT PARAMETERS
C   -----------------
C
C   XMIN, XMAX  REAL SCALARS
C         THE MINIMUM AND MAXIMUM RESPECTIVELY OF THE VECTOR X.
C
C<><><><><><><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><>
C
      DIMENSION X(*)
C
      XMIN=X(1)
      XMAX=X(1)
      DO 10 I=2,N
         IF(X(I).LT.XMIN) XMIN=X(I)
   10    IF(X(I).GT.XMAX) XMAX=X(I)
      RETURN
      END
