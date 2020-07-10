      SUBROUTINE DSLVMU(P, Y, S, X)
C
C  ***  SET  Y = S * X,  S = P X P SYMMETRIC MATRIX.  ***
C  ***  LOWER TRIANGLE OF  S  STORED ROWWISE.         ***
C
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C  ***  PARAMETER DECLARATIONS  ***
C
      INTEGER P
      DOUBLE PRECISION S(*), X(P), Y(P)
C     DIMENSION S(P*(P+1)/2)
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER I, IM1, J, K
      DOUBLE PRECISION XI
C
C  ***  NO INTRINSIC FUNCTIONS  ***
C
C  ***  EXTERNAL FUNCTION  ***
C
      DOUBLE PRECISION DDOT
C
C-----------------------------------------------------------------------
C
      J = 1
      DO 10 I = 1, P
         Y(I) = DDOT(I, S(J),1,X,1)
         J = J + I
 10      CONTINUE
C
      IF (P .LE. 1) GO TO 999
      J = 1
      DO 40 I = 2, P
         XI = X(I)
         IM1 = I - 1
         J = J + 1
         DO 30 K = 1, IM1
              Y(K) = Y(K) + S(J)*XI
              J = J + 1
 30           CONTINUE
 40      CONTINUE
C
 999  RETURN
C  ***  LAST CARD OF DSLVMU FOLLOWS  ***
      END
