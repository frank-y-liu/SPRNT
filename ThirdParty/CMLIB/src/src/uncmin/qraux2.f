      SUBROUTINE QRAUX2(NR,N,R,I,A,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C PRE-MULTIPLY R BY THE JACOBI ROTATION J(I,I+1,A,B)
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF MATRIX
C R(N,N)      <--> UPPER HESSENBERG MATRIX
C I            --> INDEX OF ROW
C A            --> SCALAR
C B            --> SCALAR
C
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
      DIMENSION R(NR,*)
      DEN=SQRT(A*A + B*B)
      C=A/DEN
      S=B/DEN
      DO 10 J=I,N
        Y=R(I,J)
        Z=R(I+1,J)
        R(I,J)=C*Y - S*Z
        R(I+1,J)=S*Y + C*Z
   10 CONTINUE
      RETURN
      END
