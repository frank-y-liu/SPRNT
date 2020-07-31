      FUNCTION PRVEC(M,U,V)
C***BEGIN PROLOGUE  PRVEC
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***REFER TO  BVSUP
C
C  This subroutine computes the inner product of a vector U
C  with the imaginary product or mate vector corresponding to V
C***ROUTINES CALLED  SDOT
C***END PROLOGUE  PRVEC
C
      DIMENSION U(*),V(*)
C***FIRST EXECUTABLE STATEMENT  PRVEC
      N=M/2
      NP=N+1
      VP=SDOT(N,U(1),1,V(NP),1)
      PRVEC=SDOT(N,U(NP),1,V(1),1) - VP
      RETURN
      END
