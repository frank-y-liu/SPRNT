      SUBROUTINE DDAWTS(NEQ,IWT,RTOL,ATOL,Y,WT,RPAR,IPAR)
C
C***BEGIN PROLOGUE  DDAWTS
C***REFER TO  DDASSL
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   830315   (YYMMDD)
C***REVISION DATE  830315   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***END PROLOGUE  DDAWTS
C-----------------------------------------------------------------------
C     this subroutine sets the error weight vector
C     wt according to wt(i)=rtol(i)*abs(y(i))+atol(i),
C     i=1,-,n.
C     rtol and atol are scalars if iwt = 0,
C     and vectors if iwt = 1.
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RTOL(*),ATOL(*),Y(*),WT(*)
      DIMENSION RPAR(*),IPAR(*)
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 20 I=1,NEQ
         IF (IWT .EQ.0) GO TO 10
           RTOLI=RTOL(I)
           ATOLI=ATOL(I)
10         WT(I)=RTOLI*DABS(Y(I))+ATOLI
20         CONTINUE
      RETURN
C-----------end of subroutine ddawts------------------------------------
      END
