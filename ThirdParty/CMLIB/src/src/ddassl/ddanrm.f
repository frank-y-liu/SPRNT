      DOUBLE PRECISION FUNCTION DDANRM(NEQ,V,WT,RPAR,IPAR)
C
C***BEGIN PROLOGUE  DDANRM
C***REFER TO  DDASSL
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   830315   (YYMMDD)
C***REVISION DATE  830315   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***END PROLOGUE  DDANRM
C-----------------------------------------------------------------------
C     this function routine computes the weighted
C     root-mean-square norm of the vector of length
C     neq contained in the array v,with weights
C     contained in the array wt of length neq.
C        ddanrm=sqrt((1/neq)*sum(v(i)/wt(i))**2)
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION V(NEQ),WT(NEQ)
      DIMENSION RPAR(*),IPAR(*)
      DDANRM = 0.0D0
      VMAX = 0.0D0
      DO 10 I = 1,NEQ
10      IF(DABS(V(I)/WT(I)) .GT. VMAX) VMAX = DABS(V(I)/WT(I))
      IF(VMAX .LE. 0.0D0) GO TO 30
      SUM = 0.0D0
      DO 20 I = 1,NEQ
20      SUM = SUM + ((V(I)/WT(I))/VMAX)**2
      DDANRM = VMAX*DSQRT(SUM/DBLE(NEQ))
30    CONTINUE
      RETURN
C------end of function ddanrm------
      END
