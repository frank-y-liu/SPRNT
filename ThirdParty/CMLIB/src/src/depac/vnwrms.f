      REAL FUNCTION VNWRMS(N,V,W)
C***BEGIN PROLOGUE  VNWRMS
C***REFER TO  DEBDF
C
C   VNWRMS computes a weighted root-mean-square vector norm for the
C   integrator package DEBDF.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  VNWRMS
C
C
CLLL. OPTIMIZE
C-----------------------------------------------------------------------
C THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED ROOT-MEAN-SQUARE NORM
C OF THE VECTOR OF LENGTH N CONTAINED IN THE ARRAY V, WITH WEIGHTS
C CONTAINED IN THE ARRAY W OF LENGTH N..
C   VNWRMS = SQRT( (1/N) * SUM( V(I)/W(I) )**2 )
C-----------------------------------------------------------------------
      INTEGER N, I
      REAL V, W, SUM
      DIMENSION V(N), W(N)
C***FIRST EXECUTABLE STATEMENT  VNWRMS
      SUM = 0.0E0
      DO 10 I = 1,N
 10     SUM = SUM + (V(I)/W(I))**2
      VNWRMS = SQRT(SUM/FLOAT(N))
      RETURN
C----------------------- END OF FUNCTION VNWRMS ------------------------
      END
