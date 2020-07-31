      SUBROUTINE DCDOT(N,FM,CX,INCX,CY,INCY,DCR,DCI)
C***BEGIN PROLOGUE  DCDOT
C***DATE WRITTEN   790101   (YYMMDD)
C***REVISION DATE  830516   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D1A4
C***KEYWORDS  COMPLEX VECTORS,DOT PRODUCT
C***AUTHOR  (NONE)
C***PURPOSE  Computes the dot product of 2 complex vectors, CX and CY,
C            e.g. CX DOT CY, or, CXconjugate DOT CY. The real and imagin
C            parts of CX and CY are converted to double precision, the
C            dot product accumulation is done in double precision and
C            the output is given as 2 double precision numbers,
C            corresponding to the real and imaginary part of the result.
C***DESCRIPTION
C
C  Computes the dot product of 2 complex vectors, CX and CY, e.g.
C    CX DOT CY, or, CXconjugate DOT CY.  The real and imaginary
C    parts of CX and CY are converted to double precision, the dot
C    product accumulation is done in double precision and the output
C    is given as 2 double precision numbers, corresponding to the real
C    and imaginary part of the result.
C     Input
C      N:  Number of complex components of CX and CY.
C      FM: =+1.0   compute CX DOT CY.
C          =-1.0   compute CXconjugate DOT CY.
C      CX(N):
C      CY(N):  Complex arrays of length N.
C      INCX:(Integer)   Spacing of elements of CX to use
C      INCY:(Integer)   Spacing of elements of CY to use.
C     Output
C      DCR:(Double Precision) Real part of dot product.
C      DCI:(Double Precision) Imaginary part of dot product.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DCDOT
      INTEGER N,INCX,INCY,KX,KY,I
      COMPLEX CX(*),CY(*)
      DOUBLE PRECISION DCR,DCI,DT1,DT2,DT3,DT4,FM
C***FIRST EXECUTABLE STATEMENT  DCDOT
      DCR=0.D0
      DCI=0.D0
      IF(N.LE.0) GO TO 10
      KX=1
      KY=1
      IF(INCX.LT.0) KX=1+(1-N)*INCX
      IF(INCY.LT.0) KY=1+(1-N)*INCY
      DO 5 I=1,N
        DT1=DBLE(REAL(CX(KX)))
        DT2=DBLE(REAL(CY(KY)))
        DT3=DBLE(AIMAG(CX(KX)))
        DT4=DBLE(AIMAG(CY(KY)))
      DCR=DCR+(DT1*DT2)-FM*(DT3*DT4)
      DCI=DCI+(DT1*DT4)+FM*(DT3*DT2)
      KX=KX+INCX
      KY=KY+INCY
    5 CONTINUE
   10 CONTINUE
      RETURN
      END
