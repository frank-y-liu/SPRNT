      SUBROUTINE CSROT(N,CX,INCX,CY,INCY,C,S)
C***BEGIN PROLOGUE  CSROT
C***DATE WRITTEN   810223   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D1B10
C***KEYWORDS  BLAS,COMPLEX,LINEAR ALGEBRA,PLANE ROTATION,VECTOR
C***AUTHOR  DONGARRA, J., (ANL)
C***PURPOSE  Applies a plane rotation to complex vectors.
C***DESCRIPTION
C
C     CSROT applies the complex Givens rotation
C
C          (X)   ( C S)(X)
C          (Y) = (-S C)(Y)
C
C     N times where for I = 0,...,N-1
C
C          X = CX(1+I*INCX)
C          Y = CY(1+I*INCY)
C
C     Argument Description
C
C        N      (integer)  number of elements in each vector
C
C        CX     (complex array)  beginning of one vector
C
C        INCX   (integer)  memory spacing of successive elements
C               of vector CX
C
C        CY     (complex array)  beginning of the other vector
C
C        INCY   (integer)  memory spacing of successive elements
C               of vector CY
C
C        C      (real)  cosine term of the rotation
C
C        S      (real)  sine term of the rotation.
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CSROT
C
      COMPLEX CX(*),CY(*),CTEMP
      REAL C,S
      INTEGER I,INCX,INCY,IX,IY,N
C***FIRST EXECUTABLE STATEMENT  CSROT
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CTEMP = C*CX(IX) + S*CY(IY)
        CY(IY) = C*CY(IY) - S*CX(IX)
        CX(IX) = CTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        CTEMP = C*CX(I) + S*CY(I)
        CY(I) = C*CY(I) - S*CX(I)
        CX(I) = CTEMP
   30 CONTINUE
      RETURN
      END
