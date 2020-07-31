      REAL FUNCTION SDSDOT(N,SB,SX,INCX,SY,INCY)
C***BEGIN PROLOGUE  SDSDOT
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  D1A4
C***KEYWORDS  BLAS,INNER PRODUCT,LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  S.P. result with inner product accumulated in d.p.
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SB  single precision scalar to be added to inner product
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C
C     --Output--
C   SDSDOT  single precision dot product (zero if N .LE. 0)
C
C     Returns S.P. result with dot product accumulated in D.P.
C     SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY)
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SDSDOT
C
      REAL              SX(*),SY(*),SB
      DOUBLE PRECISION DSDOT
C***FIRST EXECUTABLE STATEMENT  SDSDOT
      DSDOT = DBLE(SB)
      IF(N .LE. 0) GO TO 30
      IF(INCX.EQ.INCY.AND.INCX.GT.0) GO TO 40
      KX = 1
      KY = 1
      IF(INCX.LT.0) KX = 1+(1-N)*INCX
      IF(INCY.LT.0) KY = 1+(1-N)*INCY
          DO 10 I = 1,N
          DSDOT = DSDOT + DBLE(SX(KX))*DBLE(SY(KY))
          KX = KX + INCX
          KY = KY + INCY
   10     CONTINUE
   30 SDSDOT = SNGL(DSDOT)
      RETURN
   40 CONTINUE
      NS = N*INCX
          DO 50 I=1,NS,INCX
          DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I))
   50     CONTINUE
      SDSDOT = SNGL(DSDOT)
      RETURN
      END
