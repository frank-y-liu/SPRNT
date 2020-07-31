 
      COMPLEX FUNCTION CDCDOT(N,CB,CX,INCX,CY,INCY)
C***BEGIN PROLOGUE  CDCDOT
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  830815   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D1A4
C***KEYWORDS  COMPLEX,DOT PRODUCT,DOUBLE PRECISION,INNER PRODUCT
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON,R. J.,(SNLA)
C           KINCAID,D. R.,(U. OF TEXAS)
C           KROGH,F. T.,(JPL)
C***PURPOSE  Complex inner product with result accumulated in d.p.
C***DESCRIPTION
C
C                  B L A S Subprogram
C  Description of Parameters
C
C  --INPUT--
C          N  number of elements in input vector(s)
C         CB  complex scalar to be added to inner product
C         CX  complex vector with N elements
C       INCX  storage spacing between elements of CX
C         CY  complex vector with N elements
C       INCY  storage spacing between elements of CY
C
C  --OUTPUT--
C     CDCDOT  complex dot product (zero if N .LE. 0)
C
C   Returns result with dot product accumulated in D.P.
C     CDCDOT = CB + sum for I = 0 to N-1 of CX(LX+I*INCY)*CY(LY+I*INCY)
C     where LX + 1 if INCX .GE. 0, else LX =(-INCX)*N, AND LY IS
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON D.L, HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CDCDOT
      INTEGER N,INCX,INCY,I,KX,KY
      COMPLEX CX(*),CY(*),CB
      DOUBLE PRECISION DSDOTR,DSDOTI,DT1,DT2,DT3,DT4
C***FIRST EXECUTABLE STATEMENT  CDCDOT
      DSDOTR=REAL(CB)
      DSDOTI=AIMAG(CB)
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
      DSDOTR=DSDOTR+(DT1*DT2)-(DT3*DT4)
      DSDOTI=DSDOTI+(DT1*DT4)+(DT3*DT2)
      KX=KX+INCX
      KY=KY+INCY
    5 CONTINUE
   10 CDCDOT=CMPLX(SNGL(DSDOTR),SNGL(DSDOTI))
      RETURN
      END
