      INTEGER FUNCTION ICAMAX(N,CX,INCX)
C***BEGIN PROLOGUE  ICAMAX
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  D1A2
C***KEYWORDS  BLAS,COMPLEX,LINEAR ALGEBRA,MAXIMUM COMPONENT,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Find largest component of complex vector
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       CX  complex vector with N elements
C     INCX  storage spacing between elements of CX
C
C     --Output--
C   ICAMAX  smallest index (zero if N .LE. 0)
C
C      Returns the index of the component of CX having the
C      largest sum of magnitudes of real and imaginary parts.
C     ICAMAX = first I, I = 1 to N, to minimize
C        ABS(REAL(CX(1-INCX+I*INCX))) + ABS(IMAG(CX(1-INCX+I*INCX)))
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ICAMAX
C
      COMPLEX CX(*)
C***FIRST EXECUTABLE STATEMENT  ICAMAX
      ICAMAX = 0
      IF(N.LE.0) RETURN
      ICAMAX = 1
      IF(N .LE. 1) RETURN
      NS = N*INCX
      II = 1
      SUMMAX = ABS(REAL(CX(1))) + ABS(AIMAG(CX(1)))
          DO 20 I=1,NS,INCX
          SUMRI = ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
          IF(SUMMAX.GE.SUMRI) GO TO 10
          SUMMAX = SUMRI
          ICAMAX = II
   10     II = II + 1
   20     CONTINUE
      RETURN
      END
