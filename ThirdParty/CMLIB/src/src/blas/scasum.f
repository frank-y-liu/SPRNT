      FUNCTION SCASUM(N,CX,INCX)
C***BEGIN PROLOGUE  SCASUM
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  D1A3A
C***KEYWORDS  ADD,BLAS,COMPLEX,LINEAR ALGEBRA,MAGNITUDE,SUM,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Sum of magnitudes of real and imaginary components of
C            complex vector
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
C   SCASUM  single precision result (zero if N .LE. 0)
C
C     Returns sums of magnitudes of real and imaginary parts of
C     components of CX.  Note that this is not the L1 norm of CX.
C     CASUM = sum from 0 to N-1 of ABS(REAL(CX(1+I*INCX))) +
C             ABS(IMAG(CX(1+I*INCX)))
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SCASUM
      COMPLEX CX(*)
C***FIRST EXECUTABLE STATEMENT  SCASUM
      SCASUM=0.
      IF(N .LE. 0) RETURN
      NS = N*INCX
          DO 10 I=1,NS,INCX
          SCASUM = SCASUM + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
   10     CONTINUE
      RETURN
      END
