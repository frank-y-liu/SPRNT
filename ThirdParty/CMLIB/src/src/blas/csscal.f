      SUBROUTINE CSSCAL(N,SA,CX,INCX)
C***BEGIN PROLOGUE  CSSCAL
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  D1A6
C***KEYWORDS  BLAS,COMPLEX,LINEAR ALGEBRA,SCALE,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Scale a complex vector
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SA  single precision scale factor
C       CX  complex vector with N elements
C     INCX  storage spacing between elements of CX
C
C     --Output--
C       CX  scaled result (unchanged if N .LE. 0)
C
C     Replace complex CX by (single precision SA) * (complex CX)
C     For I = 0 to N-1, replace CX(1+I*INCX) with  SA * CX(1+I*INCX)
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CSSCAL
C
      COMPLEX CX(*)
      REAL    SA
C***FIRST EXECUTABLE STATEMENT  CSSCAL
      IF(N .LE. 0) RETURN
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          CX(I) = SA*CX(I)
   10     CONTINUE
      RETURN
      END
