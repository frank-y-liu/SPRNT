      REAL FUNCTION SCNRM2(N,CX,INCX)
C***BEGIN PROLOGUE  SCNRM2
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  D1A3B
C***KEYWORDS  BLAS,COMPLEX,LINEAR ALGEBRA,NORM,UNITARY,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Unitary norm of complex vector
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
C   SCNRM2  single precision result (zero if N .LE. 0)
C
C     unitary norm of the complex N-vector stored in CX() with storage
C     increment INCX .
C     If N .LE. 0, return with result = 0.
C     If N .GE. 1, then INCX must be .GE. 1
C
C           C. L. Lawson, 1978 Jan 08
C
C     Four phase method     using two built-in constants that are
C     hopefully applicable to all machines.
C         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
C         CUTHI = minimum of  SQRT(V)      over all known machines.
C     where
C         EPS = smallest no. such that EPS + 1. .GT. 1.
C         U   = smallest positive no.   (underflow limit)
C         V   = largest  no.            (overflow  limit)
C
C     Brief outline of algorithm..
C
C     Phase 1    scans zero components.
C     Move to phase 2 when a component is nonzero and .LE. CUTLO
C     Move to phase 3 when a component is .GT. CUTLO
C     Move to phase 4 when a component is .GE. CUTHI/M
C     where M = N for X() real and M = 2*N for complex.
C
C     Values for CUTLO and CUTHI..
C     From the environmental parameters listed in the IMSL converter
C     document the limiting values are as follows..
C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
C                   Univac and DEC at 2**(-103)
C                   Thus CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
C                   Thus CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
C                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SCNRM2
      LOGICAL IMAG, SCALE
      INTEGER          NEXT
      REAL         CUTLO, CUTHI, HITEST, SUM, XMAX, ABSX, ZERO, ONE
      COMPLEX      CX(*)
      DATA         ZERO, ONE /0.0E0, 1.0E0/
C
      DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
C***FIRST EXECUTABLE STATEMENT  SCNRM2
      IF(N .GT. 0) GO TO 10
         SCNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                                                 BEGIN MAIN LOOP
      DO 210 I=1,NN,INCX
         ABSX = ABS(REAL(CX(I)))
         IMAG = .FALSE.
         GO TO NEXT,(30, 50, 70, 90, 110)
   30 IF( ABSX .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      SCALE = .FALSE.
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF( ABSX .EQ. ZERO) GO TO 200
      IF( ABSX .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 ASSIGN 110 TO NEXT
      SUM = (SUM / ABSX) / ABSX
  105 SCALE = .TRUE.
      XMAX = ABSX
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( ABSX .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF( ABSX .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / ABSX)**2
         XMAX = ABSX
         GO TO 200
C
  115 SUM = SUM + (ABSX/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
   85 ASSIGN 90 TO NEXT
      SCALE = .FALSE.
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
      HITEST = CUTHI/FLOAT( N )
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
   90 IF(ABSX .GE. HITEST) GO TO 100
         SUM = SUM + ABSX**2
  200 CONTINUE
C                  CONTROL SELECTION OF REAL AND IMAGINARY PARTS.
C
      IF(IMAG) GO TO 210
         ABSX = ABS(AIMAG(CX(I)))
         IMAG = .TRUE.
      GO TO NEXT,(  50, 70, 90, 110 )
C
  210 CONTINUE
C
C              END OF MAIN LOOP.
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      SCNRM2 = SQRT(SUM)
      IF(SCALE) SCNRM2 = SCNRM2 * XMAX
  300 CONTINUE
      RETURN
      END
