      SUBROUTINE DROT(N,DX,INCX,DY,INCY,DC,DS)
C***BEGIN PROLOGUE  DROT
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  D1A8
C***KEYWORDS  BLAS,GIVENS ROTATION,LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  APPLY d.p. Givens rotation
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C       DC  D.P. element of rotation matrix
C       DS  D.P. element of rotation matrix
C
C     --Output--
C       DX  rotated vector (unchanged if N .LE. 0)
C       DY  rotated vector (unchanged if N .LE. 0)
C
C     Multiply the 2 X 2 matrix  ( DC DS) times the 2 X N matrix (DX**T)
C                                (-DS DC)                        (DY**T)
C     where **T indicates transpose.  The elements of DX are in
C     DX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else
C     LX = (-INCX)*N, and similarly for DY using LY and INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DROT
C
      DOUBLE PRECISION DX,DY,DC,DS,ZERO,ONE,W,Z
      DIMENSION DX(*),DY(*)
      DATA ZERO,ONE/0.D0,1.D0/
C***FIRST EXECUTABLE STATEMENT  DROT
      IF(N .LE. 0 .OR. (DS .EQ. ZERO .AND. DC .EQ. ONE)) GO TO 40
      IF(.NOT. (INCX .EQ. INCY .AND. INCX .GT. 0)) GO TO 20
C
           NSTEPS=INCX*N
           DO 10 I=1,NSTEPS,INCX
                W=DX(I)
                Z=DY(I)
                DX(I)=DC*W+DS*Z
                DY(I)=-DS*W+DC*Z
   10           CONTINUE
           GO TO 40
C
   20 CONTINUE
           KX=1
           KY=1
C
           IF(INCX .LT. 0) KX=1-(N-1)*INCX
           IF(INCY .LT. 0) KY=1-(N-1)*INCY
C
           DO 30 I=1,N
                W=DX(KX)
                Z=DY(KY)
                DX(KX)=DC*W+DS*Z
                DY(KY)=-DS*W+DC*Z
                KX=KX+INCX
                KY=KY+INCY
   30           CONTINUE
   40 CONTINUE
C
      RETURN
      END
