      SUBROUTINE SCOPYM(N,SX,INCX,SY,INCY)
C***BEGIN PROLOGUE  SCOPYM
C***DATE WRITTEN   801001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D1A5
C***KEYWORDS  BLAS,COPY,VECTOR
C***AUTHOR  KAHANER,DAVID(NBS)
C***PURPOSE  Copy negative of real SX to real SY.
C***DESCRIPTION
C
C       Description of Parameters
C           The * Flags Output Variables
C
C       N   Number of elements in vector(s)
C      SX   Real vector with N elements
C    INCX   Storage spacing between elements of SX
C      SY*  Real negative copy of SX
C    INCY   Storage spacing between elements of SY
C
C      ***  Note that SY = -SX  ***
C
C  Copy negative of real SX to real SY.  For I=0 to N-1,
C   copy  -SX(LX+I*INCX) to SY(LY+I*INCY), where LX=1 if
C   INCX .GE. 0, else LX = (-INCX)*N, and LY is defined
C   in a similar way using INCY.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SCOPYM
      REAL SX(*),SY(*)
C***FIRST EXECUTABLE STATEMENT  SCOPYM
      IF(N.LE.0) RETURN
       IF(INCX.EQ.INCY)  IF(INCX-1) 5,20,60
    5 CONTINUE
C
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS
C
      IX=1
      IY=1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I=1,N
        SY(IY) = -SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN UP LOOP SO REMAINING VECTOR LENGTH IS MULTIPLE OF 7
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I=1,M
        SY(I) = -SX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I= MP1,N,7
        SY(I) = -SX(I)
        SY(I + 1) = -SX(I + 1)
        SY(I + 2) = -SX(I + 2)
        SY(I + 3) = -SX(I + 3)
        SY(I + 4) = -SX(I + 4)
        SY(I + 5) = -SX(I + 5)
        SY(I + 6) = -SX(I + 6)
   50 CONTINUE
      RETURN
C
C          CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = -SX(I)
   70     CONTINUE
      RETURN
      END
