      SUBROUTINE SDATRP(X,XOUT,YOUT,YPOUT,NEQ,KOLD,PHI,PSI)
C
C***BEGIN PROLOGUE  SDATRP
C***REFER TO  SDASSL
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   830315   (YYMMDD)
C***REVISION DATE  830315   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***END PROLOGUE  SDATRP
C
C-----------------------------------------------------------------------
C     the methods in subroutine dastep use polynomials
C     to approximate the solution. sdatrp approximates the
C     solution and its derivative at time xout by evaluating
C     one of these polynomials,and its derivative,there.
C     information defining this polynomial is passed from
C     dastep, so sdatrp cannot be used alone.
C
C     the parameters are%
C     x     the current time in the integration.
C     xout  the time at which the solution is desired
C     yout  the interpolated approximation to y at xout
C           (this is output)
C     ypout the interpolated approximation to yprime at xout
C           (this is output)
C     neq   number of equations
C     kold  order used on last successful step
C     phi   array of scaled divided differences of y
C     psi   array of past stepsize history
C-----------------------------------------------------------------------
C
      DIMENSION YOUT(*),YPOUT(*)
      DIMENSION PHI(NEQ,*),PSI(*)
      KOLDP1=KOLD+1
      TEMP1=XOUT-X
      DO 10 I=1,NEQ
         YOUT(I)=PHI(I,1)
10       YPOUT(I)=0.0E0
      C=1.0E0
      D=0.0E0
      GAMMA=TEMP1/PSI(1)
      DO 30 J=2,KOLDP1
         D=D*GAMMA+C/PSI(J-1)
         C=C*GAMMA
         GAMMA=(TEMP1+PSI(J-1))/PSI(J)
         DO 20 I=1,NEQ
            YOUT(I)=YOUT(I)+C*PHI(I,J)
20          YPOUT(I)=YPOUT(I)+D*PHI(I,J)
30       CONTINUE
      RETURN
C
C------end of subroutine sdatrp------
      END
