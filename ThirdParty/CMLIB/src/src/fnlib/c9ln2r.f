      COMPLEX FUNCTION C9LN2R(Z)
C***BEGIN PROLOGUE  C9LN2R
C***DATE WRITTEN   780401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C4B
C***KEYWORDS  COMPLEX,ELEMENTARY FUNCTION,LOGARITHM,RELATIVE ERROR,
C             SECOND ORDER
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Evaluates CLOG(1+Z) from second order with relative error
C            so that  CLOG(1+Z) = Z - Z**2/2 + Z**3*C9LN2R(Z).
C***DESCRIPTION
C
C Evaluate  CLOG(1+Z)  from 2-nd order with relative error accuracy so
C that     CLOG(1+Z) = Z - Z**2/2 + Z**3*C9LN2R(Z).
C
C Now  CLOG(1+Z) = 0.5*ALOG(1+2*X+CABS(Z)**2) + I*CARG(1+Z),
C where X = REAL(Z)  and  Y = AIMAG(Z).
C We find
C     Z**3 * C9LN2R(Z) = -X*CABS(Z)**2 - 0.25*CABS(Z)**4
C        + (2*X+CABS(Z)**2)**3 * R9LN2R(2*X+CABS(Z)**2)
C        + I * (CARG(1+Z) + (X-1)*Y)
C The imaginary part must be evaluated carefully as
C     (ATAN(Y/(1+X)) - Y/(1+X)) + Y/(1+X) - (1-X)*Y
C       = (Y/(1+X))**3 * R9ATN1(Y/(1+X)) + X**2*Y/(1+X)
C
C Now we divide through by Z**3 carefully.  Write
C     1/Z**3 = (X-I*Y)/CABS(Z)**3 * (1/CABS(Z)**3)
C then   C9LN2R(Z) = ((X-I*Y)/CABS(Z))**3 * (-X/CABS(Z) - CABS(Z)/4
C        + 0.5*((2*X+CABS(Z)**2)/CABS(Z))**3 * R9LN2R(2*X+CABS(Z)**2)
C        + I*Y/(CABS(Z)*(1+X)) * ((X/CABS(Z))**2 +
C          + (Y/(CABS(Z)*(1+X)))**2 * R9ATN1(Y/(1+X)) ) )
C
C If we let  XZ = X/CABS(Z)  and  YZ = Y/CABS(Z)  we may write
C     C9LN2R(Z) = (XZ-I*YZ)**3 * (-XZ - CABS(Z)/4
C        + 0.5*(2*XZ+CABS(Z))**3 * R9LN2R(2*X+CABS(Z)**2)
C        + I*YZ/(1+X) * (XZ**2 + (YZ/(1+X))**2*R9ATN1(Y/(1+X)) ))
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R9ATN1,R9LN2R
C***END PROLOGUE  C9LN2R
      COMPLEX Z, CLOG
C***FIRST EXECUTABLE STATEMENT  C9LN2R
      X = REAL (Z)
      Y = AIMAG (Z)
C
      CABSZ = CABS(Z)
      IF (CABSZ.GT.0.8125) GO TO 20
C
      C9LN2R = CMPLX (1.0/3.0, 0.0)
      IF (CABSZ.EQ.0.0) RETURN
C
      XZ = X/CABSZ
      YZ = Y/CABSZ
C
      ARG = 2.0*XZ + CABSZ
      RPART = 0.5*ARG**3*R9LN2R(CABSZ*ARG) - XZ - 0.25*CABSZ
      Y1X = YZ/(1.0+X)
      AIPART = Y1X * (XZ**2 + Y1X**2*R9ATN1(CABSZ*Y1X) )
C
      C9LN2R = CMPLX(XZ,-YZ)**3 * CMPLX(RPART,AIPART)
      RETURN
C
 20   C9LN2R = (CLOG(1.0+Z) - Z*(1.0-0.5*Z)) / Z**3
      RETURN
C
      END
