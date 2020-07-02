      FUNCTION CARG(Z)
C***BEGIN PROLOGUE  CARG
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  A4A
C***KEYWORDS  ARGUMENT,COMPLEX,COMPLEX NUMBER,ELEMENTARY FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the argument of a complex number.
C***DESCRIPTION
C
C CARG(Z) calculates the argument of the complex number Z.  Note
C that CARG returns a real result.  If Z = X+iY, then CARG is ATAN(Y/X),
C except when both X and Y are zero, in which case the result
C will be zero.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CARG
      COMPLEX Z
C***FIRST EXECUTABLE STATEMENT  CARG
      CARG = 0.0
      IF (REAL(Z).NE.0. .OR. AIMAG(Z).NE.0.) CARG =
     1  ATAN2 (AIMAG(Z), REAL(Z))
C
      RETURN
      END
