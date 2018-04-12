      COMPLEX FUNCTION CLBETA(A,B)
C***BEGIN PROLOGUE  CLBETA
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7B
C***KEYWORDS  BETA FUNCTION,COMPLETE BETA FUNCTION,COMPLEX,LOGARITHM,
C             SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  CLBETA computes the natural log of the complex valued
C            complete Beta function of complex parameters A and B.
C***DESCRIPTION
C
C CLBETA computes the natural log of the complex valued complete beta
C function of complex parameters A and B.  This is a preliminary version
C which is not accurate.
C
C Input Parameters:
C       A   complex and the real part of A positive
C       B   complex and the real part of B positive
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CLNGAM,XERROR
C***END PROLOGUE  CLBETA
      COMPLEX A, B, CLNGAM
C***FIRST EXECUTABLE STATEMENT  CLBETA
      IF (REAL(A).LE.0.0 .OR. REAL(B).LE.0.0) CALL XERROR ( 'CLBETA  REA
     1L PART OF BOTH ARGUMENTS MUST BE GT 0', 48, 1, 2)
C
      CLBETA = CLNGAM(A) + CLNGAM(B) - CLNGAM(A+B)
C
      RETURN
      END
