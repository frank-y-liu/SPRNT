      COMPLEX FUNCTION CBETA(A,B)
C***BEGIN PROLOGUE  CBETA
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7B
C***KEYWORDS  BETA FUNCTION,COMPLETE BETA FUNCTION,COMPLEX,
C             SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  CBETA computes the complete Beta function of complex
C            parameters A and B.
C***DESCRIPTION
C
C CBETA computes the complete beta function of complex parameters A
C and B.
C Input Parameters:
C       A   complex and the real part of A positive
C       B   complex and the real part of B positive
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CGAMMA,CLBETA,GAMLIM,XERROR
C***END PROLOGUE  CBETA
      COMPLEX A, B, CGAMMA, CLBETA, CEXP
      DATA XMAX / 0.0 /
C***FIRST EXECUTABLE STATEMENT  CBETA
      IF (XMAX.EQ.0.0) CALL GAMLIM (XMIN, XMAX)
C
      IF (REAL(A).LE.0.0 .OR. REAL(B).LE.0.0) CALL XERROR ( 'CBETA   REA
     1L PART OF BOTH ARGUMENTS MUST BE GT 0', 48, 1, 2)
C
      IF (REAL(A)+REAL(B).LT.XMAX) CBETA = CGAMMA(A) * (CGAMMA(B)/
     1  CGAMMA(A+B) )
      IF (REAL(A)+REAL(B).LT.XMAX) RETURN
C
      CBETA = CEXP (CLBETA(A, B))
C
      RETURN
      END
