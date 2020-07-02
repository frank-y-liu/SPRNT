      COMPLEX FUNCTION CLNREL(Z)
C***BEGIN PROLOGUE  CLNREL
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY (YYMMDD)
C   000601 Changed CLOG to generic LOG
C***CATEGORY NO.  C4B
C***KEYWORDS  COMPLEX,ELEMENTARY FUNCTION,LOGARITHM,RELATIVE ERROR
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the principal value of the complex natural
C            logarithm of 1+Z with relative error accuracy for small
C            CABS(Z).
C***DESCRIPTION
C
C CLNREL(Z) = CLOG(1+Z) with relative error accuracy near Z = 0.
C Let   RHO = CABS(Z)  and
C       R**2 = CABS(1+Z)**2 = (1+X)**2 + Y**2 = 1 + 2*X + RHO**2 .
C Now if RHO is small we may evaluate CLNREL(Z) accurately by
C       CLOG(1+Z) = CMPLX  (ALOG(R), CARG(1+Z))
C                 = CMPLX  (0.5*ALOG(R**2), CARG(1+Z))
C                 = CMPLX  (0.5*ALNREL(2*X+RHO**2), CARG(1+Z))
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALNREL,CARG,R1MACH,XERROR
C***END PROLOGUE  CLNREL
      COMPLEX Z
      EXTERNAL CARG
      DATA SQEPS /0.0/
C***FIRST EXECUTABLE STATEMENT  CLNREL
      IF (SQEPS.EQ.0.) SQEPS = SQRT (R1MACH(4))
C
      IF (CABS(1.+Z).LT.SQEPS) CALL XERROR ( 'CLNREL  ANSWER LT HALF PRE
     1CISION BECAUSE Z TOO NEAR -1', 54,    1, 1)
C
      RHO = CABS(Z)
      IF (RHO.GT.0.375) CLNREL = LOG (1.0+Z)
      IF (RHO.GT.0.375) RETURN
C
      X = REAL(Z)
      CLNREL = CMPLX (0.5*ALNREL(2.*X+RHO**2), CARG(1.0+Z))
C
      RETURN
      END
