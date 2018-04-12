      FUNCTION CSEVL(X,CS,N)
C***BEGIN PROLOGUE  CSEVL
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  C3A2
C***KEYWORDS  CHEBYSHEV,FNLIB,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Evaluate the N-term Chebyshev series CS at X.
C***DESCRIPTION
C
C Evaluate the N-term Chebyshev series CS at X.  Adapted from
C R. Broucke, Algorithm 446, C.A.C.M., 16, 254 (1973). Also see Fox
C and Parker, Chebyshev Polynomials in Numerical Analysis, Oxford Press,
C page 56.
C
C       Input Arguments --
C X    value at which the series is to be evaluated.
C CS   array of N terms of a Chebyshev series.  In eval-
C      uating CS, only half the first coefficient is summed.
C N    number of terms in array CS.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  CSEVL
C
       DIMENSION CS(*)
C***FIRST EXECUTABLE STATEMENT  CSEVL
       IF(N.LT.1) CALL XERROR( 'CSEVL   NUMBER OF TERMS LE 0', 28, 2,2)
       IF(N.GT.1000) CALL XERROR ( 'CSEVL   NUMBER OF TERMS GT 1000',
     1   31,3,2)
       IF (X.LT. -1.0 .OR. X.GT. 1.0) CALL XERROR( 'CSEVL   X OUTSIDE (-
     11,+1)', 25, 1, 1)
C
       B1=0.
       B0=0.
       TWOX=2.*X
       DO 10 I=1,N
       B2=B1
       B1=B0
       NI=N+1-I
       B0=TWOX*B1-B2+CS(NI)
 10    CONTINUE
C
       CSEVL = 0.5 * (B0-B2)
C
       RETURN
      END
