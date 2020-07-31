      SUBROUTINE SINQI(N,WSAVE)
C***BEGIN PROLOGUE  SINQI
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  830401   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  Initialize for SINQF and SINQB.
C***DESCRIPTION
C
C  Subroutine SINQI initializes the array WSAVE which is used in
C  both SINQF and SINQB.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  The method
C          is most efficient when N is a product of small primes.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15.
C          The same work array can be used for both SINQF and SINQB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.  The contents of
C          WSAVE must not be changed between calls of SINQF or SINQB.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  COSQI
C***END PROLOGUE  SINQI
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINQI
      CALL COSQI (N,WSAVE)
      RETURN
      END
