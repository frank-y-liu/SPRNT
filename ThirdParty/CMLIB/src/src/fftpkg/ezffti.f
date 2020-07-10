      SUBROUTINE EZFFTI(N,WSAVE)
C***BEGIN PROLOGUE  EZFFTI
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  830401   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  J1A1
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  Initialize EZFFTF and EZFFTB
C***DESCRIPTION
C
C  Subroutine EZFFTI initializes the array WSAVE which is used in
C  both EZFFTF and EZFFTB.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15.
C          The same work array can be used for both EZFFTF and EZFFTB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  EZFFT1
C***END PROLOGUE  EZFFTI
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  EZFFTI
      IF (N .EQ. 1) RETURN
      CALL EZFFT1 (N,WSAVE(2*N+1),WSAVE(3*N+1))
      RETURN
      END
