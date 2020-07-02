      SUBROUTINE SINTI(N,WSAVE)
C***BEGIN PROLOGUE  SINTI
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  830401   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  Initialize for SINT.
C***DESCRIPTION
C
C  Subroutine SINTI initializes the array WSAVE which is used in
C  subroutine SINT.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  The method
C          is most efficient when N+1 is a product of small primes.
C
C  Output Parameter
C
C  WSAVE   a work array with at least INT(3.5*N+16) locations.
C          Different WSAVE arrays are required for different values
C          of N.  The contents of WSAVE must not be changed between
C          calls of SINT.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  RFFTI
C***END PROLOGUE  SINTI
      DIMENSION       WSAVE(*)
      DATA PI /3.14159265358979/
C***FIRST EXECUTABLE STATEMENT  SINTI
      IF (N .LE. 1) RETURN
      NP1 = N+1
      NS2 = N/2
      DT = PI/FLOAT(NP1)
      KS = N+2
      KF = KS+NS2-1
      FK = 0.
      DO 101 K=KS,KF
         FK = FK+1.
         WSAVE(K) = 2.*SIN(FK*DT)
  101 CONTINUE
      CALL RFFTI (NP1,WSAVE(KF+1))
      RETURN
      END
