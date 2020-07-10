      SUBROUTINE COSTI(N,WSAVE)
C***BEGIN PROLOGUE  COSTI
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  830401   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  Initialize for COST.
C***DESCRIPTION
C
C  Subroutine COSTI initializes the array WSAVE which is used in
C  subroutine COST.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  The method
C          is most efficient when N-1 is a product of small primes.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15.
C          Different WSAVE arrays are required for different values
C          of N.  The contents of WSAVE must not be changed between
C          calls of COST.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  RFFTI
C***END PROLOGUE  COSTI
      DIMENSION       WSAVE(*)
      DATA PI /3.14159265358979/
C***FIRST EXECUTABLE STATEMENT  COSTI
      IF (N .LE. 3) RETURN
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
      DT = PI/FLOAT(NM1)
      FK = 0.
      DO 101 K=2,NS2
         KC = NP1-K
         FK = FK+1.
         WSAVE(K) = 2.*SIN(FK*DT)
         WSAVE(KC) = 2.*COS(FK*DT)
  101 CONTINUE
      CALL RFFTI (NM1,WSAVE(N+1))
      RETURN
      END
