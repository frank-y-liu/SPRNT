      SUBROUTINE COSQI(N,WSAVE)
C***BEGIN PROLOGUE  COSQI
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  830401   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  Initialize for COSQF and COSQB.
C***DESCRIPTION
C
C  Subroutine COSQI initializes the array WSAVE which is used in
C  both COSQF and COSQB.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the array to be transformed.  The method
C          is most efficient when N is a product of small primes.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15.
C          The same work array can be used for both COSQF and COSQB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.  The contents of
C          WSAVE must not be changed between calls of COSQF or COSQB.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  RFFTI
C***END PROLOGUE  COSQI
      DIMENSION       WSAVE(*)
      DATA PIH /1.57079632679491/
C***FIRST EXECUTABLE STATEMENT  COSQI
      DT = PIH/FLOAT(N)
      FK = 0.
      DO 101 K=1,N
         FK = FK+1.
         WSAVE(K) = COS(FK*DT)
  101 CONTINUE
      CALL RFFTI (N,WSAVE(N+1))
      RETURN
      END
