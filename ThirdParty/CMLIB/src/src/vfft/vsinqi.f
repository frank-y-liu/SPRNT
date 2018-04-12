      SUBROUTINE VSINQI(N,WSAVE)
C***BEGIN PROLOGUE  VSINQI
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  860701   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FAST FOURIER TRANSFORM, SINE TRANSFORM, ODD WAVE
C             NUMBERS, MULTIPLE SEQUENCES
C***AUTHOR  BOISVERT, R. F. (NBS)
C***PURPOSE  Initialize for VSINQF and VSINQB.
C***DESCRIPTION
C
C  Subroutine VSINQI initializes the array WSAVE which is used in
C  both VSINQF and VSINQB.  The prime factorization of N together with
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
C          The same work array can be used for both VSINQF and VSINQB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.  The contents of
C
C          WSAVE must not be changed between calls of VSINQF or VSINQB.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VCOSQI
C***END PROLOGUE  VSINQI
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  VSINQI
      CALL VCOSQI (N,WSAVE)
      RETURN
      END
