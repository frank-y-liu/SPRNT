 
      SUBROUTINE CFFTB(N,C,WSAVE)
C***BEGIN PROLOGUE  CFFTB
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  830401   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  J1A2
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  Unnormalized inverse of CFFTF.
C***DESCRIPTION
C
C  Subroutine CFFTB computes the backward complex discrete Fourier
C  transform (the Fourier synthesis).  Equivalently, CFFTB computes
C  a complex periodic sequence from its Fourier coefficients.
C  The transform is defined below at output parameter C.
C
C  A call of CFFTF followed by a call of CFFTB will multiply the
C  sequence by N.
C
C  The array WSAVE which is used by subroutine CFFTB must be
C  initialized by calling subroutine CFFTI(N,WSAVE).
C
C  Input Parameters
C
C
C  N      the length of the complex sequence C.  The method is
C         more efficient when N is the product of small primes.
C
C  C      a complex array of length N which contains the sequence
C
C  WSAVE   a real work array which must be dimensioned at least 4*N+15
C          in the program that calls CFFTB.  The WSAVE array must be
C          initialized by calling subroutine CFFTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C          The same WSAVE array can be used by CFFTF and CFFTB.
C
C  Output Parameters
C
C  C      For J=1,...,N
C
C             C(J)=the sum from K=1,...,N of
C
C                   C(K)*EXP(I*J*K*2*PI/N)
C
C                         where I=SQRT(-1)
C
C  WSAVE   contains initialization calculations which must not be
C          destroyed between calls of subroutine CFFTF or CFFTB
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CFFTB1
C***END PROLOGUE  CFFTB
      DIMENSION       C(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  CFFTB
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
