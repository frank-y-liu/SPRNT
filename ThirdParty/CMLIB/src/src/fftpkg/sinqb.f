      SUBROUTINE SINQB(N,X,WSAVE)
C***BEGIN PROLOGUE  SINQB
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  830401   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  Unnormalized inverse of SINQF.
C***DESCRIPTION
C
C  Subroutine SINQB computes the fast Fourier transform of quarter
C  wave data.  That is, SINQB computes a sequence from its
C  representation in terms of a sine series with odd wave numbers.
C  the transform is defined below at output parameter X.
C
C  SINQF is the unnormalized inverse of SINQB since a call of SINQB
C  followed by a call of SINQF will multiply the input sequence X
C  by 4*N.
C
C  The array WSAVE which is used by subroutine SINQB must be
C  initialized by calling subroutine SINQI(N,WSAVE).
C
C
C  Input Parameters
C
C  N       the length of the array X to be transformed.  The method
C          is most efficient when N is a product of small primes.
C
C  X       an array which contains the sequence to be transformed
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15
C          in the program that calls SINQB.  The WSAVE array must be
C          initialized by calling subroutine SINQI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C
C  Output Parameters
C
C  X       for I=1,...,N
C
C               X(I)= the sum from K=1 to K=N of
C
C                 4*X(K)*SIN((2k-1)*I*PI/(2*N))
C
C               a call of SINQB followed by a call of
C               SINQF will multiply the sequence X by 4*N.
C               Therefore SINQF is the unnormalized inverse
C               of SINQB.
C
C  WSAVE   contains initialization calculations which must not
C          be destroyed between calls of SINQB or SINQF.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  COSQB
C***END PROLOGUE  SINQB
      DIMENSION       X(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINQB
      IF (N .GT. 1) GO TO 101
      X(1) = 4.*X(1)
      RETURN
  101 NS2 = N/2
      DO 102 K=2,N,2
         X(K) = -X(K)
  102 CONTINUE
      CALL COSQB (N,X,WSAVE)
      DO 103 K=1,NS2
         KC = N-K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
  103 CONTINUE
      RETURN
      END
