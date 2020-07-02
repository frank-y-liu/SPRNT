      SUBROUTINE VCOSQB(M,N,X,XT,MDIMX,WSAVE)
C***BEGIN PROLOGUE  VCOSQB
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  860701   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FAST FOURIER TRANSFORM, COSINE TRANSFORM, ODD WAVE
C             NUMBERS, MULTIPLE SEQUENCES
C***AUTHOR  BOISVERT, R. F. (NBS)
C***PURPOSE  Normalized inverse of VCOSQF.
C***DESCRIPTION
C
C  Subroutine VCOSQB computes the backward fast Fourier cosine transform
C  of M quarter wave sequences.  That is, cosine series representations
C  with only odd wave numbers.  The transform is defined below at output
C  parameter X.
C
C  The array WSAVE which is used by subroutine VCOSQB must be
C  initialized by calling subroutine VCOSQI(N,WSAVE).
C
C
C  Input Parameters
C
C  M       the number of sequences to be transformed.
C
C  N       the length of the sequences to be transformed.  The method
C          is most efficient when N is a product of small primes.
C
C  X       an array of size at least X(MDIMX,N) which contains the
C          the sequences to be transformed.  The sequences are stored
C          in the ROWS of X.  Thus, the Jth sequence is stored in
C          X(J,I), I=1,..,N.
C
C  XT      a work array of size at least XT(MDIMX,N).
C
C  MDIMX   the first dimension of the array X exactly as it appears in
C          the calling program.
C
C  WSAVE   a work array which must be dimensioned at least 2*N+15
C          in the program that calls VCOSQB.  The WSAVE array must be
C          initialized by calling subroutine VCOSQI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.
C
C  Output Parameters
C
C  X       For I=1,...,N and J=1,...,M
C
C               X(I)= the sum from K=1 to K=N of
C
C                 4*X(K)*COS((2*K-1)*(I-1)*PI/(2*N)) /SQRT(4*N)
C
C  WSAVE   contains initialization calculations which must not
C          be destroyed between calls of VCOSQF or VCOSQB.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VCOSQF followed immediately by a call of
C           of VCOSQB will return the original sequences X.  Thus,
C           VCOSQB is the correctly normalized inverse VCOSQF.
C
C  -----------------------------------------------------------------
C
C  VCOSQB is a straightforward extension of the subprogram COSQB to
C  handle M simultaneous sequences.  COSQB was orginally developed
C  by P. N. Swarztrauber of NCAR.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  VCOSQB
      DIMENSION       X(MDIMX,*), XT(MDIMX,*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  VCOSQB
      IF (M .LE. 0)  GO TO 900
      IF (N .GT. 2)  GO TO 300
      IF (N .EQ. 2)  GO TO 200
      GO TO 900
C
C  CASE  N = 2
C
  200 CONTINUE
      SCALE = 2.0E0*SQRT(0.50E0)
      DO 210 J=1,M
         X1 = SCALE*(X(J,1)+X(J,2))
         X(J,2) = X(J,1)-X(J,2)
         X(J,1) = X1
  210 CONTINUE
      GO TO 900
C
C  CASE N .GT. 2
C
C     ... PREPROCESSING
C
  300 CONTINUE
      NS2 = (N+1)/2
      NP2 = N+2
      DO 310 I=3,N,2
         DO 310 J=1,M
            XIM1 = X(J,I-1)+X(J,I)
            X(J,I) = X(J,I)-X(J,I-1)
            X(J,I-1) = XIM1
  310 CONTINUE
      DO 320 J=1,M
         X(J,1) = X(J,1)+X(J,1)
  320 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) THEN
         DO 330 J=1,M
            X(J,N) = X(J,N)+X(J,N)
  330    CONTINUE
      ENDIF
C
C     ... REAL, PERIODIC TRANSFORM
C
      CALL VRFFTB (M,N,X,XT,MDIMX,WSAVE(N+1))
C
C     ... POSTPROCESSING
C
      DO 340 K=2,NS2
         KC = NP2-K
         DO 340 J=1,M
            XT(J,K) = WSAVE(K-1)*X(J,KC)+WSAVE(KC-1)*X(J,K)
            XT(J,KC) = WSAVE(K-1)*X(J,K)-WSAVE(KC-1)*X(J,KC)
  340 CONTINUE
      IF (MODN .EQ. 0) THEN
         DO 350 J=1,M
            X(J,NS2+1) = WSAVE(NS2)*(X(J,NS2+1)+X(J,NS2+1))
  350    CONTINUE
      ENDIF
      DO 360 K=2,NS2
         KC = NP2-K
         DO 360 J=1,M
            X(J,K) = XT(J,K)+XT(J,KC)
            X(J,KC) = XT(J,K)-XT(J,KC)
  360 CONTINUE
      DO 370 J=1,M
         X(J,1) = X(J,1)+X(J,1)
  370 CONTINUE
C
C     ... NORMALIZATION
C
      SCALE = 0.5
      DO 380 I=1,N
         DO 380 J=1,M
            X(J,I) = SCALE*X(J,I)
  380 CONTINUE
C
C  EXIT
C
  900 CONTINUE
      RETURN
      END
