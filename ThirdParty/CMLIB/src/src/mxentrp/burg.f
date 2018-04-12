 
C
C
C
C
      SUBROUTINE BURG(LX,X,F,B,LA,A,M,S,Y,TABLE)
C
C
C***BEGIN PROLOGUE  BURG
C***DATE WRITTEN  1979 ?
C***REVISION DATE  11 OCTOBER 1983
C***CATEGORY NO.  L10F
C***KEYWORDS  SPECTRAL ANALYSIS, MAXIMUM ENTROPY METHOD
C***AUTHOR   MANUAL T. SILVIA
C            U.S. NAVAL UNDERWATER SYSTEMS CENTER
C            NEWPORT, RHODE ISLAND
C***PURPOSE  BURG CALCULATES THE MAXIMUM ENTROPY SPECTRUM
C            OF A REAL, UNIVARIATE TIME SERIES (OF EQUALLY
C            SPACED SAMPLE POINTS).
C***DESCRIPTION
C
C   BURG computes the coefficients of a finite length causal
C   forward or backward prediction filter, and uses both the
C   forward and backward predictions in a symmetric manner to
C   generate the Maximum Entropy Spectrum by means of a Toeplitz
C   recursion.
C
C   This program was copied from the book "Deconvolution of
C   Geophysical Time Series in the Exploration for Oil and
C   Natural Gas", by Manuel T. Silvia and Enders A. Robinson,
C   Elsevier Scientific Publishing Company, New York, 1979.
C   A more complete description of the algorithm and program
C   can be found in that work.  A few very minor modifications
C   were made in October, 1983 by Bert W. Rust of the Scientific
C   Computing Division, National Bureau of Standards.
C
C
C   I N P U T
C   ---------
C
C   LX      Length of the input time series.
C
C   X       Vector (of length LX) containing the input time series.
C
C   LA      Length of prediction filter to use in calculating
C           the spectrum.  LA should be less than LX.
C
C   M       A small integer which determines the number of frequency
C           points at which the spectrum is estimated.  The spectrum
C           will be estimated at  1 + 2**(M-1)  equally spaced
C           frequencies, with the first frequency being 0.0 and the
C           final frequency being the Nyquist frequency (reciprocal of
C           twice the sample interval).  Thus if the sample interval
C           is 1 time unit and M = 9, then the spectral estimate
C           will be calculated at 257 equally spaced frequencies
C           between  0.0  and  0.5.
C
C
C   O U T P U T
C   -----------
C
C   S       Vector of length  (2**(M-1) + 1)  containing the
C           Maximum Entropy spectral estimates.
C
C
C   M I S C E L L A N E O U S
C   -------------------------
C
C   F , B   Working vectors of length LX.
C
C   A , Y   Working vectors of length  2**M.
C
C   TABLE   Working vector of length  2**(M-2) + 1.
C
C
C***REFERENCE  M.T. Silvia and E.A. Robinson, "Deconvolution of
C                Geophysical Time Series in the Exploration for
C                Oil and Natural Gas," Elsevier, New York, 1979.
C
C***ROUTINES CALLED  COSQTB, FFTBRG
C***END PROLOGUE  BURG
C
C
C     NP=2**M CONTROLS SPECTRAL RESOLUTION
      DIMENSION X(LX),F(LX),B(LX),A(2**M),S(2**(M-1)+1),Y(2**M),
     1          TABLE((2**M)/4+1)
      NP=2**M
      V =0.0
         DO 1 I=1,LX
         V=V+X(I)**2
         F(I)=X(I)
    1    B(I)=X(I)
      V=V/FLOAT(LX)
         DO 5 N=1,LA
         SN=0.
         SD=0.
         L=N+1
            DO 2 J=L,LX
            SN=SN+F(J)*B(J-1)
    2       SD=SD+F(J)*F(J)+B(J-1)*B(J-1)
         G=2.*SN/SD
         V=V*(1.-G*G)
         A(N)=G
         IF(N.EQ.1) GO TO 4
         NH=N/2
            DO 3 K=1,NH
            KK=N-K
            AA=A(KK)-G*A(K)
            A(K)=A(K)-G*A(KK)
    3       A(KK)=AA
    4    M1=N+1
         BS=B(N)
            DO 5 K=M1,LX
            FF=F(K)
            BB=BS
            BS=B(K)
            F(K)=FF-G*BB
    5    B(K)=BB-G*FF
      DO 6 I=1,LA
      K=LA-I+1
    6 A(K+1)=-A(K)
      A(1)=1.0
         DO 7 I =1,NP
    7    Y(I)=0.
      LP=LA+2
         DO 8 I=LP,NP
    8    A(I)=0
      CALL COSQTB(M,TABLE)
      XX=(1.0/ALOG(2.0))*ALOG(FLOAT(LA+1))
         LL=IFIX(XX)
      IF(XX/FLOAT(LL).GT.1.) LL=LL+1
      CALL FFTBRG(A,Y,TABLE,M,LL,-1)
      NF=NP/2+1
      DO 9 I=1,NF
    9 S(I)=V/(A(I)*A(I)+Y(I)*Y(I))
      RETURN
      END
