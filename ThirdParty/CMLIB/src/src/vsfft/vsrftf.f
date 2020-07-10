      SUBROUTINE VSRFTF(F,L,M,N,FT,WSAVE)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION F(L,M*N),FT(L,N,M),WSAVE(M+15)
C
C     RE-ORDER INPUT
C
      DO 100 K=1,N
      DO 100 J=1,M
      K1=M*(K-1)+J
      DO 100 I=1,L
  100 FT(I,K,J)=F(I,K1)
C
C     REAL, PERIODIC TRANSFORM
C
      CALL VRFFTF(L*N,M,FT,F,L*N,WSAVE)
      DO 220 K=1,N
      K1=N*(M-1)+K
      DO 200 I=1,L
  200 F(I,K)=FT(I,K,1)
      IF (2*(M/2) .NE. M) GO TO 220
      DO 210 I=1,L
  210 F(I,K1)=FT(I,K,M)
  220 CONTINUE
      DO 300 K=1,N
      DO 300 J=2,M-1,2
      K1=N*(J-1)+K
      K2=K1+N
      DO 300 I=1,L
      F(I,K1)=FT(I,K,J)
  300 F(I,K2)=FT(I,K,J+1)
      RETURN
      END
