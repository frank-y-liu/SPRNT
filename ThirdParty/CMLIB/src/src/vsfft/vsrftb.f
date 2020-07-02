      SUBROUTINE VSRFTB(F,L,M,N,FT,WSAVE)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION F(L,N*M),FT(L,N,M),WSAVE(M+15)
C
C     RE-ORDER INPUT
C
      DO 120 K=1,N
      K1=M*(K-1)+1
      K2=M*K
      DO 100 I=1,L
  100 FT(I,K,1)=F(I,K1)
      IF (2*(M/2) .NE. M) GO TO 120
      DO 110 I=1,L
  110 FT(I,K,M)=F(I,K2)
  120 CONTINUE
      DO 200 J=2,M-1,2
      DO 200 K=1,N
      K1=M*(K-1)+J
      K2=K1+1
      DO 200 I=1,L
      FT(I,K,J)=F(I,K1)
  200 FT(I,K,J+1)=F(I,K2)
C
C     REAL, PERIODIC TRANSFORM
C
      CALL VRFFTB(L*N,M,FT,F,L*N,WSAVE)
      DO 300 J=1,M
      DO 300 K=1,N
      K1=N*(J-1)+K
      DO 300 I=1,L
  300 F(I,K1)=FT(I,K,J)
      RETURN
      END
