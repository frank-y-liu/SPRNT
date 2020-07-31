      SUBROUTINE VSSINF(F,L,M,N,FT,C1,C2,WORK)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION F(L,M*N),FT(L,N,M),C1(M),C2(M),WORK(M+15)
C
C     PREPROCESSING
C
      SCALE=SQRT(2.)
      DO 120 K=1,N
         K1=M*(K-1)+1
         K2=M*K
         DO 100 I=1,L
  100       FT(I,K,1)=SCALE*F(I,K1)
      IF (2*(M/2) .NE. M) GO TO 120
      DO 110 I=1,L
  110       FT(I,K,M)=-SCALE*F(I,K2)
  120 CONTINUE
      SCALE=.5*SCALE
      DO 200 J=2,M-1,2
      DO 200 K=1,N
         K1=M*(K-1)+J
         K2=K1+1
         DO 200 I=1,L
            FT(I,K,J)   = SCALE*(F(I,K2)-F(I,K1))
  200       FT(I,K,J+1) =-SCALE*(F(I,K2)+F(I,K1))
C
C     REAL,PERIODIC SYNTHESIS
C
      CALL VRFFTB(L*N,M,FT,F,L*N,WORK)
C
C     POSTPROCESSING
C
      DO 300 J=2,M
      DO 300 K=1,N
         K1=N*(J-2)+K
         DO 300 I=1,L
  300       F(I,K1)=C1(J)*FT(I,K,J)+C2(J)*FT(I,K,M+2-J)
      DO 400 K=1,N
         K1=N*(M-1)+K
         DO 400 I=1,L
  400       F(I,K1)=FT(I,K,1)
      RETURN
      END
