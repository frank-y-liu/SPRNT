      SUBROUTINE VSSINB(F,L,M,N,FT,C1,C2,WORK)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION F(L,M*N),FT(L,N,M),C1(M),C2(M),WORK(M+15)
C
C     PREPROCESSING
C
      DO 100 J=2,M
      DO 100 K=1,N
         K1=M*(K-1)+J-1
         K2=M*K+1-J
         DO 100 I=1,L
  100       FT(I,K,J)=C1(J)*F(I,K1)-C2(J)*F(I,K2)
      DO 200 K=1,N
         K1=M*K
         DO 200 I=1,L
  200       FT(I,K,1) = .5*F(I,K1)
C
C     REAL,PERIODIC ANALYSIS
C
      CALL VRFFTF(L*N,M,FT,F,L*N,WORK)
C
C     POSTPROCESSING
C
      SCALE=SQRT(2.)
      DO 320 K=1,N
         K1=N*(M-1)+K
         DO 300 I=1,L
  300       F(I,K) = SCALE*FT(I,K,1)
      IF (2*(M/2) .NE. M) GO TO 320
      DO 310 I=1,L
  310       F(I,K1) = -SCALE*FT(I,K,M)
  320 CONTINUE
      DO 400 J=2,M-1,2
      DO 400 K=1,N
         K1=N*(J-1)+K
         K2=K1+N
         DO 400 I=1,L
            F(I,K1) = -SCALE*(FT(I,K,J+1)+FT(I,K,J))
  400       F(I,K2)=-SCALE*(FT(I,K,J+1)-FT(I,K,J))
      RETURN
      END
