      SUBROUTINE VSCOSQ(F,L,M,N,FT,C1,C2,C3,C4,WORK)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION F(L,M*N),FT(L,N,M),C1(M),C2(M),C3(M),C4(M),WORK(M+15)
C
C     PREPROCESSING
C
      DO 120 K=1,N
         K1=M*(K-1)+1
         K2=M*K
         DO 100 I=1,L
  100       FT(I,K,1)=F(I,K1)
      IF (2*(M/2) .NE. M) GO TO 120
      DO 110 I=1,L
  110       FT(I,K,M)=-F(I,K2)
  120 CONTINUE
      DO 200 J=2,M-1,2
      JBY2=J/2
      DO 200 K=1,N
         K1=M*(K-1)+J
         K2=K1+1
         DO 200 I=1,L
            FT(I,K,J)   = F(I,K2)*C1(JBY2)+F(I,K1)*C2(JBY2)
  200       FT(I,K,J+1) = -F(I,K2)*C2(JBY2)+F(I,K1)*C1(JBY2)
C
C     REAL,PERIODIC SYNTHESIS
C
      CALL VRFFTB(L*N,M,FT,F,L*N,WORK)
C
C     POSTPROCESSING
C
      DO 300 J=1,M
      DO 300 K=1,N
         K1=N*(J-1)+K
         DO 300 I=1,L
  300       F(I,K1)=C4(J)*FT(I,K,J)+C3(J)*FT(I,K,M+1-J)
      RETURN
      END
