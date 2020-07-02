      SUBROUTINE VSCSQI(N,C1,C2,C3,C4,WSAVE)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      ENTRY VSSNQI(N,C1,C2,C3,C4,WSAVE)
      DIMENSION C1(N),C2(N),C3(N),C4(N),WSAVE(N+15)
      PI=PIMACH(1.0)
      DX=PI/N
      SCALE=SQRT(.5)
C
C     GENERATE A(I)+-B(I)
C
      DO 100 I=1,(N-1)/2
         C=COS(I*DX)
         S=SIN(I*DX)
         C1(I)=.5*(S+C)
  100    C2(I)=.5*(C-S)
C
      DX=PI/(2*N)
      DO 200 I=1,N
         C=COS((I-.5)*DX)
         S=SIN((I-.5)*DX)
         C3(I)=SCALE*(C+S)
200      C4(I)=SCALE*(C-S)
C
C     INITIALIZE VRFFTPK ROUTINE
C
      CALL VRFFTI(N,WSAVE)
      RETURN
      END
