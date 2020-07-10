      SUBROUTINE VSSINI(N,C1,C2,WSAVE)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C     Modified by R.F. Boisvert, 6 Apr 95 --- new file, was vcosi entry point
C
      DIMENSION C1(N),C2(N),WSAVE(N+15)
      PI=PIMACH(1.0)
      DX=PI/(2*N)
C
C     GENERATE A(I)+-B(I)
C
      DO 100 I=1,N
         C=COS((I-1)*DX)
         S=SIN((I-1)*DX)
         C1(I)=.5*(S+C)
  100    C2(I)=.5*(S-C)
C
C     INITIALIZE VRFFTPK ROUTINE
C
      CALL VRFFTI(N,WSAVE)
      RETURN
      END
