      SUBROUTINE U11US(A,MDA,M,N,UB,DB,MODE,NP,KRANK,KSURE,H,W,EB,IR,
     1   IC)
C***BEGIN PROLOGUE  U11US
C***REFER TO  ULSIA
C***ROUTINES CALLED  ISAMAX,SAXPY,SDOT,SNRM2,SSCAL,SSWAP,XERROR
C***DESCRIPTION
C
C       This routine performs an LQ factorization of the
C       matrix A using Householder transformations. Row
C       and column pivots are chosen to reduce the growth
C       of round-off and to help detect possible rank
C       deficiency.
*
*-------------------------------------------------------------------------
* Modified by Ron Boisvert on 10 June 1998
* Calls to SSWAP with integer arrays not legal; problems when integer data
* size is different than real data (in autodoubling compilers, for example)
*-------------------------------------------------------------------------
*
C***END PROLOGUE  U11US
      DIMENSION A(MDA,N),UB(M),DB(M),H(M),W(M),EB(M)
      INTEGER IC(N),IR(M)
C
C        INITIALIZATION
C
C***FIRST EXECUTABLE STATEMENT  U11US
      J=0
      KRANK=M
      DO 10 I=1,N
      IC(I)=I
   10 CONTINUE
      DO 12 I=1,M
      IR(I)=I
   12 CONTINUE
C
C        DETERMINE REL AND ABS ERROR VECTORS
C
C
C
C        CALCULATE ROW LENGTH
C
      DO 30 I=1,M
      H(I)=SNRM2(N,A(I,1),MDA)
      W(I)=H(I)
   30 CONTINUE
C
C         INITIALIZE ERROR BOUNDS
C
      DO  40 I=1,M
      EB(I)=AMAX1(DB(I),UB(I)*H(I))
      UB(I)=EB(I)
      DB(I)=0.0
   40 CONTINUE
C
C          DISCARD SELF DEPENDENT ROWS
C
      I=1
   50 IF(EB(I).GE.H(I)) GO TO 60
      IF(I.EQ.KRANK) GO TO 70
      I=I+1
      GO TO 50
C
C          MATRIX REDUCTION
C
   60 CONTINUE
      KK=KRANK
      KRANK=KRANK-1
      IF(MODE.EQ.0) RETURN
      IF(I.GT.NP) GO TO  64
      CALL XERROR( 'FIRST NP ROWS ARE LINEARLY DEPENDENT',36,
     18,0)
      KRANK=I-1
      RETURN
   64 CONTINUE
      IF(I.GT.KRANK) GO TO 70
      CALL SSWAP(1,EB(I),1,EB(KK),1)
      CALL SSWAP(1,UB(I),1,UB(KK),1)
      CALL SSWAP(1,W(I),1,W(KK),1)
      CALL SSWAP(1,H(I),1,H(KK),1)
*      CALL SSWAP(1,IR(I),1,IR(KK),1)
      iitemp = ir(i)
      ir(i) = ir(kk)
      ir(kk) = iitemp
      CALL SSWAP(N,A(I,1),MDA,A(KK,1),MDA)
      GO TO 50
C
C           TEST FOR ZERO RANK
C
   70 IF(KRANK.GT.0) GO TO 80
      KRANK=0
      KSURE=0
      RETURN
   80 CONTINUE
C
C        M A I N    L O O P
C
  110 CONTINUE
      J=J+1
      JP1=J+1
      JM1=J-1
      KZ=KRANK
      IF(J.LE.NP) KZ=J
C
C        EACH ROW HAS NN=N-J+1 COMPONENTS
C
      NN=N-J+1
C
C         UB DETERMINES ROW PIVOT
C
  115 IMIN=J
      IF(H(J).EQ.0.) GO TO 170
      RMIN=UB(J)/H(J)
      DO 120 I=J,KZ
      IF(UB(I).GE.H(I)*RMIN) GO TO 120
      RMIN=UB(I)/H(I)
      IMIN=I
  120 CONTINUE
C
C     TEST FOR RANK DEFICIENCY
C
      IF(RMIN.LT.1.0) GO TO 200
      TT=(EB(IMIN)+ABS(DB(IMIN)))/H(IMIN)
      IF(TT.GE.1.0) GO TO 170
C     COMPUTE EXACT UB
      DO 125 I=1,JM1
      W(I)=A(IMIN,I)
  125 CONTINUE
      L=JM1
  130 W(L)=W(L)/A(L,L)
      IF(L.EQ.1) GO TO 150
      LM1=L-1
      DO 140 I=L,JM1
      W(LM1)=W(LM1)-A(I,LM1)*W(I)
  140 CONTINUE
      L=LM1
      GO TO 130
  150 TT=EB(IMIN)
      DO 160 I=1,JM1
      TT=TT+ABS(W(I))*EB(I)
  160 CONTINUE
      UB(IMIN)=TT
      IF(UB(IMIN)/H(IMIN).GE.1.0) GO TO 170
      GO TO 200
C
C        MATRIX REDUCTION
C
  170 CONTINUE
      KK=KRANK
      KRANK=KRANK-1
      KZ=KRANK
      IF(MODE.EQ.0) RETURN
      IF(J.GT.NP) GO TO 172
      CALL XERROR( 'FIRST NP ROWS ARE LINEARLY DEPENDENT',36,
     18,0)
      KRANK=J-1
      RETURN
  172 CONTINUE
      IF(IMIN.GT.KRANK) GO TO 180
*      CALL SSWAP(1,IR(IMIN),1,IR(KK),1)
      iitemp = ir(imin)
      ir(imin) = ir(kk)
      ir(kk) = iitemp
      CALL SSWAP(N,A(IMIN,1),MDA,A(KK,1),MDA)
      CALL SSWAP(1,EB(IMIN),1,EB(KK),1)
      CALL SSWAP(1,UB(IMIN),1,UB(KK),1)
      CALL SSWAP(1,DB(IMIN),1,DB(KK),1)
      CALL SSWAP(1,W(IMIN),1,W(KK),1)
      CALL SSWAP(1,H(IMIN),1,H(KK),1)
  180 IF(J.GT.KRANK) GO TO 300
      GO TO 115
C
C        ROW PIVOT
C
  200 IF(IMIN.EQ.J) GO TO 230
      CALL SSWAP(1,H(J),1,H(IMIN),1)
      CALL SSWAP(N,A(J,1),MDA,A(IMIN,1),MDA)
      CALL SSWAP(1,EB(J),1,EB(IMIN),1)
      CALL SSWAP(1,UB(J),1,UB(IMIN),1)
      CALL SSWAP(1,DB(J),1,DB(IMIN),1)
      CALL SSWAP(1,W(J),1,W(IMIN),1)
*      CALL SSWAP(1,IR(J),1,IR(IMIN),1)
       iitemp = ir(j)
       ir(j) = ir(imin)
       ir(imin) = iitemp
C
C        COLUMN PIVOT
C
  230 CONTINUE
      JMAX=ISAMAX(NN,A(J,J),MDA)
      JMAX=JMAX+J-1
      IF(JMAX.EQ.J) GO TO 240
      CALL SSWAP(M,A(1,J),1,A(1,JMAX),1)
*      CALL SSWAP(1,IC(J),1,IC(JMAX),1)
      iitemp = ic(j)
      ic(j) = ic(jmax)
      ic(jmax) = iitemp
  240 CONTINUE
C
C     APPLY HOUSEHOLDER TRANSFORMATION
C
      TN=SNRM2(NN,A(J,J),MDA)
      IF(TN.EQ.0.0) GO TO 170
      IF(A(J,J).NE.0.0) TN=SIGN(TN,A(J,J))
      CALL SSCAL(NN,1.0/TN,A(J,J),MDA)
      A(J,J)=A(J,J)+1.0
      IF(J.EQ.M) GO TO 250
      DO 248 I=JP1,M
      BB=-SDOT(NN,A(J,J),MDA,A(I,J),MDA)/A(J,J)
      CALL SAXPY(NN,BB,A(J,J),MDA,A(I,J),MDA)
      IF(I.LE.NP) GO TO 248
      IF(H(I).EQ.0.0) GO TO 248
      TT=1.0-(ABS(A(I,J))/H(I))**2
      TT=AMAX1(TT,0.0)
      T=TT
      TT=1.0+.05*TT*(H(I)/W(I))**2
      IF(TT.EQ.1.0) GO TO 244
      H(I)=H(I)*SQRT(T)
      GO TO 246
  244 CONTINUE
      H(I)=SNRM2(N-J,A(I,J+1),MDA)
      W(I)=H(I)
  246 CONTINUE
  248 CONTINUE
  250 CONTINUE
      H(J)=A(J,J)
      A(J,J)=-TN
C
C
C          UPDATE UB, DB
C
      UB(J)=UB(J)/ABS(A(J,J))
      DB(J)=(SIGN(EB(J),DB(J))+DB(J))/A(J,J)
      IF(J.EQ.KRANK) GO TO 300
      DO 260 I=JP1,KRANK
      UB(I)=UB(I)+ABS(A(I,J))*UB(J)
      DB(I)=DB(I)-A(I,J)*DB(J)
  260 CONTINUE
      GO TO 110
C
C        E N D    M A I N    L O O P
C
  300 CONTINUE
C
C        COMPUTE KSURE
C
      KM1=KRANK-1
      DO 318 I=1,KM1
      IS=0
      KMI=KRANK-I
      DO 315 II=1,KMI
      IF(UB(II).LE.UB(II+1)) GO TO 315
      IS=1
      TEMP=UB(II)
      UB(II)=UB(II+1)
      UB(II+1)=TEMP
  315 CONTINUE
      IF(IS.EQ.0) GO TO 320
  318 CONTINUE
  320 CONTINUE
      KSURE=0
      SUM=0.0
      DO 328 I=1,KRANK
      R2=UB(I)*UB(I)
      IF(R2+SUM.GE.1.0) GO TO 330
      SUM=SUM+R2
      KSURE=KSURE+1
  328 CONTINUE
  330 CONTINUE
C
C     IF SYSTEM IS OF REDUCED RANK AND MODE = 2
C     COMPLETE THE DECOMPOSITION FOR SHORTEST LEAST SQUARES SOLUTION
C
      IF(KRANK.EQ.M .OR. MODE.LT.2) GO TO 360
      MMK=M-KRANK
      KP1=KRANK+1
      I=KRANK
  340 TN=SNRM2(MMK,A(KP1,I),1)/A(I,I)
      TN=A(I,I)*SQRT(1.0+TN*TN)
      CALL SSCAL(MMK,1.0/TN,A(KP1,I),1)
      W(I)=A(I,I)/TN+1.0
      A(I,I)=-TN
      IF(I.EQ.1) GO TO 350
      IM1=I-1
      DO 345 II=1,IM1
      TT=-SDOT(MMK,A(KP1,II),1,A(KP1,I),1)/W(I)
      TT=TT-A(I,II)
      CALL SAXPY(MMK,TT,A(KP1,I),1,A(KP1,II),1)
      A(I,II)=A(I,II)+TT*W(I)
  345 CONTINUE
      I=I-1
      GO TO 340
  350 CONTINUE
  360 CONTINUE
      RETURN
      END
