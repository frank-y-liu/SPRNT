      SUBROUTINE FFTBRG(X,Y,TABLE,M,LL,ISN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     FFTBRG IS IN PLACE DFT COMPUTATION USING SANDE ALGORITHM
C        AND MARKEL PRUNING MODIFICATION.
C
C     INPUTS ARE:
C        X     =ARRAY OF LENGTH N. REAL PART OF THE INPUT.
C        Y     =ARRAY OF LENGTH N. IMAGINARY PART OF THE INPUT.
C        TABLE =ARRAY OF LENGTH N/4 +1. QUARTER-LENGTH COS.
C        M     =INTEGER. (N=2**M)
C        LL    =NUMBER OF STAGES IN WHICH NO PRUNING IS ALLOWABLE.
C               THE ACTUAL NUMBER OF DATA POINTS.
C        ISN   =DFT OR IDFT INDICATOR
C              = -1 FOR DFT
C              = +1 FOR IDFT(INVERSE FFT)
C
C     OUTPUTS ARE:
C        X     =REAL PART OF OUTPUT
C        Y     =IMAGINARY PART OF OUTPUT
C
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DIMENSION X(2),Y(2),TABLE(*),L(12)
      EQUIVALENCE (L12,L(1)),(L11,L(2)),(L10,L(3)),(L9,L(4)),(L8,L(5)),
     1 (L7,L(6)),(L6,L(7)),(L5,L(8)),(L4,L(9)),(L3,L(10)),(L2,L(11)),
     2 (L1,L(12))
C
      N=2**M
      ND4=N/4
      ND4P1=ND4+1
      ND4P2=ND4+2
      ND2P2=ND4+ND4P2
      LLL=2**LL
         DO 8 LO=1,M
            LMX=2**(M-LO)
            LMM=LMX
            LIX=2*LMX
            ISCL=N/LIX
C     TEST FOR PRUNING
      IF(LO-M+LL) 1,2,2
    1 LMM=LLL
    2       DO 8 LM=1,LMM
               IARG=(LM-1)*ISCL+1
               IF(IARG.LE.ND4P1) GO TO 4
               K1=ND2P2-IARG
               C=-TABLE(K1)
               K3=IARG-ND4
               S=ISN*TABLE(K3)
               GO TO 6
    4       C=TABLE(IARG)
               K2=ND4P2-IARG
               S=ISN*TABLE(K2)
    6 CONTINUE
               DO 8 LI=LIX,N,LIX
                  J1=LI-LIX+LM
                  J2=J1+LMX
                  T1=X(J1)-X(J2)
                  T2=Y(J1)-Y(J2)
                  X(J1)=X(J1)+X(J2)
                  Y(J1)=Y(J1)+Y(J2)
                  X(J2)=C*T1-S*T2
                  Y(J2)=C*T2+S*T1
    8 CONTINUE
C
C     PERFORM BIT REVERSAL
C
      DO 40 J=1,12
         L(J)=1
         IF (J-M) 31,31,40
   31    L(J)=2**(M+1-J)
   40 CONTINUE
         JN=1
         DO 60 J1=1,L1
         DO 60 J2=J1,L2,L1
         DO 60 J3=J2,L3,L2
         DO 60 J4=J3,L4,L3
         DO 60 J5=J4,L5,L4
         DO 60 J6=J5,L6,L5
         DO 60 J7=J6,L7,L6
         DO 60 J8=J7,L8,L7
         DO 60 J9=J8,L9,L8
         DO 60 J10=J9,L10,L9
         DO 60 J11=J10,L11,L10
         DO 60 JR=J11,L12,L11
            IF (JN-JR) 51,51,54
   51    R=X(JN)
         X(JN)=X(JR)
         X(JR)=R
         FI=Y(JN)
         Y(JN)=Y(JR)
         Y(JR)=FI
   54    IF(ISN) 53,53,52
   52    X(JR)=X(JR)/FLOAT(N)
         Y(JR)=Y(JR)/FLOAT(N)
   53    JN=JN+1
   60 CONTINUE
      RETURN
      END
