      SUBROUTINE NXPRT(PRTCNT, S, M)
C
C***  SUBROUTINE TO COMPUTE THE NEXT S PARTITION
C
      INTEGER S, M(S), PRTCNT, I, MSUM
      IF (PRTCNT.GT.0) GO TO 20
      DO 10 I=1,S
        M(I) = 0
   10 CONTINUE
      PRTCNT = 1
      RETURN
   20 PRTCNT = PRTCNT + 1
      MSUM = M(1)
      IF (S.EQ.1) GO TO 60
      DO 50 I=2,S
        MSUM = MSUM + M(I)
        IF (M(1).LE.M(I)+1) GO TO 40
        M(1) = MSUM - (I-1)*(M(I)+1)
        DO 30 L=2,I
          M(L) = M(I) + 1
   30   CONTINUE
        RETURN
   40   M(I) = 0
   50 CONTINUE
   60 M(1) = MSUM + 1
      RETURN
      END
