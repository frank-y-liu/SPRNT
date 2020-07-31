      REAL FUNCTION WHT(S, INTRPS, M, K, MODOFM, D, MAXRDM, MOMPRD)
C***  SUBROUTINE TO CALCULATE WEIGHT FOR PARTITION M
C
      INTEGER S, M(S), K(S), D, MAXRDM, MI, KI, M1, K1, MODOFM
      REAL INTRPS(S), ZERO, MOMPRD(MAXRDM,MAXRDM)
      ZERO = 0
      DO 10 I=1,S
        INTRPS(I) = ZERO
        K(I) = 0
   10 CONTINUE
      M1 = M(1) + 1
      K1 = D - MODOFM + M1
   20 INTRPS(1) = MOMPRD(M1,K1)
      IF (S.EQ.1) GO TO 40
      DO 30 I=2,S
        MI = M(I) + 1
        KI = K(I) + MI
        INTRPS(I) = INTRPS(I) + MOMPRD(MI,KI)*INTRPS(I-1)
        INTRPS(I-1) = ZERO
        K1 = K1 - 1
        K(I) = K(I) + 1
        IF (K1.GE.M1) GO TO 20
        K1 = K1 + K(I)
        K(I) = 0
   30 CONTINUE
   40 WHT = INTRPS(S)
      RETURN
      END
