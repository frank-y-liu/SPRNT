      SUBROUTINE SDNTP (H,K,N,NQ,T,TOUT,YH,Y)
C***BEGIN PROLOGUE  SDNTP
C***REFER TO  SDRIV3
C   Subroutine SDNTP interpolates the K-th derivative of Y at TOUT,
C   using the data in the YH array.  If K has a value greater than NQ,
C   the NQ-th derivative is calculated.
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  SDNTP
      REAL FACTOR, H, R, T, TOUT, Y(*), YH(N,*)
C***FIRST EXECUTABLE STATEMENT  SDNTP
      KUSED = MIN(K, NQ)
      IF (KUSED .EQ. 0) THEN
        DO 10 I = 1,N
 10       Y(I) = YH(I,NQ+1)
        R = ((TOUT - T)/H)
        DO 20 JJ = 1,NQ
          J = NQ + 1 - JJ
          DO 20 I = 1,N
 20         Y(I) = YH(I,J) + R*Y(I)
      ELSE
        FACTOR = 1.E0
        DO 40 KK = 1,KUSED
 40       FACTOR = FACTOR*REAL(NQ+1-KK)
        DO 50 I = 1,N
 50       Y(I) = FACTOR*YH(I,NQ+1)
        IF (KUSED .NE. NQ) THEN
          R = ((TOUT - T)/H)
          DO 80 JJ = KUSED+1,NQ
            J = K + 1 + NQ - JJ
            FACTOR = 1.E0
            DO 60 KK = 1,KUSED
 60           FACTOR = FACTOR*REAL(J-KK)
            DO 70 I = 1,N
 70           Y(I) = FACTOR*YH(I,J) + R*Y(I)
 80         CONTINUE
        END IF
        DO 100 I = 1,N
 100      Y(I) = Y(I)*H**(-KUSED)
      END IF
      END
