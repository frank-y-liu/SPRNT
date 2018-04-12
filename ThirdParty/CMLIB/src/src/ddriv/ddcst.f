      SUBROUTINE DDCST (MAXORD,MINT,ISWFLG,EL,TQ)
C***BEGIN PROLOGUE  DDCST
C***REFER TO  DDRIV3
C  DDCST is called by DDNTL and sets coefficients used by the core
C  integrator DDSTP.  The array EL determines the basic method.
C  The array TQ is involved in adjusting the step size in relation
C  to truncation error.  EL and TQ depend upon MINT, and are calculated
C  for orders 1 to MAXORD(.LE. 12).  For each order NQ, the coefficients
C  EL are calculated from the generating polynomial:
C    L(T) = EL(1,NQ) + EL(2,NQ)*T + ... + EL(NQ+1,NQ)*T**NQ.
C  For the implicit Adams methods, L(T) is given by
C    dL/dT = (1+T)*(2+T)* ... *(NQ-1+T)/K,   L(-1) = 0,
C    where      K = factorial(NQ-1).
C  For the Gear methods,
C    L(T) = (1+T)*(2+T)* ... *(NQ+T)/K,
C    where      K = factorial(NQ)*(1 + 1/2 + ... + 1/NQ).
C  For each order NQ, there are three components of TQ.
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDCST
      DOUBLE PRECISION EL(13,12), FACTRL(12), GAMMA(14), SUM, TQ(3,12)
C***FIRST EXECUTABLE STATEMENT  DDCST
      FACTRL(1) = 1.D0
      IF (MAXORD .GE. 2) THEN
        DO 10 I = 2,MAXORD
 10       FACTRL(I) = DBLE(I)*FACTRL(I-1)
      END IF
C                                             COMPUTE ADAMS COEFFICIENTS
      IF (MINT .EQ. 1) THEN
        GAMMA(1) = 1.D0
        DO 40 I = 1,MAXORD+1
          SUM = 0.D0
          DO 30 J = 1,I
 30         SUM = SUM - GAMMA(J)/DBLE(I-J+2)
 40       GAMMA(I+1) = SUM
        EL(1,1) = 1.D0
        EL(2,1) = 1.D0
        EL(2,2) = 1.D0
        EL(3,2) = 1.D0
        IF (MAXORD .GE. 3) THEN
          DO 60 J = 3,MAXORD
            EL(2,J) = FACTRL(J-1)
            DO 50 I = 3,J
 50           EL(I,J) = DBLE(J-1)*EL(I,J-1) + EL(I-1,J-1)
 60         EL(J+1,J) = 1.D0
        END IF
        IF (MAXORD .GE. 2) THEN
          DO 80 J = 2,MAXORD
            EL(1,J) = EL(1,J-1) + GAMMA(J)
            EL(2,J) = 1.D0
            DO 80 I = 3,J+1
 80           EL(I,J) = EL(I,J)/(DBLE(I-1)*FACTRL(J-1))
        END IF
        DO 100 J = 1,MAXORD
          TQ(1,J) = -1.D0/(FACTRL(J)*GAMMA(J))
          TQ(2,J) = -1.D0/GAMMA(J+1)
 100      TQ(3,J) = -1.D0/GAMMA(J+2)
C                                              COMPUTE GEAR COEFFICIENTS
      ELSE IF (MINT .EQ. 2) THEN
        EL(1,1) = 1.D0
        EL(2,1) = 1.D0
        IF (MAXORD .GE. 2) THEN
          DO 130 J = 2,MAXORD
            EL(1,J) = FACTRL(J)
            DO 120 I = 2,J
 120          EL(I,J) = DBLE(J)*EL(I,J-1) + EL(I-1,J-1)
 130        EL(J+1,J) = 1.D0
          SUM = 1.D0
          DO 150 J = 2,MAXORD
            SUM = SUM + 1.D0/DBLE(J)
            DO 150 I = 1,J+1
 150          EL(I,J) = EL(I,J)/(FACTRL(J)*SUM)
        END IF
        DO 170 J = 1,MAXORD
          IF (J .GT. 1) TQ(1,J) = 1.D0/FACTRL(J-1)
          TQ(2,J) = DBLE(J+1)/EL(1,J)
 170      TQ(3,J) = DBLE(J+2)/EL(1,J)
      END IF
C                          Compute constants used in the stiffness test.
C                          These are the ratio of TQ(2,NQ) for the Gear
C                          methods to those for the Adams methods.
      IF (ISWFLG .EQ. 3) THEN
        MXRD = MIN(MAXORD, 5)
        IF (MINT .EQ. 2) THEN
          GAMMA(1) = 1.D0
          DO 190 I = 1,MXRD
            SUM = 0.D0
            DO 180 J = 1,I
 180          SUM = SUM - GAMMA(J)/DBLE(I-J+2)
 190        GAMMA(I+1) = SUM
        END IF
        IF (MXRD .GE. 2) THEN
          SUM = 1.D0
          DO 200 I = 2,MXRD
            SUM = SUM + 1.D0/DBLE(I)
 200        EL(1+I,1) = -DBLE(I+1)*SUM*GAMMA(I+1)
        END IF
      END IF
      END
