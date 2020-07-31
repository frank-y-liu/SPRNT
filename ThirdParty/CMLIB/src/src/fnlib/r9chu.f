      FUNCTION R9CHU(A,B,Z)
C***BEGIN PROLOGUE  R9CHU
C***DATE WRITTEN   770801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C11
C***KEYWORDS  CONFLUENT,HYPERGEOMETRIC,
C             LOGARITHM CONFLUENT HYPERGEOMETRIC FUNCTION,
C             SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Evaluates for large Z  Z**A * CHU(A,B,Z).
C***DESCRIPTION
C
C Evaluate for large Z  Z**A * U(A,B,Z)  where U is the logarithmic
C confluent hypergeometric function.  A rational approximation due to Y.
C L. Luke is used.  When U is not in the asymptotic region, i.e., when A
C or B is large compared with Z, considerable significance loss occurs.
C A warning is provided when the computed result is less than half
C precision.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH,XERROR
C***END PROLOGUE  R9CHU
      DIMENSION AA(4), BB(4)
      DATA EPS, SQEPS / 2*0.0 /
C***FIRST EXECUTABLE STATEMENT  R9CHU
      IF (EPS.NE.0.0) GO TO 10
      EPS = 4.0*R1MACH(4)
      SQEPS = SQRT (R1MACH(4))
C
 10   BP = 1.0 + A - B
      AB = A*BP
      CT2 = 2.0*(Z-AB)
      SAB = A + BP
C
      BB(1) = 1.0
      AA(1) = 1.0
C
      CT3 = SAB + 1.0 + AB
      BB(2) = 1.0 + 2.0*Z/CT3
      AA(2) = 1.0 + CT2/CT3
C
      ANBN = CT3 + SAB + 3.0
      CT1 = 1.0 + 2.0*Z/ANBN
      BB(3) = 1.0 + 6.0*CT1*Z/CT3
      AA(3) = 1.0 + 6.0*AB/ANBN + 3.0*CT1*CT2/CT3
C
      DO 30 I=4,300
        X2I1 = 2*I - 3
        CT1 = X2I1/(X2I1-2.0)
        ANBN = ANBN + X2I1 + SAB
        CT2 = (X2I1 - 1.0) / ANBN
        C2 = X2I1*CT2 - 1.0
        D1Z = X2I1*2.0*Z/ANBN
C
        CT3 = SAB*CT2
        G1 = D1Z + CT1*(C2+CT3)
        G2 = D1Z - C2
        G3 = CT1*(1.0 - CT3 - 2.0*CT2)
C
        BB(4) = G1*BB(3) + G2*BB(2) + G3*BB(1)
        AA(4) = G1*AA(3) + G2*AA(2) + G3*AA(1)
        IF (ABS(AA(4)*BB(1)-AA(1)*BB(4)).LT.EPS*ABS(BB(4)*BB(1)))
     1    GO TO 40
C
C IF OVERFLOWS OR UNDERFLOWS PROVE TO BE A PROBLEM, THE STATEMENTS
C BELOW COULD BE ALTERED TO INCORPORATE A DYNAMICALLY ADJUSTED SCALE
C FACTOR.
C
        DO 20 J=1,3
          BB(J) = BB(J+1)
          AA(J) = AA(J+1)
 20     CONTINUE
 30   CONTINUE
      CALL XERROR ( 'R9CHU   NO CONVERGENCE IN 300 TERMS', 35, 1, 2)
C
 40   R9CHU = AA(4)/BB(4)
C
      IF (R9CHU.LT.SQEPS .OR. R9CHU.GT.1.0/SQEPS) CALL XERROR (   'R9CHU
     1   ANSWER LESS THAN HALF PRECISION', 39, 2, 1)
C
      RETURN
      END
