      FUNCTION R9GMIT(A,X,ALGAP1,SGNGAM,ALX)
C***BEGIN PROLOGUE  R9GMIT
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  GAMMA FUNCTION,INCOMPLETE GAMMA FUNCTION,SMALL X,
C             SPECIAL FUNCTION,TRICOMI-S
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes Tricomi's incomplete Gamma function for small X.
C***DESCRIPTION
C
C Compute Tricomi's incomplete gamma function for small X.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALNGAM,R1MACH,XERROR
C***END PROLOGUE  R9GMIT
      DATA EPS, BOT / 2*0.0 /
C***FIRST EXECUTABLE STATEMENT  R9GMIT
      IF (EPS.EQ.0.0) EPS = 0.5*R1MACH(3)
      IF (BOT.EQ.0.0) BOT = ALOG(R1MACH(1))
C
      IF (X.LE.0.0) CALL XERROR ( 'R9GMIT  X SHOULD BE GT 0', 24, 1, 2)
C
      MA = A + 0.5
      IF (A.LT.0.0) MA = A - 0.5
      AEPS = A - FLOAT(MA)
C
      AE = A
      IF (A.LT.(-0.5)) AE = AEPS
C
      T = 1.0
      TE = AE
      S = T
      DO 20 K=1,200
        FK = K
        TE = -X*TE/FK
        T = TE/(AE+FK)
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 30
 20   CONTINUE
      CALL XERROR ( 'R9GMIT  NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SER
     1IES', 54, 2, 2)
C
 30   IF (A.GE.(-0.5)) ALGS = -ALGAP1 + ALOG(S)
      IF (A.GE.(-0.5)) GO TO 60
C
      ALGS = -ALNGAM(1.0+AEPS) + ALOG(S)
      S = 1.0
      M = -MA - 1
      IF (M.EQ.0) GO TO 50
      T = 1.0
      DO 40 K=1,M
        T = X*T/(AEPS-FLOAT(M+1-K))
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 50
 40   CONTINUE
C
 50   R9GMIT = 0.0
      ALGS = -FLOAT(MA)*ALOG(X) + ALGS
      IF (S.EQ.0.0 .OR. AEPS.EQ.0.0) GO TO 60
C
      SGNG2 = SGNGAM*SIGN(1.0,S)
      ALG2 = -X - ALGAP1 + ALOG(ABS(S))
C
      IF (ALG2.GT.BOT) R9GMIT = SGNG2*EXP(ALG2)
      IF (ALGS.GT.BOT) R9GMIT = R9GMIT + EXP(ALGS)
      RETURN
C
 60   R9GMIT = EXP(ALGS)
      RETURN
C
      END
