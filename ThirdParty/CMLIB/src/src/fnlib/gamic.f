      FUNCTION GAMIC(A,X)
C***BEGIN PROLOGUE  GAMIC
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION,GAMMA FUNCTION,
C             SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the complementary incomplete Gamma function.
C***DESCRIPTION
C
C Evaluate the complementary incomplete gamma function
C
C GAMIC = integral from T = X to infinity of EXP(-T) * T**(A-1.)  .
C
C GAMIC is evaluated for arbitrary real values of A and for non-negative
C values of X (even though GAMIC is defined for X .LT. 0.0), except that
C for X = 0 and A .LE. 0.0, GAMIC is undefined.  GAMIC, A, and X
C are single precision.
C
C      A slight deterioration of 2 or 3 digits accuracy will occur when
C GAMIC is very large or very small in absolute value, because log-
C arithmic variables are used.  Also, if the parameter A is very close
C to a negative integer (but not a negative integer), there is a loss
C of accuracy, which is reported if the result is less than half
C machine precision.
C
C Ref. -- W. Gautschi, An Evaluation Procedure for Incomplete Gamma
C Functions, ACM Trans. Math. Software.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALGAMS,ALNGAM,R1MACH,R9GMIC,R9GMIT,R9LGIC,R9LGIT,
C                    XERCLR,XERROR
C***END PROLOGUE  GAMIC
      DATA EPS, SQEPS, ALNEPS, BOT / 4*0.0 /
C***FIRST EXECUTABLE STATEMENT  GAMIC
      IF (EPS.NE.0.) GO TO 10
      EPS = 0.5*R1MACH(3)
      SQEPS = SQRT(R1MACH(4))
      ALNEPS = -ALOG(R1MACH(3))
      BOT = ALOG(R1MACH(1))
C
 10   IF (X.LT.0.0) CALL XERROR ( 'GAMIC   X IS NEGATIVE', 21, 2, 2)
C
      IF (X.GT.0.0) GO TO 20
      IF (A.LE.0.0) CALL XERROR (  'GAMIC   X = 0 AND A LE 0 SO GAMIC IS
     1 UNDEFINED', 46, 3, 2)
C
      GAMIC = EXP (ALNGAM(A+1.0) - ALOG(A))
      RETURN
C
 20   ALX = ALOG(X)
      SGA = 1.0
      IF (A.NE.0.0) SGA = SIGN (1.0, A)
      MA = A + 0.5*SGA
      AEPS = A - FLOAT(MA)
C
      IZERO = 0
      IF (X.GE.1.0) GO TO 60
C
      IF (A.GT.0.5 .OR. ABS(AEPS).GT.0.001) GO TO 50
      FM = -MA
      E = 2.0
      IF (FM.GT.1.0) E = 2.0*(FM+2.0)/(FM*FM-1.0)
      E = E - ALX*X**(-0.001)
      IF (E*ABS(AEPS).GT.EPS) GO TO 50
C
      GAMIC = R9GMIC (A, X, ALX)
      RETURN
C
 50   CALL ALGAMS (A+1.0, ALGAP1, SGNGAM)
      GSTAR = R9GMIT (A, X, ALGAP1, SGNGAM, ALX)
      IF (GSTAR.EQ.0.0) IZERO = 1
      IF (GSTAR.NE.0.0) ALNGS = ALOG (ABS(GSTAR))
      IF (GSTAR.NE.0.0) SGNGS = SIGN (1.0, GSTAR)
      GO TO 70
C
 60   IF (A.LT.X) GAMIC = EXP (R9LGIC(A, X, ALX))
      IF (A.LT.X) RETURN
C
      SGNGAM = 1.0
      ALGAP1 = ALNGAM (A+1.0)
      SGNGS = 1.0
      ALNGS = R9LGIT (A, X, ALGAP1)
C
C EVALUATION OF GAMIC(A,X) IN TERMS OF TRICOMI-S INCOMPLETE GAMMA FN.
C
 70   H = 1.0
      IF (IZERO.EQ.1) GO TO 80
C
      T = A*ALX + ALNGS
      IF (T.GT.ALNEPS) GO TO 90
      IF (T.GT.(-ALNEPS)) H = 1.0 - SGNGS*EXP(T)
C
      IF (ABS(H).LT.SQEPS) CALL XERCLR
      IF (ABS(H).LT.SQEPS) CALL XERROR ( 'GAMIC   RESULT LT HALF PRECISI
     1ON', 32, 1, 1)
C
 80   SGNG = SIGN (1.0, H) * SGA * SGNGAM
      T = ALOG(ABS(H)) + ALGAP1 - ALOG(ABS(A))
      IF (T.LT.BOT) CALL XERCLR
      GAMIC = SGNG * EXP(T)
      RETURN
C
 90   SGNG = -SGNGS * SGA * SGNGAM
      T = T + ALGAP1 - ALOG(ABS(A))
      IF (T.LT.BOT) CALL XERCLR
      GAMIC = SGNG * EXP(T)
      RETURN
C
      END
