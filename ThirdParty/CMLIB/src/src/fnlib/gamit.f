      FUNCTION GAMIT(A,X)
C***BEGIN PROLOGUE  GAMIT
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  GAMMA FUNCTION,INCOMPLETE GAMMA FUNCTION,SPECIAL FUNCTION,
C             TRICOMI-S
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes Tricomi's form of the incomplete Gamma function.
C***DESCRIPTION
C
C Evaluate Tricomi's incomplete gamma function defined by
C
C GAMIT = X**(-A)/GAMMA(A) * integral T = 0 to X of EXP(-T) * T**(A-1.)
C
C for A .GT. 0.0 and by analytic
C continuation for A .LE. 0.0.  GAMMA(X) is the complete
C gamma function of X.  GAMIT is evaluated for arbitrary real values of
C A and for non-negative values of X (even though GAMIT is defined for
C X .LT. 0.0), except that for X = 0 and A .LE. 0.0, GAMIT is infinite,
C a fatal error.
C
C      A slight deterioration of 2 or 3 digits accuracy will occur when
C GAMIT is very large or very small in absolute value, because log-
C arithmic variables are used.  Also, if the parameter A is very close
C to a negative integer (but not a negative integer), there is a loss
C of accuracy, which is reported if the result is less than half
C machine precision.
C
C Ref. -- W. Gautschi, An Evaluation Procedure for Incomplete Gamma
C Functions, ACM Trans. Math. Software, Vol. 5, No. 4, December 1979.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALGAMS,ALNGAM,GAMR,R1MACH,R9GMIT,R9LGIC,R9LGIT,
C                    XERCLR,XERROR
C***END PROLOGUE  GAMIT
      DATA ALNEPS, SQEPS, BOT / 3*0.0 /
C***FIRST EXECUTABLE STATEMENT  GAMIT
      IF (ALNEPS.NE.0.0) GO TO 10
      ALNEPS = -ALOG(R1MACH(3))
      SQEPS = SQRT(R1MACH(4))
      BOT = ALOG(R1MACH(1))
C
 10   IF (X.LT.0.0) CALL XERROR ( 'GAMIT   X IS NEGATIVE', 21, 2, 2)
C
      IF (X.NE.0.0) ALX = ALOG(X)
      SGA = 1.0
      IF (A.NE.0.0) SGA = SIGN (1.0, A)
      AINTA = AINT (A+0.5*SGA)
      AEPS = A - AINTA
C
      IF (X.GT.0.0) GO TO 20
      GAMIT = 0.0
      IF (AINTA.GT.0.0 .OR. AEPS.NE.0.0) GAMIT = GAMR(A+1.0)
      RETURN
C
 20   IF (X.GT.1.0) GO TO 40
      IF (A.GE.(-0.5) .OR. AEPS.NE.0.0) CALL ALGAMS (A+1.0, ALGAP1,
     1  SGNGAM)
      GAMIT = R9GMIT (A, X, ALGAP1, SGNGAM, ALX)
      RETURN
C
 40   IF (A.LT.X) GO TO 50
      T = R9LGIT (A, X, ALNGAM(A+1.0))
      IF (T.LT.BOT) CALL XERCLR
      GAMIT = EXP(T)
      RETURN
C
 50   ALNG = R9LGIC (A, X, ALX)
C
C EVALUATE GAMIT IN TERMS OF ALOG(GAMIC(A,X))
C
      H = 1.0
      IF (AEPS.EQ.0.0 .AND. AINTA.LE.0.0) GO TO 60
      CALL ALGAMS (A+1.0, ALGAP1, SGNGAM)
      T = ALOG(ABS(A)) + ALNG - ALGAP1
      IF (T.GT.ALNEPS) GO TO 70
      IF (T.GT.(-ALNEPS)) H = 1.0 - SGA*SGNGAM*EXP(T)
      IF (ABS(H).GT.SQEPS) GO TO 60
      CALL XERCLR
      CALL XERROR ( 'GAMIT   RESULT LT HALF PRECISION', 32, 1, 1)
C
 60   T = -A*ALX + ALOG(ABS(H))
      IF (T.LT.BOT) CALL XERCLR
      GAMIT = SIGN (EXP(T), H)
      RETURN
C
 70   T = T - A*ALX
      IF (T.LT.BOT) CALL XERCLR
      GAMIT = -SGA*SGNGAM*EXP(T)
      RETURN
C
      END
