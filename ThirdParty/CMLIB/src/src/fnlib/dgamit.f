      DOUBLE PRECISION FUNCTION DGAMIT(A,X)
C***BEGIN PROLOGUE  DGAMIT
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  COMPLEMENTARY,COMPLEMENTARY INCOMPLETE GAMMA FUNCTION,
C             DOUBLE PRECISION,GAMMA FUNCTION,SPECIAL FUNCTION,TRICOMI
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Calculates Tricomi's form of the incomplete Gamma function.
C***DESCRIPTION
C
C Evaluate Tricomi's incomplete gamma function defined by
C
C DGAMIT = X**(-A)/GAMMA(A) * integral T = 0 to X of EXP(-T) * T**(A-1.)
C
C for A .GT. 0.0 and by analytic
C continuation for A .LE. 0.0.  Gamma(X) is the complete
C gamma function of X.  DGAMIT is evaluated for arbitrary real values of
C A and for non-negative values of X (even though DGAMIT is defined for
C X .LT. 0.0), except that for X = 0 and A .LE. 0.0, DGAMIT is infinite,
C a fatal error.  The function and both arguments are double precision.
C
C      A slight deterioration of 2 or 3 digits accuracy will occur when
C DGAMIT is very large or very small in absolute value, because log-
C arithmic variables are used.  Also, if the parameter A is very close
C to a negative integer (but not a negative integer), there is a loss
C of accuracy, which is reported if the result is less than half
C machine precision.
C
C Ref. -- W. Gautschi, An Evaluation Procedure for Incomplete Gamma
C Functions, ACM Trans. Math. Software, Vol. 5, No. 4, December 1979.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,D9GMIT,D9LGIC,D9LGIT,DGAMR,DINT,DLGAMS,
C                    DLNGAM,XERCLR,XERROR
C***END PROLOGUE  DGAMIT
      DOUBLE PRECISION A, X, AEPS, AINTA, ALGAP1, ALNEPS, ALNG, ALX,
     1  BOT, H, SGA, SGNGAM, SQEPS, T, D1MACH, DGAMR, D9GMIT, D9LGIT,
     2  DLNGAM, D9LGIC, DINT
      DATA ALNEPS, SQEPS, BOT / 3*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DGAMIT
      IF (ALNEPS.NE.0.D0) GO TO 10
      ALNEPS = -DLOG (D1MACH(3))
      SQEPS = DSQRT (D1MACH(4))
      BOT = DLOG (D1MACH(1))
C
 10   IF (X.LT.0.D0) CALL XERROR ( 'DGAMIT  X IS NEGATIVE', 21, 2, 2)
C
      IF (X.NE.0.D0) ALX = DLOG (X)
      SGA = 1.0D0
      IF (A.NE.0.D0) SGA = DSIGN (1.0D0, A)
      AINTA = DINT (A + 0.5D0*SGA)
      AEPS = A - AINTA
C
      IF (X.GT.0.D0) GO TO 20
      DGAMIT = 0.0D0
      IF (AINTA.GT.0.D0 .OR. AEPS.NE.0.D0) DGAMIT = DGAMR(A+1.0D0)
      RETURN
C
 20   IF (X.GT.1.D0) GO TO 30
      IF (A.GE.(-0.5D0) .OR. AEPS.NE.0.D0) CALL DLGAMS (A+1.0D0, ALGAP1,
     1  SGNGAM)
      DGAMIT = D9GMIT (A, X, ALGAP1, SGNGAM, ALX)
      RETURN
C
 30   IF (A.LT.X) GO TO 40
      T = D9LGIT (A, X, DLNGAM(A+1.0D0))
      IF (T.LT.BOT) CALL XERCLR
      DGAMIT = DEXP (T)
      RETURN
C
 40   ALNG = D9LGIC (A, X, ALX)
C
C EVALUATE DGAMIT IN TERMS OF DLOG (DGAMIC (A, X))
C
      H = 1.0D0
      IF (AEPS.EQ.0.D0 .AND. AINTA.LE.0.D0) GO TO 50
C
      CALL DLGAMS (A+1.0D0, ALGAP1, SGNGAM)
      T = DLOG (DABS(A)) + ALNG - ALGAP1
      IF (T.GT.ALNEPS) GO TO 60
C
      IF (T.GT.(-ALNEPS)) H = 1.0D0 - SGA * SGNGAM * DEXP(T)
      IF (DABS(H).GT.SQEPS) GO TO 50
C
      CALL XERCLR
      CALL XERROR ( 'DGAMIT  RESULT LT HALF PRECISION', 32, 1, 1)
C
 50   T = -A*ALX + DLOG(DABS(H))
      IF (T.LT.BOT) CALL XERCLR
      DGAMIT = DSIGN (DEXP(T), H)
      RETURN
C
 60   T = T - A*ALX
      IF (T.LT.BOT) CALL XERCLR
      DGAMIT = -SGA * SGNGAM * DEXP(T)
      RETURN
C
      END
