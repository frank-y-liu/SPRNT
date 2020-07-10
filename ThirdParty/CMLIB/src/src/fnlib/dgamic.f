      DOUBLE PRECISION FUNCTION DGAMIC(A,X)
C***BEGIN PROLOGUE  DGAMIC
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  COMPLEMENTARY,COMPLEMENTARY INCOMPLETE GAMMA FUNCTION,
C             DOUBLE PRECISION,GAMMA,GAMMA FUNCTION,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Calculates the d.p. complementary incomplete Gamma function
C***DESCRIPTION
C
C Evaluate the complementary incomplete gamma function
C
C DGAMIC = integral from T = X to infinity of EXP(-T) * T**(A-1.)  .
C
C DGAMIC is evaluated for arbitrary real values of A and for non-
C negative values of X (even though DGAMIC is defined for X .LT. 0.0),
C except that for X = 0 and A .LE. 0.0, DGAMIC is undefined.  The
C function and both arguments are double precision.
C
C      A slight deterioration of 2 or 3 digits accuracy will occur when
C DGAMIC is very large or very small in absolute value, because log-
C arithmic variables are used.  Also, if the parameter A is very close
C to a negative INTEGER (but not a negative integer), there is a loss
C of accuracy, which is reported if the result is less than half
C machine precision.
C
C Ref. -- W. Gautschi, An Evaluation Procedure for Incomplete Gamma
C Functions, ACM Trans. Math. Software.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,D9GMIC,D9GMIT,D9LGIC,D9LGIT,DINT,DLGAMS,
C                    DLNGAM,XERCLR,XERROR
C***END PROLOGUE  DGAMIC
      DOUBLE PRECISION A, X, AEPS, AINTA, ALGAP1, ALNEPS, ALNGS, ALX,
     1  BOT, E, EPS, GSTAR, H, SGA, SGNG, SGNGAM, SGNGS, SQEPS, T,
     2  D1MACH, DLNGAM, D9GMIC, D9GMIT, D9LGIC, D9LGIT, DINT
      DATA EPS, SQEPS, ALNEPS, BOT / 4*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DGAMIC
      IF (EPS.NE.0.D0) GO TO 10
      EPS = 0.5D0*D1MACH(3)
      SQEPS = DSQRT (D1MACH(4))
      ALNEPS = -DLOG (D1MACH(3))
      BOT = DLOG (D1MACH(1))
C
 10   IF (X.LT.0.D0) CALL XERROR ( 'DGAMIC  X IS NEGATIVE', 21, 2, 2)
C
      IF (X.GT.0.D0) GO TO 20
      IF (A.LE.0.D0) CALL XERROR ( 'DGAMIC  X = 0 AND A LE 0 SO DGAMIC I
     1S UNDEFINED', 47, 3, 2)
C
      DGAMIC = DEXP (DLNGAM(A+1.D0) - DLOG(A))
      RETURN
C
 20   ALX = DLOG (X)
      SGA = 1.0D0
      IF (A.NE.0.D0) SGA = DSIGN (1.0D0, A)
      AINTA = DINT (A + 0.5D0*SGA)
      AEPS = A - AINTA
C
      IZERO = 0
      IF (X.GE.1.0D0) GO TO 40
C
      IF (A.GT.0.5D0 .OR. DABS(AEPS).GT.0.001D0) GO TO 30
      E = 2.0D0
      IF (-AINTA.GT.1.D0) E = 2.D0*(-AINTA+2.D0)/(AINTA*AINTA-1.0D0)
      E = E - ALX * X**(-0.001D0)
      IF (E*DABS(AEPS).GT.EPS) GO TO 30
C
      DGAMIC = D9GMIC (A, X, ALX)
      RETURN
C
 30   CALL DLGAMS (A+1.0D0, ALGAP1, SGNGAM)
      GSTAR = D9GMIT (A, X, ALGAP1, SGNGAM, ALX)
      IF (GSTAR.EQ.0.D0) IZERO = 1
      IF (GSTAR.NE.0.D0) ALNGS = DLOG (DABS(GSTAR))
      IF (GSTAR.NE.0.D0) SGNGS = DSIGN (1.0D0, GSTAR)
      GO TO 50
C
 40   IF (A.LT.X) DGAMIC = DEXP (D9LGIC(A, X, ALX))
      IF (A.LT.X) RETURN
C
      SGNGAM = 1.0D0
      ALGAP1 = DLNGAM (A+1.0D0)
      SGNGS = 1.0D0
      ALNGS = D9LGIT (A, X, ALGAP1)
C
C EVALUATION OF DGAMIC(A,X) IN TERMS OF TRICOMI-S INCOMPLETE GAMMA FN.
C
 50   H = 1.D0
      IF (IZERO.EQ.1) GO TO 60
C
      T = A*ALX + ALNGS
      IF (T.GT.ALNEPS) GO TO 70
      IF (T.GT.(-ALNEPS)) H = 1.0D0 - SGNGS*DEXP(T)
C
      IF (DABS(H).LT.SQEPS) CALL XERCLR
      IF (DABS(H).LT.SQEPS) CALL XERROR ( 'DGAMIC  RESULT LT HALF PRECIS
     1ION', 32, 1, 1)
C
 60   SGNG = DSIGN (1.0D0, H) * SGA * SGNGAM
      T = DLOG(DABS(H)) + ALGAP1 - DLOG(DABS(A))
      IF (T.LT.BOT) CALL XERCLR
      DGAMIC = SGNG * DEXP(T)
      RETURN
C
 70   SGNG = -SGNGS * SGA * SGNGAM
      T = T + ALGAP1 - DLOG(DABS(A))
      IF (T.LT.BOT) CALL XERCLR
      DGAMIC = SGNG * DEXP(T)
      RETURN
C
      END
