      FUNCTION POCH(A,X)
C***BEGIN PROLOGUE  POCH
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C1,C7A
C***KEYWORDS  POCHHAMMER,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Evaluates a generalization of Pochhammer's symbol.
C***DESCRIPTION
C
C Evaluate a generalization of Pochhammer's symbol
C (A)-sub-X = GAMMA(A+X)/GAMMA(A).  For X a non-negative integer,
C POCH(A,X) is just Pochhammer's symbol.  A and X are single precision.
C This is a preliminary version.  Error handling when POCH(A,X) is
C less than half precision is probably incorrect.  Grossly incorrect
C arguments are not handled properly.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALGAMS,ALNREL,FAC,GAMMA,GAMR,R9LGMC,XERROR
C***END PROLOGUE  POCH
      EXTERNAL GAMMA
      DATA PI / 3.1415926535 89793238 E0 /
C***FIRST EXECUTABLE STATEMENT  POCH
      AX = A + X
      IF (AX.GT.0.0) GO TO 30
      IF (AINT(AX).NE.AX) GO TO 30
C
      IF (A.GT.0.0 .OR. AINT(A).NE.A) CALL XERROR ( 'POCH    A+X IS NON-
     1POSITIVE INTEGER BUT A IS NOT', 48, 2, 2)
C
C WE KNOW HERE THAT BOTH A+X AND A ARE NON-POSITIVE INTEGERS.
C
      POCH = 1.0
      IF (X.EQ.0.0) RETURN
C
      N = X
      IF (AMIN1(A+X,A).LT.(-20.0)) GO TO 20
C
      POCH = (-1.0)**N * FAC(-INT(A))/FAC(-INT(A)-N)
      RETURN
C
 20   POCH = (-1.0)**N * EXP ((A-0.5)*ALNREL(X/(A-1.0))
     1  + X*ALOG(-A+1.0-X) - X + R9LGMC(-A+1.) - R9LGMC(-A-X+1.) )
      RETURN
C
C HERE WE KNOW A+X IS NOT ZERO OR A NEGATIVE INTEGER.
C
 30   POCH = 0.0
      IF (A.LE.0.0 .AND. AINT(A).EQ.A) RETURN
C
      N = ABS(X)
      IF (FLOAT(N).NE.X .OR. N.GT.20) GO TO 50
C
C X IS A SMALL NON-POSITIVE INTEGER, PRESUMMABLY A COMMON CASE.
C
      POCH = 1.0
      IF (N.EQ.0) RETURN
      DO 40 I=1,N
        POCH = POCH * (A+FLOAT(I-1))
 40   CONTINUE
      RETURN
C
 50   ABSAX = ABS(A+X)
      ABSA = ABS(A)
      IF (AMAX1(ABSAX,ABSA).GT.20.0) GO TO 60
      POCH = GAMMA(A+X)*GAMR(A)
      RETURN
C
 60   IF (ABS(X).GT.0.5*ABSA) GO TO 70
C
C HERE ABS(X) IS SMALL AND BOTH ABS(A+X) AND ABS(A) ARE LARGE.  THUS,
C A+X AND A MUST HAVE THE SAME SIGN.  FOR NEGATIVE A, WE USE
C GAMMA(A+X)/GAMMA(A) = GAMMA(-A+1)/GAMMA(-A-X+1) *
C SIN(PI*A)/SIN(PI*(A+X))
C
      B = A
      IF (B.LT.0.0) B = -A - X + 1.0
      POCH = EXP ((B-0.5)*ALNREL(X/B) + X*ALOG(B+X) - X +
     1  R9LGMC(B+X) - R9LGMC(B) )
      IF (A.LT.0.0 .AND. POCH.NE.0.0) POCH = POCH/(COS(PI*X) +
     1  COT(PI*A)*SIN(PI*X))
      RETURN
C
 70   CALL ALGAMS (A+X, ALNGAX, SGNGAX)
      CALL ALGAMS (A, ALNGA, SGNGA)
      POCH = SGNGAX * SGNGA * EXP(ALNGAX-ALNGA)
C
      RETURN
      END
