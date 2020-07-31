      FUNCTION BETAI(X,PIN,QIN)
C***BEGIN PROLOGUE  BETAI
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7F
C***KEYWORDS  BETA FUNCTION,INCOMPLETE BETA FUNCTION,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the incomplete Beta function.
C***DESCRIPTION
C
C BETAI(X,PIN,QIN) calculates the incomplete beta function
C ratio based on Bosten and Battiste, Remark on Algorithm 179,
C Comm. ACM, Vol. 17, p. 153, (1974).
C
C X     value to which function is to be integrated.  X must be in (0,1)
C PIN   input (1st) parameter (must be greater than 0)
C QIN   input (2nd) parameter (must be greater than 0)
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALBETA,R1MACH,XERROR
C***END PROLOGUE  BETAI
      DATA EPS, ALNEPS, SML, ALNSML / 4*0.0 /
C***FIRST EXECUTABLE STATEMENT  BETAI
      IF (EPS.NE.0.) GO TO 10
      EPS = R1MACH(3)
      ALNEPS = ALOG(EPS)
      SML = R1MACH(1)
      ALNSML = ALOG(SML)
C
 10   IF (X.LT.0. .OR. X.GT.1.0) CALL XERROR (  'BETAI   X IS NOT IN THE
     1 RANGE (0,1)', 35, 1, 2)
      IF (PIN.LE.0. .OR. QIN.LE.0.) CALL XERROR ( 'BETAI   P AND/OR Q IS
     1 LE ZERO', 29, 2, 2)
C
      Y = X
      P = PIN
      Q = QIN
      IF (Q.LE.P .AND. X.LT.0.8) GO TO 20
      IF (X.LT.0.2) GO TO 20
      Y = 1.0 - Y
      P = QIN
      Q = PIN
C
 20   IF ((P+Q)*Y/(P+1.).LT.EPS) GO TO 80
C
C EVALUATE THE INFINITE SUM FIRST.
C TERM WILL EQUAL Y**P/BETA(PS,P) * (1.-PS)I * Y**I / FAC(I)
C
      PS = Q - AINT(Q)
      IF (PS.EQ.0.) PS = 1.0
      XB = P*ALOG(Y) -  ALBETA(PS, P) - ALOG(P)
      BETAI = 0.0
      IF (XB.LT.ALNSML) GO TO 40
C
      BETAI = EXP (XB)
      TERM = BETAI*P
      IF (PS.EQ.1.0) GO TO 40
C
      N = AMAX1 (ALNEPS/ALOG(Y), 4.0)
      DO 30 I=1,N
        TERM = TERM*(FLOAT(I)-PS)*Y/FLOAT(I)
        BETAI = BETAI + TERM/(P+FLOAT(I))
 30   CONTINUE
C
C NOW EVALUATE THE FINITE SUM, MAYBE.
C
 40   IF (Q.LE.1.0) GO TO 70
C
      XB = P*ALOG(Y) + Q*ALOG(1.0-Y) - ALBETA(P,Q) - ALOG(Q)
      IB = AMAX1 (XB/ALNSML, 0.0)
      TERM = EXP (XB - FLOAT(IB)*ALNSML)
      C = 1.0/(1.0-Y)
      P1 = Q*C/(P+Q-1.)
C
      FINSUM = 0.0
      N = Q
      IF (Q.EQ.FLOAT(N)) N = N - 1
      DO 50 I=1,N
        IF (P1.LE.1.0 .AND. TERM/EPS.LE.FINSUM) GO TO 60
        TERM = (Q-FLOAT(I-1))*C*TERM/(P+Q-FLOAT(I))
C
        IF (TERM.GT.1.0) IB = IB - 1
        IF (TERM.GT.1.0) TERM = TERM*SML
C
        IF (IB.EQ.0) FINSUM = FINSUM + TERM
 50   CONTINUE
C
 60   BETAI = BETAI + FINSUM
 70   IF (Y.NE.X .OR. P.NE.PIN) BETAI = 1.0 - BETAI
      BETAI = AMAX1 (AMIN1 (BETAI, 1.0), 0.0)
      RETURN
C
 80   BETAI = 0.0
      XB = P*ALOG(AMAX1(Y,SML)) - ALOG(P) - ALBETA(P,Q)
      IF (XB.GT.ALNSML .AND. Y.NE.0.) BETAI = EXP (XB)
      IF (Y.NE.X .OR. P.NE.PIN) BETAI = 1.0 - BETAI
      RETURN
C
      END
