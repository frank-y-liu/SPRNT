      DOUBLE PRECISION FUNCTION DBETAI(X,PIN,QIN)
C***BEGIN PROLOGUE  DBETAI
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY (YYMMDD)
C   000601 Changed MAX to MAX1        (RFB)
C***CATEGORY NO.  C7F
C***KEYWORDS  BETA,BETA FUNCTION,DOUBLE PRECISION,
C             INCOMPLETE BETA FUNCTION,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Calculates the d.p. incomplete Beta function.
C***DESCRIPTION
C
C DBETAI(X,PIN,QIN) calculates the double precision incomplete beta
C function ratio for double precision arguments X, PIN, and QIN.  It is
C based on Bosten and Battiste, Remark on Algorithm 179, Comm. ACM,
C V 17, P 153, (1974).
C
C             Input Arguments --
C X      upper limit of integration.  X must be in (0,1) inclusive.
C PIN    first beta distribution parameter.  PIN must be .GT. 0.0.
C QIN    second beta distribution parameter.  QIN must be .GT. 0.0.
C DBETAI the incomplete beta function ratio is the probability that a
C        random variable from a beta distribution having parameters
C        PIN and QIN will be less than or equal to X.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DINT,DLBETA,XERROR
C***END PROLOGUE  DBETAI
      DOUBLE PRECISION X, PIN, QIN, ALNEPS, ALNSML, C, EPS, FINSUM, P,
     1  PS, Q, SML, TERM, XB, XI, Y, DINT, D1MACH, DLBETA
      DATA EPS, ALNEPS, SML, ALNSML / 4*0.0D0 /
C***FIRST EXECUTABLE STATEMENT  DBETAI
      IF (EPS.NE.0.0D0) GO TO 10
      EPS = D1MACH(3)
      ALNEPS = DLOG (EPS)
      SML = D1MACH(1)
      ALNSML = DLOG (SML)
C
 10   IF (X.LT.0.D0 .OR. X.GT.1.D0) CALL XERROR ( 'DBETAI  X IS NOT IN T
     1HE RANGE (0,1)', 35, 1, 2)
      IF (PIN.LE.0.D0 .OR. QIN.LE.0.D0) CALL XERROR (  'DBETAI  P AND/OR
     1 Q IS LE ZERO', 29, 2, 2)
C
      Y = X
      P = PIN
      Q = QIN
      IF (Q.LE.P .AND. X.LT.0.8D0) GO TO 20
      IF (X.LT.0.2D0) GO TO 20
      Y = 1.0D0 - Y
      P = QIN
      Q = PIN
C
 20   IF ((P+Q)*Y/(P+1.D0).LT.EPS) GO TO 80
C
C EVALUATE THE INFINITE SUM FIRST.  TERM WILL EQUAL
C Y**P/BETA(PS,P) * (1.-PS)-SUB-I * Y**I / FAC(I) .
C
      PS = Q - DINT(Q)
      IF (PS.EQ.0.D0) PS = 1.0D0
      XB = P*DLOG(Y) - DLBETA(PS,P) - DLOG(P)
      DBETAI = 0.0D0
      IF (XB.LT.ALNSML) GO TO 40
C
      DBETAI = DEXP (XB)
      TERM = DBETAI*P
      IF (PS.EQ.1.0D0) GO TO 40
      N = MAX (SNGL(ALNEPS/DLOG(Y)), 4.0)
      DO 30 I=1,N
        XI = I
        TERM = TERM * (XI-PS)*Y/XI
        DBETAI = DBETAI + TERM/(P+XI)
 30   CONTINUE
C
C NOW EVALUATE THE FINITE SUM, MAYBE.
C
 40   IF (Q.LE.1.0D0) GO TO 70
C
      XB = P*DLOG(Y) + Q*DLOG(1.0D0-Y) - DLBETA(P,Q) - DLOG(Q)
      IB = MAX (SNGL(XB/ALNSML), 0.0)
      TERM = DEXP (XB - DBLE(FLOAT(IB))*ALNSML )
      C = 1.0D0/(1.D0-Y)
      P1 = Q*C/(P+Q-1.D0)
C
      FINSUM = 0.0D0
      N = Q
      IF (Q.EQ.DBLE(FLOAT(N))) N = N - 1
      DO 50 I=1,N
        IF (P1.LE.1.0D0 .AND. TERM/EPS.LE.FINSUM) GO TO 60
        XI = I
        TERM = (Q-XI+1.0D0)*C*TERM/(P+Q-XI)
C
        IF (TERM.GT.1.0D0) IB = IB - 1
        IF (TERM.GT.1.0D0) TERM = TERM*SML
C
        IF (IB.EQ.0) FINSUM = FINSUM + TERM
 50   CONTINUE
C
 60   DBETAI = DBETAI + FINSUM
 70   IF (Y.NE.X .OR. P.NE.PIN) DBETAI = 1.0D0 - DBETAI
      DBETAI = DMAX1 (DMIN1 (DBETAI, 1.0D0), 0.0D0)
      RETURN
C
 80   DBETAI = 0.0D0
      XB = P*DLOG(DMAX1(Y,SML)) - DLOG(P) - DLBETA(P,Q)
      IF (XB.GT.ALNSML .AND. Y.NE.0.0D0) DBETAI = DEXP(XB)
      IF (Y.NE.X .OR. P.NE.PIN) DBETAI = 1.0D0 - DBETAI
C
      RETURN
      END
