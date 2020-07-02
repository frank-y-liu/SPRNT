      SUBROUTINE EXINT(X,N,KODE,M,TOL,EN,IERR)
C***BEGIN PROLOGUE  EXINT
C***DATE WRITTEN   800501   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  C5
C***KEYWORDS  EXPONENTIAL INTEGRAL,SPECIAL FUNCTION
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  EXINT computes M member sequences of exponential integrals
C            E(N+K,X), K=0,1,...,M-1 for N .GE. 1 and X .GE. 0.
C***DESCRIPTION
C
C     Written by D. E. Amos, Sandia Laboratories, Albuquerque,
C     New Mexico, 87185
C
C     Reference
C         Computation of Exponential Integrals by D. E. Amos, ACM
C         Trans. Math Software, 6, 1980
C
C     Abstract
C         EXINT computes M member sequences of exponential integrals
C         E(N+K,X), K=0,1,...,M-1 for N .GE. 1 and X .GE. 0.  The power
C         series is implemented for X .LE. XCUT and the confluent
C         hypergeometric representation
C
C                     E(A,X) = EXP(-X)*(X**(A-1))*U(A,A,X)
C
C         is computed for X > XCUT.  Since sequences are computed in a
C         stable fashion by recurring away from X, A is selected as the
C         integer closest to X within the constraint N .LE. A .LE. N+M-1
C         For the U computation,  A is further modified to be the
C         nearest even integer.  Indices are carried forward or
C         backward by the two term recursion relation
C
C                     K*E(K+1,X) + X*E(K,X) = EXP(-X)
C
C         once E(A,X) is computed.  The U function is computed by means
C         of the backward recursive Miller algorithm applied to the
C         three term contiguous relation for U(A+K,A,X), K=0,1,...
C         This produces accurate ratios and determines U(A+K,A,X), and
C         hence E(A,X), to within a multiplicative constant C.
C         Another contiguous relation applied to C*U(A,A,X) and
C         C*U(A+1,A,X) gets C*U(A+1,A+1,X), a quantity proportional to
C         E(A+1,X).  The normalizing constant C is obtained from the
C         two term recursion relation above with K=A.
C
C         EXINT calls I1MACH, R1MACH, PSIXN, XERROR
C
C     Description of Arguments
C
C         Input
C           X       X .GT. 0.0 for N=1 and  X .GE. 0.0 for N .GE. 2
C           N       order of the first member of the sequence, N .GE. 1
C           KODE    a selection parameter for scaled values
C                   KODE=1   returns        E(N+K,X), K=0,1,...,M-1.
C                       =2   returns EXP(X)*E(N+K,X), K=0,1,...,M-1.
C           M       number of exponential integrals in the sequence,
C                   M .GE. 1
C           TOL     relative accuracy wanted, ETOL .LE. TOL .LE. 0.1
C                   ETOL = single precision unit roundoff = R1MACH(4)
C
C         Output
C           EN      a vector of dimension at least M containing values
C                   EN(K) = E(N+K-1,X) or EXP(X)*E(N+K-1,X), K=1,M
C                   depending on KODE
C           IERR    underflow indicator
C                   IERR=0   a normal return
C                       =1   X exceeds XLIM and an underflow occurs.
C                            EN(K)=0.0E0 , K=1,M returned on KODE=1
C
C     Error Conditions
C         An improper input parameter is a fatal error.
C         Underflow is a non fatal error.  Zero answers are returned.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  I1MACH,PSIXN,R1MACH,XERROR
C***END PROLOGUE  EXINT
C
      REAL             A,AA,AAMS,AH,AK,AT,B,BK,BT,CC,CNORM,CT,EM,EMX,EN,
     1                 ETOL,FNM,FX,PT,P1,P2,S,TOL,TX,X,XCUT,XLIM,XTOL,Y,
     2                 YT,Y1,Y2
      REAL             R1MACH,PSIXN
      INTEGER I,IC,ICASE,ICT,IERR,IK,IND,IX,I1M,JSET,K,KK,KN,KODE,KS,M,
     1        ML,MU,N,ND,NM
      INTEGER I1MACH
      DIMENSION EN(*), A(99), B(99), Y(2)
C***FIRST EXECUTABLE STATEMENT  EXINT
      ETOL=R1MACH(4)
      I1M=-I1MACH(12)
      XLIM=2.302E0*(FLOAT(I1M)*R1MACH(5)-3.0E0)
      IF (N.LT.1) GO TO 280
      IF (KODE.LT.1 .OR. KODE.GT.2) GO TO 290
      IF (M.LT.1) GO TO 300
      IF(TOL.LT.ETOL .OR. TOL.GT.0.1E0) GO TO 310
C
      IERR = 0
      XCUT=2.0E0
      IF(ETOL.LT.1.0E-8) XCUT=1.0E0
      IF (X.GT.XCUT) GO TO 100
      IF(X.LT.0.0E0) GO TO 320
      IF(X.EQ.0.0E0 .AND. N.EQ.1) GO TO 330
      IF(X.EQ.0.0E0 .AND. N.GT.1) GO TO 80
C
C     SERIES FOR E(N,X) FOR X.LE.XCUT
C
      TX=X+0.5E0
      IX=INT(TX)
C     ICASE=1 MEANS INTEGER CLOSEST TO X IS 2 AND N=1
C     ICASE=2 MEANS INTEGER CLOSEST TO X IS 0,1, OR 2 AND N.GE.2
      ICASE = 2
      IF (IX.GT.N) ICASE = 1
      NM = N - ICASE + 1
      ND = NM + 1
      IND = 3 - ICASE
      MU = M - IND
      ML = 1
      KS = ND
      FNM=FLOAT(NM)
      S = 0.0E0
      XTOL=3.0E0*TOL
      IF (ND.EQ.1) GO TO 10
      XTOL=0.3333E0*TOL
      S=1.0E0/FNM
   10 CONTINUE
      AA=1.0E0
      AK=1.0E0
      IC=35
      IF(X.LT.ETOL) IC=1
      DO 50 I=1,IC
        AA = -AA*X/AK
        IF (I.EQ.NM) GO TO 30
        S = S - AA/(AK-FNM)
        IF (ABS(AA).LE.XTOL*ABS(S)) GO TO 20
        AK=AK+1.0E0
        GO TO 50
   20   CONTINUE
        IF (I.LT.2) GO TO 40
        IF (ND-2.GT.I .OR. I.GT.ND-1) GO TO 60
        AK = AK + 1.0E0
        GO TO 50
   30   S = S + AA*(-ALOG(X)+PSIXN(ND))
        XTOL=3.0E0*TOL
   40   AK = AK + 1.0E0
   50 CONTINUE
      IF(IC.NE.1) GO TO 340
   60 IF (ND.EQ.1) S = S + (-ALOG(X)+PSIXN(1))
      IF (KODE.EQ.2) S = S*EXP(X)
      EN(1) = S
      EMX=1.0E0
      IF (M.EQ.1) GO TO 70
      EN(IND) = S
      AA=FLOAT(KS)
      IF (KODE.EQ.1) EMX = EXP(-X)
      GO TO (240, 260), ICASE
   70 IF (ICASE.EQ.2) RETURN
      IF (KODE.EQ.1) EMX = EXP(-X)
      EN(1) = (EMX-S)/X
      RETURN
   80 CONTINUE
      DO 90 I=1,M
        EN(I)=1.0E0/FLOAT(N+I-2)
   90 CONTINUE
      RETURN
C
C     BACKWARD RECURSIVE MILLER ALGORITHM FOR
C              E(N,X)=EXP(-X)*(X**(N-1))*U(N,N,X)
C     WITH RECURSION AWAY FROM N=INTEGER CLOSEST TO X.
C     U(A,B,X) IS THE SECOND CONFLUENT HYPERGEOMETRIC FUNCTION
C
  100 CONTINUE
      EMX=1.0E0
      IF (KODE.EQ.2) GO TO 130
      IF (X.LE.XLIM) GO TO 120
      IERR = 1
      DO 110 I=1,M
        EN(I)=0.0E0
  110 CONTINUE
      RETURN
  120 EMX = EXP(-X)
  130 CONTINUE
      IX=INT(X+0.5E0)
      KN = N + M - 1
      IF (KN.LE.IX) GO TO 140
      IF (N.LT.IX .AND. IX.LT.KN) GO TO 170
      IF (N.GE.IX) GO TO 160
      GO TO 360
  140 ICASE = 1
      KS = KN
      ML = M - 1
      MU = -1
      IND = M
      IF (KN.GT.1) GO TO 180
  150 KS = 2
      ICASE = 3
      GO TO 180
  160 ICASE = 2
      IND = 1
      KS = N
      MU = M - 1
      IF (N.GT.1) GO TO 180
      IF (KN.EQ.1) GO TO 150
      IX = 2
  170 ICASE = 1
      KS = IX
      ML = IX - N
      IND = ML + 1
      MU = KN - IX
  180 CONTINUE
      IK = KS/2
      AH=FLOAT(IK)
      JSET=1+KS-(IK+IK)
C     START COMPUTATION FOR
C              EN(IND) = C*U( A , A ,X)    JSET=1
C              EN(IND) = C*U(A+1,A+1,X)    JSET=2
C     FOR AN EVEN INTEGER A.
      IC = 0
      AA = AH + AH
      AAMS=AA-1.0E0
      AAMS=AAMS*AAMS
      TX=X+X
      FX=TX+TX
      AK=AH
      XTOL=TOL
      IF(TOL.LE.1.0E-3) XTOL=20.0E0*TOL
      CT=AAMS+FX*AH
      EM=(AH+1.0E0)/((X+AA)*XTOL*SQRT(CT))
      BK=AA
      CC=AH*AH
C     FORWARD RECURSION FOR P(IC),P(IC+1) AND INDEX IC FOR BACKWARD
C     RECURSION
      P1=0.0E0
      P2=1.0E0
  200 CONTINUE
      IF(IC.EQ.99) GO TO 350
      IC=IC+1
      AK=AK+1.0E0
      AT=BK/(BK+AK+CC+FLOAT(IC))
      BK=BK+AK+AK
      A(IC)=AT
      BT=(AK+AK+X)/(AK+1.0E0)
      B(IC)=BT
      PT=P2
      P2=BT*P2-AT*P1
      P1=PT
      CT=CT+FX
      EM=EM*AT*(1.0E0-TX/CT)
      IF(EM*(AK+1.0E0) .GT. P1*P1) GO TO 200
      ICT = IC
      KK=IC+1
      BT=TX/(CT+FX)
      Y2=(BK/(BK+CC+FLOAT(KK)))*(P1/P2)*(1.0E0-BT+0.375E0*BT*BT)
      Y1 = 1.0E0
C     BACKWARD RECURRENCE FOR
C              Y1=             C*U( A ,A,X)
C              Y2= C*(A/(1+A/2))*U(A+1,A,X)
      DO 220 K=1,ICT
        KK=KK-1
        YT = Y1
        Y1=(B(KK)*Y1-Y2)/A(KK)
        Y2 = YT
  220 CONTINUE
C     THE CONTIGUOUS RELATION
C              X*U(B,C+1,X)=(C-B)*U(B,C,X)+U(B-1,C,X)
C     WITH  B=A+1 , C=A IS USED FOR
C              Y(2) = C * U(A+1,A+1,X)
C     X IS INCORPORATED INTO THE NORMALIZING RELATION FOR CNORM.
      Y(1) = Y1
      Y(2)=Y1-Y2*(AH+1.0E0)/AA
      CNORM = EMX/(AA*Y(2)+X*Y(1))
      IF (ICASE.EQ.3) GO TO 230
      EN(IND) = CNORM*Y(JSET)
      IF (M.EQ.1) RETURN
      AA=FLOAT(KS)
      GO TO (240, 260), ICASE
C
C     RECURSION SECTION  N*E(N+1,X) + X*E(N,X)=EMX
C
  230 EN(1) = (EMX-CNORM*Y(1))/X
      RETURN
  240 K = IND - 1
      DO 250 I=1,ML
        AA=AA-1.0E0
        EN(K) = (EMX-AA*EN(K+1))/X
        K = K - 1
  250 CONTINUE
      IF (MU.LE.0) RETURN
      AA=FLOAT(KS)
  260 K = IND
      DO 270 I=1,MU
        EN(K+1) = (EMX-X*EN(K))/AA
        AA=AA+1.0E0
        K = K + 1
  270 CONTINUE
      RETURN
C
C
  280 CALL XERROR( 'IN EXINT,  N NOT GREATER THAN 0',31,2,1)
      RETURN
  290 CALL XERROR( 'IN EXINT,  KODE NOT 1 OR 2',26,2,1)
      RETURN
  300 CALL XERROR( 'IN EXINT,  M NOT GREATER THAN 0',31,2,1)
      RETURN
  310 CALL XERROR( 'IN EXINT,  TOL NOT WITHIN LIMITS',32,2,1)
      RETURN
  320 CALL XERROR( 'IN EXINT,  X IS NOT ZERO OR POSITIVE',36,2,1)
      RETURN
  330 CALL XERROR( 'IN EXINT,  THE EXPONENTIAL INTEGRAL IS NOT DEFINED F
     1OR X=0 AND N=1',66,2,1)
      RETURN
  340 CALL XERROR( 'IN EXINT,  RELATIVE ERROR TEST FOR SERIES TERMINATIO
     1N NOT MET IN 36 TERMS',73,3,1)
      RETURN
  350 CALL XERROR( 'IN EXINT,  TERMINATION TEST FOR MILLER ALGORITHM NOT
     1 MET IN 99 STEPS',68,3,1)
      RETURN
  360 CALL XERROR( 'IN EXINT,  AN ERROR IN PLACING INT(X+0.5) WITH RESPE
     1CT TO N AND N+M-1 OCCURRED FOR X.GT.XCUT',92,8,1)
      RETURN
      END
