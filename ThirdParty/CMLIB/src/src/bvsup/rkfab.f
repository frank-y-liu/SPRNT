      SUBROUTINE RKFAB(NCOMP,XPTS,NXPTS,NFC,IFLAG,Z,MXNON,P,NTP,IP,YHP,
     1   NIV,U,V,W,S,STOWA,G,WORK,IWORK,NFCC)
C***BEGIN PROLOGUE  RKFAB
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***REFER TO  BVSUP
C
C **********************************************************************
C
C     Subroutine RKFAB integrates the initial value equations using
C     the variable-step RUNGE-KUTTA-FEHLBERG integration scheme or
C     the variable-order ADAMS method and orthonormalization
C     determined by a linear dependence test.
C
C **********************************************************************
C***ROUTINES CALLED  BVDER,DEABM,DERKF,REORT,STOR1
C***COMMON BLOCKS    ML15TO,ML17BW,ML18JR,ML8SZ
C***END PROLOGUE  RKFAB
C
      DIMENSION P(NTP,MXNON),IP(NFCC,MXNON),U(NCOMP,NFC,NXPTS),
     1          V(NCOMP,NXPTS),W(NFCC,MXNON),Z(MXNON),YHP(NCOMP,*),
     2          XPTS(NXPTS),S(*),STOWA(*),WORK(*),IWORK(*),
     3          G(NCOMP)
C
C **********************************************************************
C
      COMMON /ML8SZ/ NFCD,NCOMPD,INHOMO,IGOFX,XSAV,C,IVP
      COMMON /ML15TO/XBEG,X,INFO(15),KOP,XOP,ISTKOP,MNSWOT,
     1               NSWOT,KNSWOT,LOTJP,PWCND,TND,XOT,PX,XEND
      COMMON /ML18JR/NXPTSD,NIC,RE,AE,NOPG,MXNOND,NDISK,NTAPE,NEQ,
     1              INDPVT,INTEG,TOL,NPS,NTPD,NEQIVP,NUMORT,NFCCD,ICOCO
      COMMON /ML17BW/ KKKZPW,NEEDW,NEEDIW,K1,K2,K3,K4,K5,K6,K7,K8,K9,
     1                K10,K11,L1,L2,KKKINT,LLLINT
C
      EXTERNAL BVDER
C
C **********************************************************************
C  INITIALIZATION OF COUNTERS AND VARIABLES.
C
C***FIRST EXECUTABLE STATEMENT  RKFAB
      KOD = 1
      NON = 1
      X = XBEG
      JON = 1
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 1
      INFO(4) = 1
      WORK(1) = XEND
      IF (NOPG .EQ. 0)  GO TO 1
      INFO(3) = 0
      IF (X .EQ. Z(1))  JON = 2
    1 NFCP1 = NFC + 1
C
C **********************************************************************
C *****BEGINNING OF INTEGRATION LOOP AT OUTPUT POINTS.******************
C **********************************************************************
C
      DO 110 KOPP = 2,NXPTS
      KOP=KOPP
C
    5 XOP = XPTS(KOP)
      IF (NDISK .EQ. 0)  KOD = KOP
C
C     STEP BY STEP INTEGRATION LOOP BETWEEN OUTPUT POINTS.
C
   10 XXOP = XOP
      IF (NOPG .EQ. 0)   GO TO 15
      IF (XEND.GT.XBEG.AND.XOP.GT.Z(JON)) XXOP=Z(JON)
      IF (XEND.LT.XBEG.AND.XOP.LT.Z(JON)) XXOP=Z(JON)
C
C **********************************************************************
   15 GO TO (20,25),INTEG
C     DERKF INTEGRATOR
C
   20 CALL DERKF(BVDER,NEQ,X,YHP,XXOP,INFO,RE,AE,IDID,WORK,KKKINT,
     1           IWORK,LLLINT,G,IPAR)
      GO TO 28
C     DEABM INTEGRATOR
C
   25 CALL DEABM(BVDER,NEQ,X,YHP,XXOP,INFO,RE,AE,IDID,WORK,KKKINT,
     1           IWORK,LLLINT,G,IPAR)
   28 IF(IDID .GE. 1) GO TO 30
      INFO(1) = 1
      IF(IDID .EQ. -1) GO TO 15
      IFLAG = 20 - IDID
      RETURN
C
C **********************************************************************
C     GRAM-SCHMIDT ORTHOGONALIZATION TEST FOR ORTHONORMALIZATION
C     (TEMPORARILY USING U AND V IN THE TEST)
C
   30 IF (NOPG .EQ. 0)  GO TO 35
      IF (XXOP .NE. Z(JON))  GO TO 100
      JFLAG=2
      GO TO 40
   35 JFLAG=1
      IF (INHOMO .EQ. 3  .AND.  X .EQ. XEND) JFLAG=3
C
   40 IF (NDISK .EQ. 0) NON=NUMORT+1
      CALL REORT(NCOMP,U(1,1,KOD),V(1,KOD),YHP,NIV,
     1           W(1,NON),S,P(1,NON),IP(1,NON),STOWA,JFLAG)
C
      IF (JFLAG .NE. 30) GO TO 45
      IFLAG=30
      RETURN
C
   45 IF (JFLAG .EQ. 10) GO TO 5
C
      IF (JFLAG .NE. 0)  GO TO 100
C
C **********************************************************************
C     STORE ORTHONORMALIZED VECTORS INTO SOLUTION VECTORS.
C
      IF (NUMORT .LT. MXNON)  GO TO 65
      IF (X .EQ. XEND) GO TO 65
      IFLAG = 13
      RETURN
C
   65 NUMORT = NUMORT + 1
      CALL STOR1(YHP,U(1,1,KOD),YHP(1,NFCP1),V(1,KOD),1,
     1           NDISK,NTAPE)
C
C **********************************************************************
C     STORE ORTHONORMALIZATION INFORMATION, INITIALIZE
C     INTEGRATION FLAG, AND CONTINUE INTEGRATION TO THE NEXT
C     ORTHONORMALIZATION POINT OR OUTPUT POINT.
C
      Z(NUMORT) = X
      IF (INHOMO .EQ. 1  .AND.  NPS .EQ. 0)  C = S(NFCP1) * C
      IF (NDISK .EQ. 0)  GO TO 90
      IF (INHOMO .EQ. 1)  WRITE (NTAPE) (W(J,1), J = 1,NFCC)
      WRITE(NTAPE) (IP(J,1), J = 1,NFCC),(P(J,1), J = 1,NTP)
   90 INFO(1) = 0
      JON = JON + 1
      IF (NOPG .EQ. 1  .AND.  X .NE. XOP)  GO TO 10
C
C **********************************************************************
C     CONTINUE INTEGRATION IF WE ARE NOT AT AN OUTPUT POINT.
C
  100 IF (IDID .EQ. 1)  GO TO 15
C
C     STORAGE OF HOMOGENEOUS SOLUTIONS IN U AND THE PARTICULAR
C     SOLUTION IN V AT THE OUTPUT POINTS.
C
      CALL STOR1(U(1,1,KOD),YHP,V(1,KOD),YHP(1,NFCP1),0,NDISK,NTAPE)
  110 CONTINUE
C **********************************************************************
C **********************************************************************
C
      IFLAG = 0
      RETURN
      END
