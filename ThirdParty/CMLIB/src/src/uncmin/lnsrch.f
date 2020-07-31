      SUBROUTINE LNSRCH(N,X,F,G,P,XPLS,FPLS,FCN,MXTAKE,
     +   IRETCD,STEPMX,STEPTL,SX,IPR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C PURPOSE
C -------
C FIND A NEXT NEWTON ITERATE BY LINE SEARCH.
C
C PARAMETERS
C ----------
C N            --> DIMENSION OF PROBLEM
C X(N)         --> OLD ITERATE:   X[K-1]
C F            --> FUNCTION VALUE AT OLD ITERATE, F(X)
C G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE
C P(N)         --> NON-ZERO NEWTON STEP
C XPLS(N)     <--  NEW ITERATE X[K]
C FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS)
C FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION
C IRETCD      <--  RETURN CODE
C MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED
C STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE
C STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C
C INTERNAL VARIABLES
C ------------------
C SLN              NEWTON LENGTH
C RLN              RELATIVE LENGTH OF NEWTON STEP
C
      INTEGER N,IRETCD
      DIMENSION SX(N)
      DIMENSION X(N),G(N),P(N)
      DIMENSION XPLS(N)
      LOGICAL MXTAKE
C
      IPR=IPR
      MXTAKE=.FALSE.
      IRETCD=2
C$    WRITE(IPR,954)
C$    WRITE(IPR,955) (P(I),I=1,N)
      TMP=0.0
      DO 5 I=1,N
        TMP=TMP+SX(I)*SX(I)*P(I)*P(I)
    5 CONTINUE
      SLN=SQRT(TMP)
      IF(SLN.LE.STEPMX) GO TO 10
C
C NEWTON STEP LONGER THAN MAXIMUM ALLOWED
        SCL=STEPMX/SLN
        CALL SCLMUL(N,SCL,P,P)
        SLN=STEPMX
C$      WRITE(IPR,954)
C$      WRITE(IPR,955) (P(I),I=1,N)
   10 CONTINUE
      SLP=DDOT(N,G,1,P,1)
      RLN=0.
      DO 15 I=1,N
        RLN=MAX(RLN,ABS(P(I))/MAX(ABS(X(I)),1./SX(I)))
   15 CONTINUE
      RMNLMB=STEPTL/RLN
      ALMBDA=1.0
C$    WRITE(IPR,952) SLN,SLP,RMNLMB,STEPMX,STEPTL
C
C LOOP
C CHECK IF NEW ITERATE SATISFACTORY.  GENERATE NEW LAMBDA IF NECESSARY.
C
  100 CONTINUE
      IF(IRETCD.LT.2) RETURN
      DO 105 I=1,N
        XPLS(I)=X(I) + ALMBDA*P(I)
  105 CONTINUE
      CALL FCN(N,XPLS,FPLS)
C$    WRITE(IPR,950) ALMBDA
C$    WRITE(IPR,951)
C$    WRITE(IPR,955) (XPLS(I),I=1,N)
C$    WRITE(IPR,953) FPLS
      IF(FPLS.GT. F+SLP*1.E-4*ALMBDA) GO TO 130
C     IF(FPLS.LE. F+SLP*1.E-4*ALMBDA)
C     THEN
C
C SOLUTION FOUND
C
        IRETCD=0
        IF(ALMBDA.EQ.1.0 .AND. SLN.GT. .99*STEPMX) MXTAKE=.TRUE.
        GO TO 100
C
C SOLUTION NOT (YET) FOUND
C
C     ELSE
  130   IF(ALMBDA .GE. RMNLMB) GO TO 140
C       IF(ALMBDA .LT. RMNLMB)
C       THEN
C
C NO SATISFACTORY XPLS FOUND SUFFICIENTLY DISTINCT FROM X
C
          IRETCD=1
          GO TO 100
C       ELSE
C
C CALCULATE NEW LAMBDA
C
  140     IF(ALMBDA.NE.1.0) GO TO 150
C         IF(ALMBDA.EQ.1.0)
C         THEN
C
C FIRST BACKTRACK: QUADRATIC FIT
C
            TLMBDA=-SLP/(2.*(FPLS-F-SLP))
            GO TO 170
C         ELSE
C
C ALL SUBSEQUENT BACKTRACKS: CUBIC FIT
C
  150       T1=FPLS-F-ALMBDA*SLP
            T2=PFPLS-F-PLMBDA*SLP
            T3=1.0/(ALMBDA-PLMBDA)
            A=T3*(T1/(ALMBDA*ALMBDA) - T2/(PLMBDA*PLMBDA))
            B=T3*(T2*ALMBDA/(PLMBDA*PLMBDA)
     +           - T1*PLMBDA/(ALMBDA*ALMBDA) )
            DISC=B*B-3.0*A*SLP
            IF(DISC.LE. B*B) GO TO 160
C           IF(DISC.GT. B*B)
C           THEN
C
C ONLY ONE POSITIVE CRITICAL POINT, MUST BE MINIMUM
C
              TLMBDA=(-B+SIGN(1.0D0,A)*SQRT(DISC))/(3.0*A)
              GO TO 165
C           ELSE
C
C BOTH CRITICAL POINTS POSITIVE, FIRST IS MINIMUM
C
  160         TLMBDA=(-B-SIGN(1.0D0,A)*SQRT(DISC))/(3.0*A)
C           ENDIF
  165       IF(TLMBDA.GT. .5*ALMBDA) TLMBDA=.5*ALMBDA
C         ENDIF
  170     PLMBDA=ALMBDA
          PFPLS=FPLS
          IF(TLMBDA.GE. ALMBDA*.1) GO TO 180
C         IF(TLMBDA.LT.ALMBDA/10.)
C         THEN
            ALMBDA=ALMBDA*.1
            GO TO 190
C         ELSE
  180       ALMBDA=TLMBDA
C         ENDIF
C       ENDIF
C     ENDIF
  190 GO TO 100
  950 FORMAT(18H LNSRCH    ALMBDA=,E20.13)
  951 FORMAT(29H LNSRCH    NEW ITERATE (XPLS))
  952 FORMAT(18H LNSRCH    SLN   =,E20.13/
     +       18H LNSRCH    SLP   =,E20.13/
     +       18H LNSRCH    RMNLMB=,E20.13/
     +       18H LNSRCH    STEPMX=,E20.13/
     +       18H LNSRCH    STEPTL=,E20.13)
  953 FORMAT(19H LNSRCH    F(XPLS)=,E20.13)
  954 FORMAT(26H0LNSRCH    NEWTON STEP (P))
  955 FORMAT(14H LNSRCH       ,5(E20.13,3X))
      END