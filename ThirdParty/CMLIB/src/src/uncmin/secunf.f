      SUBROUTINE SECUNF(NR,N,X,G,A,UDIAG,XPLS,GPLS,EPSM,ITNCNT,
     +     RNF,IAGFLG,NOUPDT,S,Y,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C UPDATE HESSIAN BY THE BFGS UNFACTORED METHOD
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> OLD ITERATE, X[K-1]
C G(N)         --> GRADIENT OR APPROXIMATE AT OLD ITERATE
C A(N,N)      <--> ON ENTRY: APPROXIMATE HESSIAN AT OLD ITERATE
C                    IN UPPER TRIANGULAR PART (AND UDIAG)
C                  ON EXIT:  UPDATED APPROX HESSIAN AT NEW ITERATE
C                    IN LOWER TRIANGULAR PART AND DIAGONAL
C                  [LOWER TRIANGULAR PART OF SYMMETRIC MATRIX]
C UDIAG        --> ON ENTRY: DIAGONAL OF HESSIAN
C XPLS(N)      --> NEW ITERATE, X[K]
C GPLS(N)      --> GRADIENT OR APPROXIMATE AT NEW ITERATE
C EPSM         --> MACHINE EPSILON
C ITNCNT       --> ITERATION COUNT
C RNF          --> RELATIVE NOISE IN OPTIMIZATION FUNCTION FCN
C IAGFLG       --> =1 IF ANALYTIC GRADIENT SUPPLIED, =0 OTHERWISE
C NOUPDT      <--> BOOLEAN: NO UPDATE YET
C                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C S(N)         --> WORKSPACE
C Y(N)         --> WORKSPACE
C T(N)         --> WORKSPACE
C
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
      DIMENSION X(N),G(N),XPLS(N),GPLS(N)
      DIMENSION A(NR,*)
      DIMENSION UDIAG(N)
      DIMENSION S(N),Y(N),T(N)
      LOGICAL NOUPDT,SKPUPD
C
C COPY HESSIAN IN UPPER TRIANGULAR PART AND UDIAG TO
C LOWER TRIANGULAR PART AND DIAGONAL
C
      DO 5 J=1,N
        A(J,J)=UDIAG(J)
        IF(J.EQ.N) GO TO 5
        JP1=J+1
        DO 4 I=JP1,N
          A(I,J)=A(J,I)
    4   CONTINUE
    5 CONTINUE
C
      IF(ITNCNT.EQ.1) NOUPDT=.TRUE.
      DO 10 I=1,N
        S(I)=XPLS(I)-X(I)
        Y(I)=GPLS(I)-G(I)
   10 CONTINUE
      DEN1=DDOT(N,S,1,Y,1)
      SNORM2=DNRM2(N,S,1)
      YNRM2=DNRM2(N,Y,1)
      IF(DEN1.LT.SQRT(EPSM)*SNORM2*YNRM2) GO TO 100
C     IF(DEN1.GE.SQRT(EPSM)*SNORM2*YNRM2)
C     THEN
        CALL MVMLTS(NR,N,A,S,T)
        DEN2=DDOT(N,S,1,T,1)
        IF(.NOT. NOUPDT) GO TO 50
C       IF(NOUPDT)
C       THEN
C
C         H <-- [(S+)Y/(S+)HS]H
C
          GAM=DEN1/DEN2
          DEN2=GAM*DEN2
          DO 30 J=1,N
            T(J)=GAM*T(J)
            DO 20 I=J,N
              A(I,J)=GAM*A(I,J)
   20       CONTINUE
   30     CONTINUE
          NOUPDT=.FALSE.
C       ENDIF
   50   SKPUPD=.TRUE.
C
C       CHECK UPDATE CONDITION ON ROW I
C
        DO 60 I=1,N
          TOL=RNF*MAX(ABS(G(I)),ABS(GPLS(I)))
          IF(IAGFLG.EQ.0) TOL=TOL/SQRT(RNF)
          IF(ABS(Y(I)-T(I)).LT.TOL) GO TO 60
C         IF(ABS(Y(I)-T(I)).GE.TOL)
C         THEN
            SKPUPD=.FALSE.
            GO TO 70
C         ENDIF
   60   CONTINUE
   70   IF(SKPUPD) GO TO 100
C       IF(.NOT.SKPUPD)
C       THEN
C
C         BFGS UPDATE
C
          DO 90 J=1,N
            DO 80 I=J,N
              A(I,J)=A(I,J)+Y(I)*Y(J)/DEN1-T(I)*T(J)/DEN2
   80       CONTINUE
   90     CONTINUE
C       ENDIF
C     ENDIF
  100 RETURN
      END
