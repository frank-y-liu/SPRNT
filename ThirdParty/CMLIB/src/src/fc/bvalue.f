      FUNCTION BVALUE(T,A,N,K,X,IDERIV)
C***BEGIN PROLOGUE  BVALUE
C***REFER TO  FC
C
C Calculates value at *X* of *IDERIV*-th deriv. of spline from B-repr.
C***ROUTINES CALLED  INTERV
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***END PROLOGUE  BVALUE
      DIMENSION T(*),A(*)
      DIMENSION AJ(20),DP(20),DM(20)
C***FIRST EXECUTABLE STATEMENT  BVALUE
      BVALUE = 0.
      KMIDER = K - IDERIV
      IF (KMIDER .LE. 0)               GO TO 99
C
C  *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
C      (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
      KM1 = K-1
      CALL INTERV (T(K), N+1-KM1, X, I, MFLAG )
      I = I + KM1
      IF (MFLAG)                       99,20,9
    9 IF (X .GT. T(I))                 GO TO 99
   10 IF (I .EQ. K)                    GO TO 99
      I = I - 1
      IF (X .EQ. T(I))                 GO TO 10
C
C  *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES .
   20 IMK = I-K
      DO 21 J=1,K
   21    AJ(J) = A(IMK+J)
      IF (IDERIV .EQ. 0)               GO TO 30
   22 DO 23 J=1,IDERIV
         KMJ = K-J
         FKMJ = FLOAT(KMJ)
         DO 23 JJ=1,KMJ
            IHI = I + JJ
   23       AJ(JJ) = (AJ(JJ+1) - AJ(JJ))/(T(IHI) - T(IHI-KMJ))*FKMJ
C
C  *** COMPUTE VALUE AT *X* IN (T(I),T(I+1)) OF IDERIV-TH DERIVATIVE,
C      GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
   30 IF (IDERIV .EQ. KM1)             GO TO 39
      IP1 = I+1
      DO 32 J=1,KMIDER
         DP(J) = T(I+J) - X
   32    DM(J) = X - T(IP1-J)
      IDERP1 = IDERIV+1
      DO 33 J=IDERP1,KM1
         KMJ = K-J
         ILO = KMJ
         DO 33 JJ=1,KMJ
            AJ(JJ) = (AJ(JJ+1)*DM(ILO) + AJ(JJ)*DP(JJ))/(DM(ILO)+DP(JJ))
   33       ILO = ILO - 1
   39 BVALUE = AJ(1)
C
   99                                  RETURN
      END
