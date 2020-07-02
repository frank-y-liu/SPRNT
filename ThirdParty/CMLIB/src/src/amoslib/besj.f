      SUBROUTINE BESJ(X,ALPHA,N,Y,NZ)
C***BEGIN PROLOGUE  BESJ
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  C10A3
C***KEYWORDS  BESSEL FUNCTION,J BESSEL FUNCTION,SPECIAL FUNCTION
C***AUTHOR  AMOS, D. E., (SNLA)
C           DANIEL, S. L., (SNLA)
C           WESTON, M. K., (SNLA)
C***PURPOSE  BESJ computes an N member sequence of J Bessel functions
C            J/SUB(ALPHA+K-1)/(X), K=1,..,N for non-negative ALPHA and X
C***DESCRIPTION
C
C     Written by D. E. Amos, S. L. Daniel and M. K. Weston, January 1975
C
C     References
C         SAND-75-0147
C
C         CDC 6600 subroutines IBESS and JBESS for Bessel functions
C         I(NU,X) and J(NU,X), X .GE. 0, NU .GE. 0  by D. E. Amos, S. L.
C         Daniel, M. K. Weston,  ACM Transactions on Mathematical
C         Software, Vol. 3, pp. 76-92 (1977).
C
C         Tables of Bessel Functions of Moderate or Large Orders,
C         NPL Mathematical Tables, Vol. 6, by F. W. J. Olver, Her
C         Majesty's Stationery Office, London, 1962.
C
C     Abstract
C         BESJ computes an N member sequence of J bessel functions
C         J/sub(ALPHA+K-1)/(X), K=1,...,N for non-negative ALPHA and X.
C         A combination of the power series, the asymptotic expansion
C         for X to infinity and the uniform asymptotic expansion for
C         NU to infinity are applied over subdivisions of the (NU,X)
C         plane.  For values of (NU,X) not covered by one of these
C         formulae, the order is incremented or decremented by integer
C         values into a region where one of the formulae apply. Backward
C         recursion is applied to reduce orders by integer values except
C         where the entire sequence lies in the oscillatory region.  In
C         this case forward recursion is stable and values from the
C         asymptotic expansion for X to infinity start the recursion
C         when it is efficient to do so.  Leading terms of the series
C         and uniform expansion are tested for underflow.  If a sequence
C         is requested and the last member would underflow, the result
C         is set to zero and the next lower order tried, etc., until a
C         member comes on scale or all members are set to zero.
C         Overflow cannot occur.
C
C         BESJ calls ASYJY, JAIRY, ALNGAM, R1MACH, I1MACH, XERROR
C
C     Description of Arguments
C
C         Input
C           X      - X .GE. 0.0E0
C           ALPHA  - order of first member of the sequence,
C                    ALPHA .GE. 0.0E0
C           N      - number of members in the sequence, N .GE. 1
C
C         Output
C           Y      - a vector whose first  N components contain
C                    values for J/sub(ALPHA+K-1)/(X), K=1,...,N
C           NZ     - number of components of Y set to zero due to
C                    underflow,
C                    NZ=0   , normal return, computation completed
C                    NZ .NE. 0, last NZ components of Y set to zero,
C                             Y(K)=0.0E0, K=N-NZ+1,...,N.
C
C     Error Conditions
C         Improper input arguments - a fatal error
C         Underflow  - a non-fatal error (NZ .NE. 0)
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALNGAM,ASYJY,I1MACH,JAIRY,R1MACH,XERROR
C***END PROLOGUE  BESJ
C
      EXTERNAL JAIRY
      INTEGER I,IALP,IDALP,IFLW,IN,INLIM,IS,I1,I2,K,KK,KM,KT,N,NN,
     1        NS,NZ
      INTEGER I1MACH
      REAL       AK,AKM,ALPHA,ANS,AP,ARG,COEF,DALPHA,DFN,DTM,EARG,
     1           ELIM1,ETX,FIDAL,FLGJY,FN,FNF,FNI,FNP1,FNU,FNULIM,
     2           GLN,PDF,PIDT,PP,RDEN,RELB,RTTP,RTWO,RTX,RZDEN,
     3           S,SA,SB,SXO2,S1,S2,T,TA,TAU,TB,TEMP,TFN,TM,TOL,
     4           TOLLN,TRX,TX,T1,T2,WK,X,XO2,XO2L,Y
      REAL R1MACH, ALNGAM
      DIMENSION Y(*), TEMP(3), FNULIM(2), PP(4), WK(7)
      DATA RTWO,PDF,RTTP,PIDT                    / 1.34839972492648E+00,
     1 7.85398163397448E-01, 7.97884560802865E-01, 1.57079632679490E+00/
      DATA  PP(1),  PP(2),  PP(3),  PP(4)        / 8.72909153935547E+00,
     1 2.65693932265030E-01, 1.24578576865586E-01, 7.70133747430388E-04/
      DATA INLIM           /      150            /
      DATA FNULIM(1), FNULIM(2) /      100.0E0,     60.0E0     /
C***FIRST EXECUTABLE STATEMENT  BESJ
      NZ = 0
      KT = 1
C     I1MACH(14) REPLACES I1MACH(11) IN A DOUBLE PRECISION CODE
C     I1MACH(15) REPLACES I1MACH(12) IN A DOUBLE PRECISION CODE
      TA = R1MACH(3)
      TOL = AMAX1(TA,1.0E-15)
      I1 = I1MACH(11) + 1
      I2 = I1MACH(12)
      TB = R1MACH(5)
      ELIM1 = 2.303E0*(FLOAT(-I2)*TB-3.0E0)
C     TOLLN = -LN(TOL)
      TOLLN = 2.303E0*TB*FLOAT(I1)
      TOLLN = AMIN1(TOLLN,34.5388E0)
      IF (N-1) 720, 10, 20
   10 KT = 2
   20 NN = N
      IF (X) 730, 30, 80
   30 IF (ALPHA) 710, 40, 50
   40 Y(1) = 1.0E0
      IF (N.EQ.1) RETURN
      I1 = 2
      GO TO 60
   50 I1 = 1
   60 DO 70 I=I1,N
        Y(I) = 0.0E0
   70 CONTINUE
      RETURN
   80 CONTINUE
      IF (ALPHA.LT.0.0E0) GO TO 710
C
      IALP = INT(ALPHA)
      FNI = FLOAT(IALP+N-1)
      FNF = ALPHA - FLOAT(IALP)
      DFN = FNI + FNF
      FNU = DFN
      XO2 = X*0.5E0
      SXO2 = XO2*XO2
C
C     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
C     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE
C     APPLIED.
C
      IF (SXO2.LE.(FNU+1.0E0)) GO TO 90
      TA = AMAX1(20.0E0,FNU)
      IF (X.GT.TA) GO TO 120
      IF (X.GT.12.0E0) GO TO 110
      XO2L = ALOG(XO2)
      NS = INT(SXO2-FNU) + 1
      GO TO 100
   90 FN = FNU
      FNP1 = FN + 1.0E0
      XO2L = ALOG(XO2)
      IS = KT
      IF (X.LE.0.50E0) GO TO 330
      NS = 0
  100 FNI = FNI + FLOAT(NS)
      DFN = FNI + FNF
      FN = DFN
      FNP1 = FN + 1.0E0
      IS = KT
      IF (N-1+NS.GT.0) IS = 3
      GO TO 330
  110 ANS = AMAX1(36.0E0-FNU,0.0E0)
      NS = INT(ANS)
      FNI = FNI + FLOAT(NS)
      DFN = FNI + FNF
      FN = DFN
      IS = KT
      IF (N-1+NS.GT.0) IS = 3
      GO TO 130
  120 CONTINUE
      RTX = SQRT(X)
      TAU = RTWO*RTX
      TA = TAU + FNULIM(KT)
      IF (FNU.LE.TA) GO TO 480
      FN = FNU
      IS = KT
C
C     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY
C
  130 CONTINUE
      I1 = IABS(3-IS)
      I1 = MAX0(I1,1)
      FLGJY = 1.0E0
      CALL ASYJY(JAIRY,X,FN,FLGJY,I1,TEMP(IS),WK,IFLW)
      IF(IFLW.NE.0) GO TO 380
      GO TO (320, 450, 620), IS
  310 TEMP(1) = TEMP(3)
      KT = 1
  320 IS = 2
      FNI = FNI - 1.0E0
      DFN = FNI + FNF
      FN = DFN
      IF(I1.EQ.2) GO TO 450
      GO TO 130
C
C     SERIES FOR (X/2)**2.LE.NU+1
C
  330 CONTINUE
      GLN = ALNGAM(FNP1)
      ARG = FN*XO2L - GLN
      IF (ARG.LT.(-ELIM1)) GO TO 400
      EARG = EXP(ARG)
  340 CONTINUE
      S = 1.0E0
      IF (X.LT.TOL) GO TO 360
      AK = 3.0E0
      T2 = 1.0E0
      T = 1.0E0
      S1 = FN
      DO 350 K=1,17
        S2 = T2 + S1
        T = -T*SXO2/S2
        S = S + T
        IF (ABS(T).LT.TOL) GO TO 360
        T2 = T2 + AK
        AK = AK + 2.0E0
        S1 = S1 + FN
  350 CONTINUE
  360 CONTINUE
      TEMP(IS) = S*EARG
      GO TO (370, 450, 610), IS
  370 EARG = EARG*FN/XO2
      FNI = FNI - 1.0E0
      DFN = FNI + FNF
      FN = DFN
      IS = 2
      GO TO 340
C
C     SET UNDERFLOW VALUE AND UPDATE PARAMETERS
C
  380 Y(NN) = 0.0E0
      NN = NN - 1
      FNI = FNI - 1.0E0
      DFN = FNI + FNF
      FN = DFN
      IF (NN-1) 440, 390, 130
  390 KT = 2
      IS = 2
      GO TO 130
  400 Y(NN) = 0.0E0
      NN = NN - 1
      FNP1 = FN
      FNI = FNI - 1.0E0
      DFN = FNI + FNF
      FN = DFN
      IF (NN-1) 440, 410, 420
  410 KT = 2
      IS = 2
  420 IF (SXO2.LE.FNP1) GO TO 430
      GO TO 130
  430 ARG = ARG - XO2L + ALOG(FNP1)
      IF (ARG.LT.(-ELIM1)) GO TO 400
      GO TO 330
  440 NZ = N - NN
      RETURN
C
C     BACKWARD RECURSION SECTION
C
  450 CONTINUE
      NZ = N - NN
      IF (KT.EQ.2) GO TO 470
C     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
      Y(NN) = TEMP(1)
      Y(NN-1) = TEMP(2)
      IF (NN.EQ.2) RETURN
      TRX = 2.0E0/X
      DTM = FNI
      TM = (DTM+FNF)*TRX
      K = NN + 1
      DO 460 I=3,NN
        K = K - 1
        Y(K-2) = TM*Y(K-1) - Y(K)
        DTM = DTM - 1.0E0
        TM = (DTM+FNF)*TRX
  460 CONTINUE
      RETURN
  470 Y(1) = TEMP(2)
      RETURN
C
C     ASYMPTOTIC EXPANSION FOR X TO INFINITY WITH FORWARD RECURSION IN
C     OSCILLATORY REGION X.GT.MAX(20, NU), PROVIDED THE LAST MEMBER
C     OF THE SEQUENCE IS ALSO IN THE REGION.
C
  480 CONTINUE
      IN = INT(ALPHA-TAU+2.0E0)
      IF (IN.LE.0) GO TO 490
      IDALP = IALP - IN - 1
      KT = 1
      GO TO 500
  490 CONTINUE
      IDALP = IALP
      IN = 0
  500 IS = KT
      FIDAL = FLOAT(IDALP)
      DALPHA = FIDAL + FNF
      ARG = X - PIDT*DALPHA - PDF
      SA = SIN(ARG)
      SB = COS(ARG)
      COEF = RTTP/RTX
      ETX = 8.0E0*X
  510 CONTINUE
      DTM = FIDAL + FIDAL
      DTM = DTM*DTM
      TM = 0.0E0
      IF (FIDAL.EQ.0.0E0 .AND. ABS(FNF).LT.TOL) GO TO 520
      TM = 4.0E0*FNF*(FIDAL+FIDAL+FNF)
  520 CONTINUE
      TRX = DTM - 1.0E0
      T2 = (TRX+TM)/ETX
      S2 = T2
      RELB = TOL*ABS(T2)
      T1 = ETX
      S1 = 1.0E0
      FN = 1.0E0
      AK = 8.0E0
      DO 530 K=1,13
        T1 = T1 + ETX
        FN = FN + AK
        TRX = DTM - FN
        AP = TRX + TM
        T2 = -T2*AP/T1
        S1 = S1 + T2
        T1 = T1 + ETX
        AK = AK + 8.0E0
        FN = FN + AK
        TRX = DTM - FN
        AP = TRX + TM
        T2 = T2*AP/T1
        S2 = S2 + T2
        IF (ABS(T2).LE.RELB) GO TO 540
        AK = AK + 8.0E0
  530 CONTINUE
  540 TEMP(IS) = COEF*(S1*SB-S2*SA)
      IF(IS.EQ.2) GO TO 560
  550 FIDAL = FIDAL + 1.0E0
      DALPHA = FIDAL + FNF
      IS = 2
      TB = SA
      SA = -SB
      SB = TB
      GO TO 510
C
C     FORWARD RECURSION SECTION
C
  560 IF (KT.EQ.2) GO TO 470
      S1 = TEMP(1)
      S2 = TEMP(2)
      TX = 2.0E0/X
      TM = DALPHA*TX
      IF (IN.EQ.0) GO TO 580
C
C     FORWARD RECUR TO INDEX ALPHA
C
      DO 570 I=1,IN
        S = S2
        S2 = TM*S2 - S1
        TM = TM + TX
        S1 = S
  570 CONTINUE
      IF (NN.EQ.1) GO TO 600
      S = S2
      S2 = TM*S2 - S1
      TM = TM + TX
      S1 = S
  580 CONTINUE
C
C     FORWARD RECUR FROM INDEX ALPHA TO ALPHA+N-1
C
      Y(1) = S1
      Y(2) = S2
      IF (NN.EQ.2) RETURN
      DO 590 I=3,NN
        Y(I) = TM*Y(I-1) - Y(I-2)
        TM = TM + TX
  590 CONTINUE
      RETURN
  600 Y(1) = S2
      RETURN
C
C     BACKWARD RECURSION WITH NORMALIZATION BY
C     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.
C
  610 CONTINUE
C     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
      AKM = AMAX1(3.0E0-FN,0.0E0)
      KM = INT(AKM)
      TFN = FN + FLOAT(KM)
      TA = (GLN+TFN-0.9189385332E0-0.0833333333E0/TFN)/(TFN+0.5E0)
      TA = XO2L - TA
      TB = -(1.0E0-1.5E0/TFN)/TFN
      AKM = TOLLN/(-TA+SQRT(TA*TA-TOLLN*TB)) + 1.5E0
      IN = KM + INT(AKM)
      GO TO 660
  620 CONTINUE
C     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
      GLN = WK(3) + WK(2)
      IF (WK(6).GT.30.0E0) GO TO 640
      RDEN = (PP(4)*WK(6)+PP(3))*WK(6) + 1.0E0
      RZDEN = PP(1) + PP(2)*WK(6)
      TA = RZDEN/RDEN
      IF (WK(1).LT.0.10E0) GO TO 630
      TB = GLN/WK(5)
      GO TO 650
  630 TB=(1.259921049E0+(0.1679894730E0+0.0887944358E0*WK(1))*WK(1))
     1 /WK(7)
      GO TO 650
  640 CONTINUE
      TA = 0.5E0*TOLLN/WK(4)
      TA=((0.0493827160E0*TA-0.1111111111E0)*TA+0.6666666667E0)*TA*WK(6)
      IF (WK(1).LT.0.10E0) GO TO 630
      TB = GLN/WK(5)
  650 IN = INT(TA/TB+1.5E0)
      IF (IN.GT.INLIM) GO TO 310
  660 CONTINUE
      DTM = FNI + FLOAT(IN)
      TRX = 2.0E0/X
      TM = (DTM+FNF)*TRX
      TA = 0.0E0
      TB = TOL
      KK = 1
  670 CONTINUE
C
C     BACKWARD RECUR UNINDEXED
C
      DO 680 I=1,IN
        S = TB
        TB = TM*TB - TA
        TA = S
        DTM = DTM - 1.0E0
        TM = (DTM+FNF)*TRX
  680 CONTINUE
C     NORMALIZATION
      IF (KK.NE.1) GO TO 690
      TA = (TA/TB)*TEMP(3)
      TB = TEMP(3)
      KK = 2
      IN = NS
      IF (NS.NE.0) GO TO 670
  690 Y(NN) = TB
      NZ = N - NN
      IF (NN.EQ.1) RETURN
      K = NN - 1
      Y(K) = TM*TB - TA
      IF (NN.EQ.2) RETURN
      DTM = DTM - 1.0E0
      TM = (DTM+FNF)*TRX
      KM = K - 1
C
C     BACKWARD RECUR INDEXED
C
      DO 700 I=1,KM
        Y(K-1) = TM*Y(K) - Y(K+1)
        DTM = DTM - 1.0E0
        TM = (DTM+FNF)*TRX
        K = K - 1
  700 CONTINUE
      RETURN
C
C
C
  710 CONTINUE
      CALL XERROR( 'BESJ - ORDER, ALPHA, LESS THAN ZERO.', 36, 2, 1)
      RETURN
  720 CONTINUE
      CALL XERROR( 'BESJ - N LESS THAN ONE.', 23, 2, 1)
      RETURN
  730 CONTINUE
      CALL XERROR( 'BESJ - X LESS THAN ZERO.', 24, 2, 1)
      RETURN
      END
