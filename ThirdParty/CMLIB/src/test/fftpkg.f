C
C   DRIVER FOR TESTING CMLIB ROUTINES
C      COSQB    COSQF    COSQI    COST     COSTI    EZFFTB
C      EZFFTF   RFFTB    RFFTF    RFFTI    SINQB    SINQF
C      SINQI    SINT     SINTI
C
C    ONE INPUT DATA CARD IS REQUIRED
C         READ(LIN,1) KPRINT,TIMES
C    1    FORMAT(I1,E10.0)
C
C     KPRINT = 0   NO PRINTING
C              1   NO PRINTING FOR PASSED TESTS, SHORT MESSAGE
C                  FOR FAILED TESTS
C              2   PRINT SHORT MESSAGE FOR PASSED TESTS, FULLER
C                  INFORMATION FOR FAILED TESTS
C              3   PRINT COMPLETE QUICK-CHECK RESULTS
C
C                ***IMPORTANT NOTE***
C         ALL QUICK CHECKS USE ROUTINES R2MACH AND D2MACH
C         TO SET THE ERROR TOLERANCES.
C     TIMES IS A CONSTANT MULTIPLIER THAT CAN BE USED TO SCALE THE
C     VALUES OF R1MACH AND D1MACH SO THAT
C               R2MACH(I) = R1MACH(I) * TIMES   FOR I=3,4,5
C               D2MACH(I) = D1MACH(I) * TIMES   FOR I=3,4,5
C     THIS MAKES IT EASILY POSSIBLE TO CHANGE THE ERROR TOLERANCES
C     USED IN THE QUICK CHECKS.
C     IF TIMES .LE. 0.0 THEN TIMES IS DEFAULTED TO 1.0
C
C              ***END NOTE***
C
      COMMON/UNIT/LUN
      COMMON/MSG/ICNT,JTEST(38)
      COMMON/XXMULT/TIMES
      LUN=I1MACH(2)
      LIN=I1MACH(1)
      ITEST=1
C
C     READ KPRINT,TIMES PARAMETERS FROM DATA CARD..
C
      READ(LIN,1) KPRINT,TIMES
1     FORMAT(I1,E10.0)
      IF(TIMES.LE.0.) TIMES=1.
      CALL XSETUN(LUN)
      CALL XSETF(1)
      CALL XERMAX(1000)
C   TEST FFT PACKAGE
      CALL FFTQX(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C
      IF(KPRINT.GE.1.AND.ITEST.NE.1) WRITE(LUN,2)
2     FORMAT(/' ***** WARNING -- AT LEAST ONE TEST FOR SUBLIBRARY,
     1 FFTPKG HAS FAILED ***** ')
      IF(KPRINT.GE.1.AND.ITEST.EQ.1) WRITE(LUN,3)
3     FORMAT(/' ----- SUBLIBRARY FFTPKG PASSED ALL TESTS ----- ')
      END
      DOUBLE PRECISION FUNCTION D2MACH(I)
      DOUBLE PRECISION D1MACH
      COMMON/XXMULT/TIMES
      D2MACH=D1MACH(I)
      IF(I.EQ.1.OR. I.EQ.2) RETURN
      D2MACH = D2MACH * DBLE(TIMES)
      RETURN
      END
      SUBROUTINE FFTQX(KPRINT,IPASS)
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C                       VERSION 3  JUNE 1979
C
C                         A TEST DRIVER FOR
C          A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE FAST FOURIER
C           TRANSFORM OF PERIODIC AND OTHER SYMMETRIC SEQUENCES
C
C                              BY
C
C                       PAUL N SWARZTRAUBER
C
C       NATIONAL CENTER FOR ATMOSPHERIC RESEARCH  BOULDER,COLORADO 80307
C
C        WHICH IS SPONSORED BY THE NATIONAL SCIENCE FOUNDATION
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C             THIS PROGRAM TESTS THE PACKAGE OF FAST FOURIER
C     TRANSFORMS FOR BOTH COMPLEX AND REAL PERIODIC SEQUENCES AND
C     CERTIAN OTHER SYMMETRIC SEQUENCES THAT ARE LISTED BELOW.
C
C     1.   RFFTI     INITIALIZE  RFFTF AND RFFTB
C     2.   RFFTF     FORWARD TRANSFORM OF A REAL PERIODIC SEQUENCE
C     3.   RFFTB     BACKWARD TRANSFORM OF A REAL COEFFICIENT ARRAY
C
C     4.   EZFFTF    A SIMPLIFIED REAL PERIODIC FORWARD TRANSFORM
C     5.   EZFFTB    A SIMPLIFIED REAL PERIODIC BACKWARD TRANSFORM
C
C     6.   SINTI     INITIALIZE SINT
C     7.   SINT      SINE TRANSFORM OF A REAL ODD SEQUENCE
C
C     8.   COSTI     INITIALIZE COST
C     9.   COST      COSINE TRANSFORM OF A REAL EVEN SEQUENCE
C
C     10.  SINQI     INITIALIZE SINQF AND SINQB
C     11.  SINQF     FORWARD SINE TRANSFORM WITH ODD WAVE NUMBERS
C     12.  SINQB     UNNORMALIZED INVERSE OF SINQF
C
C     13.  COSQI     INITIALIZE COSQF AND COSQB
C     14.  COSQF     FORWARD COSINE TRANSFORM WITH ODD WAVE NUMBERS
C     15.  COSQB     UNNORMALIZED INVERSE OF COSQF
C
C     16.  CFFTI     INITIALIZE CFFTF AND CFFTB
C     17.  CFFTF     FORWARD TRANSFORM OF A COMPLEX PERIODIC SEQUENCE
C     18.  CFFTB     UNNORMALIZED INVERSE OF CFFTF
C
C
      COMMON/MSG/ICNT,ITEST
      INTEGER ITEST(38)
      COMMON/UNIT/LUN
      DIMENSION       ND(10)     ,X(200)     ,Y(200)     ,W(2000)    ,
     1                A(100)     ,B(100)     ,AH(100)    ,BH(100)    ,
     2                XH(200)
      DATA ND(1),ND(2),ND(3),ND(4),ND(5)/120,32,54,4,2/
      NNS = 5
      ITOP=17
      ERRMAX=SQRT(R1MACH(3))
      IF(KPRINT.GE.2)WRITE(LUN,999)
999   FORMAT('1','FFT QUICK CHECK')
      IPASS=1
      DO 150 NZ=1,NNS
      DO 100 I=1,ITOP
100   ITEST(I)=0
      ICNT=1
         N = ND(NZ)
         FN = FLOAT(N)
         TFN = FN+FN
         NP1 = N+1
         NM1 = N-1
         DO 101 J=1,N
            X(J) = FLOAT(J)/FN
            Y(J) = X(J)
  101    CONTINUE
C
C     TEST SUBROUTINES RFFTI,RFFTF AND RFFTB
C
         CALL RFFTI (N,W)
         PI = 3.14159265358979
         DT = (PI+PI)/FN
         NS2 = N/2
         CF = 1.
         IF (NS2 .LT. 2) GO TO 104
         DO 103 K=2,NS2
            SUM1 = 0.
            SUM2 = 0.
            ARG = FLOAT(K-1)*DT
            DO 102 I=1,N
               ARG1 = FLOAT(I-1)*ARG
               SUM1 = SUM1+X(I)*COS(ARG1)
               SUM2 = SUM2-X(I)*SIN(ARG1)
  102       CONTINUE
            Y(2*K-2) = CF*SUM1
            Y(2*K-1) = CF*SUM2
  103    CONTINUE
  104    SUM1 = 0.
         SUM2 = 0.
         DO 105 I=1,NM1,2
            SUM1 = SUM1+X(I)
            SUM2 = SUM2+X(I+1)
  105    CONTINUE
         Y(1) = CF*(SUM1+SUM2)
         Y(N) = CF*(SUM1-SUM2)
         CALL RFFTF (N,X,W)
         RFTF = 0.
         DO 106 I=1,N
            X(I) = X(I)/TFN
            Y(I) = Y(I)/TFN
            RFTF = AMAX1(RFTF,ABS(X(I)-Y(I)))
  106    CONTINUE
      IF(RFTF.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
         NHM = N/2-1
        CF=2.
         DO 109 I=1,N
            SUM = X(1)+FLOAT((-1)**(I-1))*X(N)
            IF (NHM .EQ. 0) GO TO 108
            ARG = FLOAT(I-1)*DT
            DO 107 K=2,NS2
               ARG1 = FLOAT(K-1)*ARG
               SUM = SUM+CF*X(2*K-2)*COS(ARG1)-CF*X(2*K-1)*SIN(ARG1)
  107       CONTINUE
  108       Y(I) = SUM
  109    CONTINUE
         CALL RFFTB (N,X,W)
         RFTB = 0.
         DO 110 I=1,N
            RFTB = AMAX1(RFTB,ABS(X(I)-Y(I)))
            X(I) = Y(I)
  110    CONTINUE
      IF(RFTB.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
         CALL RFFTB (N,Y,W)
         CALL RFFTF (N,Y,W)
         CF = 1./(FN+FN)
         CF=2.0*CF
         RFTFB = 0.
         DO 111 I=1,N
            RFTFB = AMAX1(RFTFB,ABS(CF*Y(I)-X(I)))
  111    CONTINUE
      IF(RFTFB.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
C
C     TEST SUBROUTINES SINTI AND SINT
C
         DT = PI/FN
         DO 112 I=1,NM1
            X(I) = FLOAT(I)/FN
  112    CONTINUE
         DO 114 I=1,NM1
            Y(I) = 0.
            ARG1 = FLOAT(I)*DT
            DO 113 K=1,NM1
               Y(I) = Y(I)+X(K)*SIN(FLOAT(K)*ARG1)
  113       CONTINUE
            Y(I) = 2.*Y(I)
  114    CONTINUE
         CALL SINTI (NM1,W)
         CALL SINT (NM1,X,W)
         CF = .125/FN
         SINTT = 0.
         DO 115 I=1,NM1
            X(I) = CF*X(I)
            Y(I) = CF*Y(I)
            SINTT = AMAX1(SINTT,ABS(X(I)-Y(I)))
            X(I) = FLOAT(I)/FN
            Y(I) = X(I)
  115    CONTINUE
      IF(SINTT.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
         CALL SINT (NM1,X,W)
         CALL SINT (NM1,X,W)
         SINTFB = 0.
       CF=1./(2.*N)
         DO 116 I=1,NM1
            SINTFB = AMAX1(SINTFB,ABS(CF*X(I)-Y(I)))
  116    CONTINUE
      IF(SINTFB.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
C
C     TEST SUBROUTINES COSTI AND COST
C
         DO 117 I=1,NP1
            X(I) = FLOAT(I)/FN
  117    CONTINUE
         DO 119 I=1,NP1
            Y(I) = .5*(X(1)+FLOAT((-1)**(I+1))*X(N+1))
            ARG = FLOAT(I-1)*DT
            DO 118 K=2,N
               Y(I) = Y(I)+X(K)*COS(FLOAT(K-1)*ARG)
  118       CONTINUE
            Y(I) = 2.*Y(I)
  119    CONTINUE
         CALL COSTI (NP1,W)
         CALL COST (NP1,X,W)
         COSTT = 0.
         DO 120 I=1,NP1
            X(I) = CF*X(I)
            Y(I) = CF*Y(I)
            COSTT = AMAX1(COSTT,ABS(X(I)-Y(I)))
            X(I) = Y(I)
  120    CONTINUE
      IF(COSTT.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
         CALL COST (NP1,X,W)
         CALL COST (NP1,X,W)
         COSTFB = 0.
      CF=1./(2.*N)
         DO 121 I=1,NP1
            COSTFB = AMAX1(COSTFB,ABS(CF*X(I)-Y(I)))
  121    CONTINUE
      IF(COSTFB.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
C
C     TEST SUBROUTINES SINQI,SINQF AND SINQB
C
         DO 122 I=1,N
            Y(I) = SIN(FLOAT(I))
  122    CONTINUE
         DT = PI/(FN+FN)
         DO 124 I=1,N
            X(I) = 0.
            ARG = DT*FLOAT(I)
            DO 123 K=1,N
               X(I) = X(I)+Y(K)*SIN(FLOAT(K+K-1)*ARG)
  123       CONTINUE
            X(I) = 4.*X(I)
  124    CONTINUE
         CALL SINQI (N,W)
         CALL SINQB (N,Y,W)
         SINQBT = 0.
         DO 125 I=1,N
            X(I) = CF*X(I)
            Y(I) = CF*Y(I)
            SINQBT = AMAX1(SINQBT,ABS(Y(I)-X(I)))
  125    CONTINUE
      IF(SINQBT.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
         DO 127 I=1,N
            ARG = FLOAT(I+I-1)*DT
            Y(I) = FLOAT((-1)**(I+1))*X(N)
            DO 126 K=1,NM1
               Y(I) = Y(I)+2.*X(K)*SIN(FLOAT(K)*ARG)
  126       CONTINUE
  127    CONTINUE
         CALL SINQF (N,X,W)
         SINQFT = 0.
         DO 128 I=1,N
            SINQFT = AMAX1(SINQFT,ABS(X(I)-Y(I)))
            Y(I) = X(I)
  128    CONTINUE
      IF(SINQFT.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
      CF=.5*CF
         CALL SINQF (N,Y,W)
         CALL SINQB (N,Y,W)
         SINQFB = 0.
         DO 129 I=1,N
            SINQFB = AMAX1(SINQFB,ABS(CF*Y(I)-X(I)))
  129    CONTINUE
      IF(SINQFB.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
C
C     TEST SUBROUTINES COSQI,COSQF AND COSQB
C
         DO 130 I=1,N
            Y(I) = SIN(FLOAT(I))
  130    CONTINUE
         DO 132 I=1,N
            X(I) = 0.
            ARG = FLOAT(I-1)*DT
            DO 131 K=1,N
               X(I) = X(I)+Y(K)*COS(FLOAT(K+K-1)*ARG)
  131       CONTINUE
            X(I) = 4.*X(I)
  132    CONTINUE
         CALL COSQI (N,W)
         CALL COSQB (N,Y,W)
         COSQBT = 0.
         DO 133 I=1,N
            X(I) = CF*X(I)
            Y(I) = CF*Y(I)
            COSQBT = AMAX1(COSQBT,ABS(X(I)-Y(I)))
  133    CONTINUE
      IF(COSQBT.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
         DO 135 I=1,N
            Y(I) = X(1)
            ARG = FLOAT(I+I-1)*DT
            DO 134 K=2,N
               Y(I) = Y(I)+2.*X(K)*COS(FLOAT(K-1)*ARG)
  134       CONTINUE
  135    CONTINUE
         CALL COSQF (N,X,W)
         COSQFT = 0.
         DO 136 I=1,N
            COSQFT = AMAX1(COSQFT,ABS(Y(I)-X(I)))
            Y(I) = X(I)
  136    CONTINUE
      IF(COSQFT.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
         CALL COSQB (N,X,W)
         CALL COSQF (N,X,W)
         COSQFB = 0.
         DO 137 I=1,N
            COSQFB = AMAX1(COSQFB,ABS(CF*X(I)-Y(I)))
  137    CONTINUE
      IF(COSQFB.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
C
C     TEST PROGRAMS EZFFTF,EZFFTB
C
      CALL EZFFTI(N,W)
         DO 138 I=1,N
            X(I) = SIN(FLOAT(I))
  138    CONTINUE
         TPI = 8.*ATAN(1.)
         DT = TPI/FLOAT(N)
         NS2 = N/2
         CF = 2./FLOAT(N)
         NS2M = NS2-1
         IF (NS2M .LE. 0) GO TO 141
         DO 140 K=1,NS2M
            SUM1 = 0.
            SUM2 = 0.
            ARG = FLOAT(K)*DT
            DO 139 I=1,N
               ARG1 = FLOAT(I-1)*ARG
               SUM1 = SUM1+X(I)*COS(ARG1)
               SUM2 = SUM2+X(I)*SIN(ARG1)
  139       CONTINUE
            A(K) = CF*SUM1
            B(K) = CF*SUM2
  140    CONTINUE
  141    NM1 = N-1
         SUM1 = 0.
         SUM2 = 0.
         DO 142 I=1,NM1,2
            SUM1 = SUM1+X(I)
            SUM2 = SUM2+X(I+1)
  142    CONTINUE
         AZERO = .5*CF*(SUM1+SUM2)
         A(NS2) = .5*CF*(SUM1-SUM2)
         CALL EZFFTF (N,X,AZEROH,AH,BH,W)
         DEZF1 = AMAX1(ABS(AZEROH-AZERO),ABS(A(NS2)-AH(NS2)))
         DO 143 I=1,NS2M
            DEZF1 = AMAX1(DEZF1,ABS(AH(I)-A(I)),ABS(BH(I)-B(I)))
  143    CONTINUE
      IF(DEZF1.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
         CALL EZFFTB (N,XH,AZEROH,AH,BH,W)
         DEZB1 = 0.
         DO 144 I=1,N
            DEZB1 = AMAX1(DEZB1,ABS(XH(I)-X(I)))
  144    CONTINUE
      IF(DEZB1.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
C
C     TEST FOR N ODD
C
         NP1 = N+1
         CALL EZFFTI(NP1,W)
         DT = TPI/FLOAT(NP1)
         NM1S2 = (NP1-1)/2
         DO 145 I=1,NM1S2
            A(I) = SIN(FLOAT(I))
            B(I) = COS(FLOAT(I))
  145    CONTINUE
         AZERO = .25
         DO 147 I=1,NP1
            SUM = AZERO
            ARG = FLOAT(I-1)*DT
            DO 146 K=1,NM1S2
               ARG1 = FLOAT(K)*ARG
               SUM = SUM+A(K)*COS(ARG1)+B(K)*SIN(ARG1)
  146       CONTINUE
            X(I) = SUM
  147    CONTINUE
         CALL EZFFTB (NP1,XH,AZERO,A,B,W)
         DEZB2 = 0.
         DO 148 I=1,NP1
            DEZB2 = AMAX1(DEZB2,ABS(X(I)-XH(I)))
  148    CONTINUE
      IF(DEZB2.LE.10.*ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      ICNT=ICNT+1
         CALL EZFFTF (NP1,XH,AZEROH,AH,BH,W)
         DEZF2 = ABS(AZEROH-AZERO)
         DO 149 I=1,NM1S2
            DEZF2 = AMAX1(DEZF2,ABS(AH(I)-A(I)),ABS(BH(I)-B(I)))
  149    CONTINUE
      IF(DEZF2.LE.ERRMAX) ITEST(ICNT)=1
C     IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.ITEST(ICNT).NE.1))CALL PASS
      IPSS=1
      DO 1490 I=1,ITOP
1490  IPSS=IPSS*ITEST(I)
      ICNT=NZ
      IF(KPRINT.GE.2.OR.(KPRINT.EQ.1.AND.IPSS.NE.1))CALL PASS
      IPASS=IPSS*IPASS
         IF(KPRINT.EQ.3)
     1   WRITE (LUN,1001) N,RFTF,RFTB,RFTFB,SINTT,SINTFB,COSTT,COSTFB,
     2                  SINQFT,SINQBT,SINQFB,COSQFT,COSQBT,COSQFB,DEZF1,
     3                  DEZB1,DEZF2,DEZB2
  150 CONTINUE
      IF(KPRINT.GE.1.AND.IPASS.NE.1) WRITE(LUN,99998)
      IF(KPRINT.GE.2.AND.IPASS.EQ.1) WRITE(LUN,99999)
99999 FORMAT(/' ***** FFT ROUTINES PASSED ALL TESTS ***** ')
99998 FORMAT(/' ***** FFT ROUTINES FAILED SOME TESTS ***** ')
      RETURN
C
C
 1001 FORMAT ('0N',I5,' RFFTF  ',E10.3,' RFFTB  ',E10.3,' RFFTFB ',
     1        E10.3,' SINT   ',E10.3,' SINTFB ',E10.3,' COST   ',E10.3/
     2        7X,' COSTFB ',E10.3,' SINQF  ',E10.3,' SINQB  ',E10.3,
     3        ' SINQFB ',E10.3,' COSQF  ',E10.3,' COSQB  ',E10.3/7X,
     4        ' COSQFB ',E10.3,' DEZF1  ',E10.3,' DEZB1  ',E10.3,
     5        ' DEZF2  ',E10.3,' DEZB2  ',E10.3)
C
      END
      SUBROUTINE PASS
      COMMON /UNIT/ LUN
      COMMON /MSG/ ICNT,ITEST
      DIMENSION ITEST(38)
      IF(ITEST(ICNT).EQ.0) GO TO 10
      WRITE(LUN,100) ICNT
  100 FORMAT(/'TEST NUMBER',I5,'  PASSED')
      RETURN
   10 WRITE(LUN,200) ICNT
  200 FORMAT(/' *****TEST NUMBER',I5,'  FAILED***** ')
      RETURN
      END
      REAL FUNCTION R2MACH(I)
      COMMON/XXMULT/TIMES
      R2MACH=R1MACH(I)
      IF(I.EQ.1.OR. I.EQ.2) RETURN
      R2MACH = R2MACH * TIMES
      RETURN
      END
