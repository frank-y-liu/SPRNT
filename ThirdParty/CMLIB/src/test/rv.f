C
C   DRIVER FOR TESTING CMLIB ROUTINES
C     RNOR    UNI
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
C   TEST  RNOR
      CALL RVQX1(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C   TEST UNI
      CALL RVQX2(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C
      IF(KPRINT.GE.1.AND.ITEST.NE.1) WRITE(LUN,2)
2     FORMAT(/' ***** WARNING -- AT LEAST ONE TEST FOR SUBLIBRARY RV
     1  HAS FAILED ***** ')
      IF(KPRINT.GE.1.AND.ITEST.EQ.1) WRITE(LUN,3)
3     FORMAT(/' ----- SUBLIBRARY RV PASSED ALL TESTS ----- ')
      END
      DOUBLE PRECISION FUNCTION D2MACH(I)
      DOUBLE PRECISION D1MACH
      COMMON/XXMULT/TIMES
      D2MACH=D1MACH(I)
      IF(I.EQ.1.OR. I.EQ.2) RETURN
      D2MACH = D2MACH * DBLE(TIMES)
      RETURN
      END
      REAL FUNCTION R2MACH(I)
      COMMON/XXMULT/TIMES
      R2MACH=R1MACH(I)
      IF(I.EQ.1.OR. I.EQ.2) RETURN
      R2MACH = R2MACH * TIMES
      RETURN
      END
      SUBROUTINE RVQX1(KPRINT,IPASS)
C              Test of RNOR    (Normal (0,1) generator)
      COMMON/UNIT/LUN
      EXTERNAL DERF
      REAL H(20)
      INTEGER KPRINT,IPASS
      DOUBLE PRECISION DERF, XL,XU
      IPASS=1
      DUM=RNOR(87)
      IF(KPRINT.GT.2) WRITE(LUN,*) DUM
      A=-8.
      B=8.0
      M=20
      N=20000
       DO 9 IND=1,M
    9   H(IND)=0.0
      DO 1 I=1,N
      VAL=(RNOR(0))
      IF(I.LE.5 .AND. KPRINT.GT.2) WRITE(LUN,*) VAL
      IF(I.EQ.999 .AND. KPRINT.GT.2) WRITE(LUN,*) VAL
      CALL HIST(VAL,H,M,A,B)
    1 CONTINUE
      DO 2 I=1,M
        H(I)=H(I)/FLOAT(N)
   2  CONTINUE
      IF(KPRINT.EQ.3) WRITE(LUN,*) (H(I),I=1,M)
      CHISQ=0.
      DO 3 I=1,M
      XL=(A+(I-1.)*(B-A)/FLOAT(M))/SQRT(2.D0)
      XU=(A+I*(B-A)/FLOAT(M))/SQRT(2.D0)
      IF(XU.GE.0. .AND. XL.GE.0.)THEN
         P=DERF(XU)-DERF(XL)
      ELSEIF(XU.LE.0. .AND. XL.LE.0.)THEN
         P=DERF(-XL)-DERF(-XU)
      ELSE
         P=DERF(XU)+DERF(-XL)
      ENDIF
      P=P/2.
      CHISQ=CHISQ+(H(I)-P)**2/P
      IF(KPRINT.EQ.3) WRITE(LUN,*)  P,(H(I)-P)**2*N/P,H(I)
   3  CONTINUE
      CHISQ=CHISQ*N
      IF(KPRINT.EQ.3) WRITE(LUN,*) ' CHI SQUARE=',CHISQ
      IF(CHISQ.GT.14.0)THEN
        IF(KPRINT.GE.1)WRITE(LUN,*)'* TEST FOR SUBROUTINE RNOR FAILED *'
        IPASS=2
      ELSE
        IF(KPRINT.GE.2)WRITE(LUN,*)'* TEST FOR SUBROUTINE RNOR PASSED *'
      ENDIF
   10 CONTINUE
      RETURN
      END
      SUBROUTINE HIST(VAL,H,N,A,B)
      REAL H(N)
      HL=(B-A)/FLOAT(N)
      DO 1 I=1,N
        IF(VAL.GE.A+(I-1.)*HL .AND. VAL.LT.A+I*HL) THEN
          H(I)=H(I)+1.
          GO TO 2
        ENDIF
   1  CONTINUE
    2 CONTINUE
      RETURN
      END
      SUBROUTINE RVQX2(KPRINT,IPASS)
C                   TEST OF UNI (RANDOM NUMBER GENERATOR)
      DIMENSION X(100)
      COMMON/UNIT/LUN
      INTEGER KPRINT,IPASS
      IPASS=1
      DO 2 I=1,100
    2 X(I)=0.
      JD=305
      R=UNI(305)
      DO 1 I=1,10000
      R1=UNI(0)
      CALL PART(K,10,R1)
      X(K)=X(K)+1
    1 CONTINUE
      IF(KPRINT.GT.2) WRITE(LUN,100)  (X(I),I=1,10)
      CHI=0.0
      DO 3 I=1,10
         CHI=CHI+(X(I)-1000.)**2
    3 CONTINUE
      CHI=10./10000.*CHI
      IF(KPRINT.GT.2) WRITE(LUN,*) 'CHI SQUARE =',CHI
      IF(CHI.LT.21.7 .AND. KPRINT.GE.2) WRITE(LUN,888)
      IF(CHI.GE.21.7) THEN
       WRITE(LUN,889)
       IPASS=2
      ENDIF
  100 FORMAT(E16.8)
  888 FORMAT (1X,' ** TEST FOR SUBROUTINE UNI PASSED ** ')
  889 FORMAT (1X,' ** TEST FOR SUBROUTINE UNI FAILED ** ')
      RETURN
      END
      SUBROUTINE PART(I,N,X)
      DO 1 J=1,N
      I=J
      IF((I-1)/FLOAT(N).LE.X .AND. I/FLOAT(N).GT.X) RETURN
    1 CONTINUE
      RETURN
      END