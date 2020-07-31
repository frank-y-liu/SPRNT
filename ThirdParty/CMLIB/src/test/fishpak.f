C
C   DRIVER FOR TESTING CMLIB ROUTINES
C      HWSCRT
C      HWSPLR
C      HWSCYL
C      HWSSSP
C      HWSCSP
C      GENBUN
C      BLKTRI
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
      COMMON WORK(10000)
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
C  TEST HWSCRT
      CALL QXCRT(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C  TEST HWSPLR
      CALL QXPLR(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C  TEST HWSCYL
      CALL QXCYL(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C  TEST HWSSSP
      CALL QXSSP(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C  TEST HWSCSP
      CALL QXCSP(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C  TEST GENBUN
      CALL QXGBUN(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C  TEST BLKTRI
      CALL QXBLKT(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C
      IF(KPRINT.GE.1.AND.ITEST.NE.1) WRITE(LUN,2)
2     FORMAT(/' ***** WARNING -- AT LEAST ONE TEST FOR SUBLIBRARY,
     1 FISHPAK HAS FAILED ***** ')
      IF(KPRINT.GE.1.AND.ITEST.EQ.1) WRITE(LUN,3)
3     FORMAT(/54H------------- FISHPAK PASSED ALL TESTS----------------)
      END
      SUBROUTINE QXBLKT(KPRINT,IPASS)
C FISHPAK39  FROM PORTLIB                                  01/03/80
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                        F I S H P A K                          *
C     *                                                               *
C     *                                                               *
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
C     *                                                               *
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
C     *                                                               *
C     *                  (VERSION  3 , JUNE 1979)                     *
C     *                                                               *
C     *                             BY                                *
C     *                                                               *
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
C     *                                                               *
C     *                             OF                                *
C     *                                                               *
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
C     *                                                               *
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *
C     *                                                               *
C     *                   WHICH IS SPONSORED BY                       *
C     *                                                               *
C     *              THE NATIONAL SCIENCE FOUNDATION                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C     PROGRAM TO ILLUSTRATE THE USE OF BLKTRI
C
      COMMON/UNIT/LUN
      COMMON/MSG/ICNT,ITEST(38)
      COMMON       Y(75,105)  ,AM(75)     ,BM(75)     ,CM(75)     ,
     1                AN(105)    ,BN(105)    ,CN(105)    ,W(1952)    ,
     2                S(75)      ,T(105)
C
      ERMAX=1.E-3
      IFLG = 0
      NP = 1
      N = 63
      MP = 1
      M = 50
      IDIMY = 75
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING THE
C     COEFFICIENTS AND THE ARRAY Y.
C
      DELTAS = 1./FLOAT(M+1)
      DO 101 I=1,M
         S(I) = FLOAT(I)*DELTAS
  101 CONTINUE
      DELTAT = 1./FLOAT(N+1)
      DO 102 J=1,N
         T(J) = FLOAT(J)*DELTAT
  102 CONTINUE
C
C     COMPUTE THE COEFFICIENTS AM,BM,CM CORRESPONDING TO THE S DIRECTION
C
      HDS = DELTAS/2.
      TDS = DELTAS+DELTAS
      DO 103 I=1,M
         TEMP1 = 1./(S(I)*TDS)
         TEMP2 = 1./((S(I)-HDS)*TDS)
         TEMP3 = 1./((S(I)+HDS)*TDS)
         AM(I) = TEMP1*TEMP2
         CM(I) = TEMP1*TEMP3
         BM(I) = -(AM(I)+CM(I))
  103 CONTINUE
C
C     COMPUTE THE COEFFICIENTS AN,BN,CN CORRESPONDING TO THE T DIRECTION
C
      HDT = DELTAT/2.
      TDT = DELTAT+DELTAT
      DO 104 J=1,N
         TEMP1 = 1./(T(J)*TDT)
         TEMP2 = 1./((T(J)-HDT)*TDT)
         TEMP3 = 1./((T(J)+HDT)*TDT)
         AN(J) = TEMP1*TEMP2
         CN(J) = TEMP1*TEMP3
         BN(J) = -(AN(J)+CN(J))
  104 CONTINUE
C
C     COMPUTE RIGHT SIDE OF EQUATION
C
      DO 106 J=1,N
         DO 105 I=1,M
            Y(I,J) = 3.75*S(I)*T(J)*(S(I)**4.+T(J)**4.)
  105    CONTINUE
  106 CONTINUE
C
C     INCLUDE NONHOMOGENEOUS BOUNDARY INTO RIGHT SIDE. NOTE THAT THE
C     CORNER AT J=N,I=M INCLUDES CONTRIBUTIONS FROM BOTH BOUNDARIES.
C
      DO 107 J=1,N
         Y(M,J) = Y(M,J)-CM(M)*T(J)**5.
  107 CONTINUE
      DO 108 I=1,M
         Y(I,N) = Y(I,N)-CN(N)*S(I)**5.
  108 CONTINUE
C
  109 CALL BLKTRI (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,IERROR,W)
      IFLG = IFLG+1
      IF (IFLG-1) 109,109,110
C
C     COMPUTE DISCRETIZATION ERROR
C
  110 ERR = 0.
      DO 112 J=1,N
         DO 111 I=1,M
            Z = ABS(Y(I,J)-(S(I)*T(J))**5.)
            IF (Z .GT. ERR) ERR = Z
  111    CONTINUE
  112 CONTINUE
      IPASS=1
      IF(ERR.GT.ERMAX)IPASS=0
      IF(KPRINT.EQ.0)RETURN
      IF(KPRINT.GE.3 .OR. (KPRINT.EQ.2.AND.IPASS.EQ.0))
     +  WRITE(LUN,1001) IERROR,ERR,W(1)
      RETURN
C
 1001 FORMAT ('1',20X,'SUBROUTINE BLKTRI EXAMPLE'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/
     3        18X,'DISCRETIZATION ERROR = 1.6478E-05'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 823'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/18X,'DISCRETIZATION ERROR =',E12.5/
     7        12X,'REQUIRED LENGTH OF W ARRAY =', F4.0)
C
      END
      SUBROUTINE QXCRT(KPRINT,IPASS)
C FISHPAK25  FROM PORTLIB                                  01/03/80
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                        F I S H P A K                          *
C     *                                                               *
C     *                                                               *
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
C     *                                                               *
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
C     *                                                               *
C     *                  (VERSION  3 , JUNE 1979)                     *
C     *                                                               *
C     *                             BY                                *
C     *                                                               *
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
C     *                                                               *
C     *                             OF                                *
C     *                                                               *
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
C     *                                                               *
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *
C     *                                                               *
C     *                   WHICH IS SPONSORED BY                       *
C     *                                                               *
C     *              THE NATIONAL SCIENCE FOUNDATION                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C          PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HWSCRT TO SOLVE
C     THE EQUATION
C
C     (D/DX)(DU/DX) + (D/DY)(DU/DY) - 4*U
C
C     = (2 - (4 + PI**2/4)*X**2)*COS((Y+1)*PI/2)
C
C     WITH THE BOUNDARY CONDITIONS
C     ON THE RECTANGLE 0 .LT. X .LT. 2, -1 .LT. Y .LT. 3 WITH THE
C
C     U(0,Y) = 0
C                                          -1 .LE. Y .LE. 3
C     (DU/DX)(2,Y) = 4*COS((Y+1)*PI/2)
C
C     AND WITH U PERIODIC IN Y.
C          THE X-INTERVAL WILL BE DIVIDED INTO 40 PANELS AND THE
C     Y-INTERVAL WILL BE DIVIDED INTO 80 PANELS.
C
      COMMON /UNIT/LUN
      COMMON/MSG/ICNT,ITEST(38)
      COMMON       F(45,82)   ,BDB(81)    ,W(1200)    ,X(41)      ,
     1                Y(81)
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.  ALSO NOTE THAT W
C     IS DIMENSIONED 6*(N+1) + 8*(M+1).
C
      IDIMF = 45
      ERMAX=1.E-3
      A = 0.
      B = 2.
      M = 40
      MBDCND = 2
      C = -1.
      D = 3.
      N = 80
      NBDCND = 0
      ELMBDA = -4.
C
C     AUXILIARY QUANTITIES.
C
      PI = ACOS(-1.)
      PIBY2 = PI/2.
      PISQ = PI**2
      MP1 = M+1
      NP1 = N+1
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE HELMHOLTZ EQUATION.
C
      DO 101 I=1,MP1
         X(I) = FLOAT(I-1)/20.
  101 CONTINUE
      DO 102 J=1,NP1
         Y(J) = -1.+FLOAT(J-1)/20.
  102 CONTINUE
C
C     GENERATE BOUNDARY DATA.
C
      DO 103 J=1,NP1
         BDB(J) = 4.*COS((Y(J)+1.)*PIBY2)
  103 CONTINUE
C
C     BDA, BDC, AND BDD ARE DUMMY VARIABLES.
C
      DO 104 J=1,NP1
         F(1,J) = 0.
  104 CONTINUE
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO 106 I=2,MP1
         DO 105 J=1,NP1
            F(I,J) = (2.-(4.+PISQ/4.)*X(I)**2)*COS((Y(J)+1.)*PIBY2)
  105    CONTINUE
  106 CONTINUE
      CALL HWSCRT(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,ELMBDA,F,
     1             IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C                U(X,Y) = X**2*COS((Y+1)*PIBY2)
C
      ERR = 0.
      DO 108 I=1,MP1
         DO 107 J=1,NP1
            Z = ABS(F(I,J)-X(I)**2*COS((Y(J)+1.)*PIBY2))
            IF (Z .GT. ERR) ERR = Z
  107    CONTINUE
  108 CONTINUE
      IPASS=1
      IF(ERR.GT.ERMAX)IPASS=0
      IF(KPRINT.EQ.0)RETURN
      IF(KPRINT.GE.3 .OR. (KPRINT.EQ.2.AND.IPASS.EQ.0))
     +  WRITE(LUN,1001) IERROR,ERR,W(1)
      RETURN
C
 1001 FORMAT ('1',20X,'SUBROUTINE HWSCRT EXAMPLE'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/
     3        18X,'DISCRETIZATION ERROR = 5.36508E-04'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 880'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/18X,'DISCRETIZATION ERROR =',E12.5/
     7        12X,'REQUIRED LENGTH OF W ARRAY =',F4.0)
C
      END
      SUBROUTINE QXCSP(KPRINT,IPASS)
C FISHPAK29  FROM PORTLIB                                  01/03/80
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                        F I S H P A K                          *
C     *                                                               *
C     *                                                               *
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
C     *                                                               *
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
C     *                                                               *
C     *                  (VERSION  3 , JUNE 1979)                     *
C     *                                                               *
C     *                             BY                                *
C     *                                                               *
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
C     *                                                               *
C     *                             OF                                *
C     *                                                               *
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
C     *                                                               *
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *
C     *                                                               *
C     *                   WHICH IS SPONSORED BY                       *
C     *                                                               *
C     *              THE NATIONAL SCIENCE FOUNDATION                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C     PROGRAM TO ILLUSTRATE THE USE OF HWSCSP
C
      COMMON/UNIT/LUN
      COMMON/MSG/ICNT,ITEST(38)
      COMMON       F(48,33)   ,BDTF(33)   ,W(1200)     ,R(33)      ,
     1                THETA(48)
C
C     THE VALUE OF IDIMF IS THE FIRST DIMENSION OF F.  SINCE M=36, N=32,
C     L=N THEREFORE K=5 AND W(.) IS DIMENSIONED 2*(L+1)*(K-1) + 6*(M+N)
C     + MAX(4*N,6*M) + 14 = 902.
C
      ERMAX=1.E-3
      PI = ACOS(-1.)
      INTL = 0
      TS = 0.
      TF = PI/2.
      M = 36
      MBDCND = 6
      RS = 0.
      RF = 1.
      N = 32
      NBDCND = 5
      ELMBDA = 0.
      IDIMF = 48
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING THE
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE EQUATION.
C
      MP1 = M+1
      DTHETA = TF/FLOAT(M)
      DO 101 I=1,MP1
         THETA(I) = FLOAT(I-1)*DTHETA
  101 CONTINUE
      NP1 = N+1
      DR = 1./FLOAT(N)
      DO 102 J=1,NP1
         R(J) = FLOAT(J-1)*DR
  102 CONTINUE
C
C     GENERATE NORMAL DERIVATIVE DATA AT EQUATOR
C
      DO 103 J=1,NP1
         BDTF(J) = 0.
  103 CONTINUE
C
C     COMPUTE BOUNDARY DATA ON THE SURFACE OF THE SPHERE
C
      DO 104 I=1,MP1
         F(I,N+1) = COS(THETA(I))**4
  104 CONTINUE
C
C     COMPUTE RIGHT SIDE OF EQUATION
C
      DO 106 I=1,MP1
         CI4 = 12.*COS(THETA(I))**2
         DO 105 J=1,N
            F(I,J) = CI4*R(J)**2
  105    CONTINUE
  106 CONTINUE
C
      CALL HWSCSP (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,BDRS,
     1             BDRF,ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR
C
      ERR = 0.
      DO 108 I=1,MP1
         CI4 = COS(THETA(I))**4
         DO 107 J=1,N
            Z = ABS(F(I,J)-CI4*R(J)**4)
            IF (Z .GT. ERR) ERR = Z
  107    CONTINUE
  108 CONTINUE
      IW = INT(W(1))
      IPASS=1
      IF(ERR.GT.ERMAX)IPASS=0
      IF(KPRINT.EQ.0)GO TO 99999
      IF(KPRINT.GE.3 .OR. (KPRINT.EQ.2.AND.IPASS.EQ.0))
     +  WRITE(LUN,1001) IERROR,ERR,IW
C
C     THE FOLLOWING PROGRAM ILLUSTRATES THE USE OF HWSCSP TO SOLVE
C     A THREE DIMENSIONAL PROBLEM WHICH HAS LONGITUDNAL DEPENDENCE
C
99999 MBDCND = 2
      NBDCND = 1
      DPHI = PI/72.
      ELMBDA = -2.*(1.-COS(DPHI))/DPHI**2
C
C     COMPUTE BOUNDARY DATA ON THE SURFACE OF THE SPHERE
C
      DO 109 I=1,MP1
         F(I,N+1) = SIN(THETA(I))
  109 CONTINUE
C
C     COMPUTE RIGHT SIDE OF THE EQUATION
C
      DO 111 J=1,N
         DO 110 I=1,MP1
            F(I,J) = 0.
  110    CONTINUE
  111 CONTINUE
C
      CALL HWSCSP (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,BDRS,
     1             BDRF,ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR   (FOURIER COEFFICIENTS)
C
      ERR = 0
      DO 113 I=1,MP1
         SI = SIN(THETA(I))
         DO 112 J=1,NP1
            Z = ABS(F(I,J)-R(J)*SI)
            IF (Z .GT. ERR) ERR = Z
  112    CONTINUE
  113 CONTINUE
C
      IW = INT(W(1))
      IF(ERR.GT.ERMAX)IPASS=0
      IF(KPRINT.EQ.0)RETURN
      IF(KPRINT.GE.3 .OR. (KPRINT.EQ.2.AND.IPASS.EQ.0))
     +  WRITE(LUN,1002) IERROR,ERR,IW
      RETURN
C
 1001 FORMAT ('1',20X,'SUBROUTINE HWSCSP EXAMPLE 1'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/
     3        18X,'DISCRETIZATION ERROR = 7.99842E-04'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 775'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/18X,'DISCRETIZATION ERROR =',E12.5/
     7        12X,'REQUIRED LENGTH OF W ARRAY =',I4)
 1002 FORMAT ('1',20X,'SUBROUTINE HWSCSP EXAMPLE 2'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/
     3        18X,'DISCRETIZATION ERROR = 5.86824E-05'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 775'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/18X,'DISCRETIZATION ERROR =',E12.5/
     7        12X,'REQUIRED LENGTH OF W ARRAY =',I4)
C
      END
      SUBROUTINE QXCYL(KPRINT,IPASS)
C FISHPAK27  FROM PORTLIB                                  01/03/80
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                        F I S H P A K                          *
C     *                                                               *
C     *                                                               *
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
C     *                                                               *
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
C     *                                                               *
C     *                  (VERSION  3 , JUNE 1979)                     *
C     *                                                               *
C     *                             BY                                *
C     *                                                               *
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
C     *                                                               *
C     *                             OF                                *
C     *                                                               *
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
C     *                                                               *
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *
C     *                                                               *
C     *                   WHICH IS SPONSORED BY                       *
C     *                                                               *
C     *              THE NATIONAL SCIENCE FOUNDATION                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C          PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HWSCYL TO SOLVE
C     THE EQUATION
C
C     (1/R)(D/DR)(R*(DU/DR)) + (D/DZ)(DU/DZ)
C
C     = (2*R*Z)**2*(4*Z**2 + 3*R**2)
C
C     ON THE RECTANGLE 0 .LT. R .LT. 1, 0 .LT. Z .LT. 1 WITH THE
C     BOUNDARY CONDITIONS
C
C     U(0,Z) UNSPECIFIED
C                                            0 .LE. Z .LE. 1
C     (DU/DR)(1,Z) = 4*Z**4
C
C     AND
C
C     (DU/DZ)(R,0) = 0
C                                            0 .LE. R .LE. 1
C     (DU/DZ)(R,1) = 4*R**4 .
C
C          THE R-INTERVAL WILL BE DIVIDED INTO 50 PANELS AND THE
C     Z-INTERVAL WILL BE DIVIDED INTO 100 PANELS.
C
      COMMON/UNIT/LUN
      COMMON/MSG/ICNT,ITEST(38)
      COMMON       F(75,105)  ,BDA(101)   ,BDB(101)   ,BDC(51)    ,
     1                BDD(51)    ,W(1200)    ,R(51)      ,Z(101)
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.  ALSO NOTE THAT W
C     IS DIMENSIONED 6*(N+1) + 8*(M+1).
C
      IDIMF = 75
      ERMAX=1.5E-3
      A = 0.
      B = 1.
      M = 50
      MBDCND = 6
      C = 0.
      D = 1.
      N = 100
      NBDCND = 3
      ELMBDA = 0.
C
C     AUXILIARY QUANTITIES.
C
      MP1 = M+1
      NP1 = N+1
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE POISSON EQUATION.
C
      DO 101 I=1,MP1
         R(I) = FLOAT(I-1)/50.
  101 CONTINUE
      DO 102 J=1,NP1
         Z(J) = FLOAT(J-1)/100.
  102 CONTINUE
C
C     GENERATE BOUNDARY DATA.
C
      DO 103 J=1,NP1
         BDB(J) = 4.*Z(J)**4
  103 CONTINUE
      DO 104 I=1,MP1
         BDC(I) = 0.
         BDD(I) = 4.*R(I)**4
  104 CONTINUE
C
C     BDA IS A DUMMY VARIABLE.
C
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO 106 I=1,MP1
         DO 105 J=1,NP1
            F(I,J) = 4.*R(I)**2*Z(J)**2*(4.*Z(J)**2+3.*R(I)**2)
  105    CONTINUE
  106 CONTINUE
      CALL HWSCYL(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,ELMBDA,F,
     1             IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR BY MINIMIZING OVER ALL A THE FUNCTION
C     NORM(F(I,J) - A*1 - U(R(I),Z(J))).  THE EXACT SOLUTION IS
C                U(R,Z) = (R*Z)**4 + ARBITRARY CONSTANT.
C
      X = 0.
      DO 108 I=1,MP1
         DO 107 J=1,NP1
            X = X+F(I,J)-(R(I)*Z(J))**4
  107    CONTINUE
  108 CONTINUE
      X = X/FLOAT(NP1*MP1)
      DO 110 I=1,MP1
         DO 109 J=1,NP1
            F(I,J) = F(I,J)-X
  109    CONTINUE
  110 CONTINUE
      ERR = 0.
      DO 112 I=1,MP1
         DO 111 J=1,NP1
            X = ABS(F(I,J)-(R(I)*Z(J))**4)
            IF (X .GT. ERR) ERR = X
  111    CONTINUE
  112 CONTINUE
      IPASS=1
      IF(ERR.GT.ERMAX)IPASS=0
      IF(KPRINT.EQ.0)RETURN
      IF(KPRINT.GE.3 .OR. (KPRINT.EQ.2.AND.IPASS.EQ.0))
     +  WRITE(LUN,1001) IERROR,PERTRB,ERR,W(1)
      RETURN
C
 1001 FORMAT ('1',20X,'SUBROUTINE HWSCYL EXAMPLE'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/32X,'PERTRB = 2.26734E-04'/
     3        18X,'DISCRETIZATION ERROR = 3.73672E-04'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 1118'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/32X,'PERTRB =',E12.5/
     7        18X,'DISCRETIZATION ERROR =',E12.5/
     8        12X,'REQUIRED LENGTH OF W ARRAY =',F5.0)
C
      END
      SUBROUTINE QXGBUN(KPRINT,IPASS)
C FISHPAK00  FROM PORTLIB                                  07/75
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                        F I S H P A K                          *
C     *                                                               *
C     *                                                               *
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
C     *                                                               *
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
C     *                                                               *
C     *                  (VERSION  3 , JUNE 1979)                     *
C     *                                                               *
C     *                             BY                                *
C     *                                                               *
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
C     *                                                               *
C     *                             OF                                *
C     *                                                               *
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
C     *                                                               *
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *
C     *                                                               *
C     *                   WHICH IS SPONSORED BY                       *
C     *                                                               *
C     *              THE NATIONAL SCIENCE FOUNDATION                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C     PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE GENBUN
C
      COMMON/UNIT/LUN
      COMMON/MSG/ICNT,ITEST(38)
      COMMON F(25,130), A(20), B(20), C(20), W(1200), X(20), Y(120)
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMY.  ALSO NOTE THAT
C     W(.) IS DIMENSIONED 6*N + 5*M.
C
      ERMAX=1.E-2
      IDIMY = 25
      IFLG = 0
      MPEROD = 1
      M = 20
      DELTAX = 1./FLOAT(M)
      NPEROD = 0
      N = 120
      PI = ACOS(-1.)
      DELTAY = 2.*PI/FLOAT(N)
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
C     COEFFICIENTS AND RIGHT SIDE OF EQUATION.
C
      DO 100 I=1,M
         X(I) = FLOAT(I-1)*DELTAX
  100    CONTINUE
      DO 105 J=1,N
         Y(J) = -PI + FLOAT(J-1)*DELTAY
  105    CONTINUE
C
C     GENERATE COEFFICIENTS.
C
      S = (DELTAY/DELTAX)**2
      T = S*DELTAX
      A(1) = 0.
      B(1) = -2.*S
      C(1) = 2.*S
      DO 110 I=2,M
         A(I) = (1.+X(I))**2*S + (1.+X(I))*T
         C(I) = (1.+X(I))**2*S - (1.+X(I))*T
         B(I) = -2.*(1.+X(I))**2*S
  110    CONTINUE
      C(M) = 0.
C
C     GENERATE RIGHT SIDE OF EQUATION FOR I = 1 SHOWING INTRODUCTION OF
C     BOUNDARY DATA.
C
      DYSQ = DELTAY**2
      DO 115 J=1,N
         F(1,J) = DYSQ*(11. + 8./DELTAX)*SIN(Y(J))
  115    CONTINUE
C
C     GENERATE RIGHT SIDE.
C
      MM1 = M-1
      DO 125 I=2,MM1
         DO 120 J=1,N
            F(I,J) = DYSQ*3.*(1.+X(I))**4*SIN(Y(J))
  120       CONTINUE
  125    CONTINUE
C
C     GENERATE RIGHT SIDE FOR I = M SHOWING INTRODUCTION OF
C     BOUNDARY DATA.
C
      DO 130 J=1,N
         F(M,J) = DYSQ*(3.*(1.+X(M))**4 - 16.*((1.+X(M))/DELTAX)**2
     +                + 16.*(1.+X(M))/DELTAX)*SIN(Y(J))
  130    CONTINUE
      CALL GENBUN(NPEROD,N,MPEROD,M,A,B,C,IDIMY,F,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C                   U(X,Y) = (1+X)**4*SIN(Y)
C
      ERR = 0.
      DO 140 I=1,M
         DO 135 J=1,N
            Z = ABS(F(I,J)-(1.+X(I))**4*SIN(Y(J)))
            IF (Z .GT. ERR) ERR = Z
  135       CONTINUE
  140    CONTINUE
      IW=INT(W(1))
      IPASS=1
      IF(ERR.GT.ERMAX)IPASS=0
      IF(KPRINT.EQ.0)RETURN
      IF(KPRINT.GE.3 .OR. (KPRINT.EQ.2.AND.IPASS.EQ.0))
     +  WRITE(LUN,1001) IERROR,ERR,IW
      RETURN
C
 1001 FORMAT ('1',20X,'SUBROUTINE GENBUN EXAMPLE'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/
     3        18X,'DISCRETIZATION ERROR = 7.94113E-03'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 820'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/18X,'DISCRETIZATION ERROR =',E12.5/
     7        12X,'REQUIRED LENGTH OF W ARRAY =',I4)
C
      END
      SUBROUTINE QXPLR(KPRINT,IPASS)
C FISHPAK26  FROM PORTLIB                                  01/03/80
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                        F I S H P A K                          *
C     *                                                               *
C     *                                                               *
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
C     *                                                               *
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
C     *                                                               *
C     *                  (VERSION  3 , JUNE 1979)                     *
C     *                                                               *
C     *                             BY                                *
C     *                                                               *
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
C     *                                                               *
C     *                             OF                                *
C     *                                                               *
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
C     *                                                               *
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *
C     *                                                               *
C     *                   WHICH IS SPONSORED BY                       *
C     *                                                               *
C     *              THE NATIONAL SCIENCE FOUNDATION                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C          PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HWSPLR TO SOLVE
C     THE EQUATION
C
C     (1/R)(D/DR)(R*(DU/DR)) + (1/R**2)(D/DTHETA)(DU/DTHETA) = 16*R**2
C
C     ON THE QUARTER-DISK 0 .LT. R .LT. 1, 0 .LT. THETA .LT. PI/2 WITH
C     WITH THE BOUNDARY CONDITIONS
C
C     U(1,THETA) = 1 - COS(4*THETA), 0 .LE. THETA .LE. 1
C
C     AND
C
C     (DU/DTHETA)(R,0) = (DU/DTHETA)(R,PI/2) = 0,  0 .LE. R .LE. 1.
C
C     (NOTE THAT THE SOLUTION U IS UNSPECIFIED AT R = 0.)
C          THE R-INTERVAL WILL BE DIVIDED INTO 50 PANELS AND THE
C     THETA-INTERVAL WILL BE DIVIDED INTO 48 PANELS.
C
      COMMON/UNIT/LUN
      COMMON/MSG/ICNT,ITEST(38)
      COMMON       F(100,50)  ,BDC(51)    ,BDD(51)    ,W(1200)    ,
     1                R(51)      ,THETA(49)
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.  ALSO NOTE THAT W
C     IS DIMENSIONED 6*(N+1) + 8*(M+1).
C
      IDIMF = 100
      ERMAX=1.E-3
      A = 0.
      B = 1.
      M = 50
      MBDCND = 5
      C = 0.
      PI = ACOS(-1.)
      D = PI/2.
      N = 48
      NBDCND = 3
      ELMBDA = 0.
C
C     AUXILIARY QUANTITIES.
C
      MP1 = M+1
      NP1 = N+1
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE POISSON EQUATION.
C
      DO 101 I=1,MP1
         R(I) = FLOAT(I-1)/50.
  101 CONTINUE
      DO 102 J=1,NP1
         THETA(J) = FLOAT(J-1)*PI/96.
  102 CONTINUE
C
C     GENERATE BOUNDARY DATA.
C
      DO 103 I=1,MP1
         BDC(I) = 0.
         BDD(I) = 0.
  103 CONTINUE
C
C     BDA AND BDB ARE DUMMY VARIABLES.
C
      DO 104 J=1,NP1
         F(MP1,J) = 1.-COS(4.*THETA(J))
  104 CONTINUE
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO 106 I=1,M
         DO 105 J=1,NP1
            F(I,J) = 16.*R(I)**2
  105    CONTINUE
  106 CONTINUE
      CALL HWSPLR(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,ELMBDA,F,
     1             IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C                U(R,THETA) = R**4*(1 - COS(4*THETA))
C
      ERR = 0.
      DO 108 I=1,MP1
         DO 107 J=1,NP1
            Z = ABS(F(I,J)-R(I)**4*(1.-COS(4.*THETA(J))))
            IF (Z .GT. ERR) ERR = Z
  107    CONTINUE
  108 CONTINUE
      IPASS=1
      IF(ERR.GT.ERMAX)IPASS=0
      IF(KPRINT.EQ.0)RETURN
      IF(KPRINT.GE.3 .OR. (KPRINT.EQ.2.AND.IPASS.EQ.0))
     +  WRITE(LUN,1001) IERROR,ERR,W(1)
      RETURN
C
 1001 FORMAT ('1',20X,'SUBROUTINE HWSPLR EXAMPLE'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/
     3        18X,'DISCRETIZATION ERROR = 6.19134E-04'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 882'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/18X,'DISCRETIZATION ERROR =',E12.5/
     7        12X,'REQUIRED LENGTH OF W ARRAY =',F4.0)
C
      END
      SUBROUTINE QXSSP(KPRINT,IPASS)
C FISHPAK28  FROM PORTLIB                                  01/03/80
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                        F I S H P A K                          *
C     *                                                               *
C     *                                                               *
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
C     *                                                               *
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
C     *                                                               *
C     *                  (VERSION  3 , JUNE 1979)                     *
C     *                                                               *
C     *                             BY                                *
C     *                                                               *
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
C     *                                                               *
C     *                             OF                                *
C     *                                                               *
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
C     *                                                               *
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *
C     *                                                               *
C     *                   WHICH IS SPONSORED BY                       *
C     *                                                               *
C     *              THE NATIONAL SCIENCE FOUNDATION                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C
C     PROGRAM TO ILLUSTRATE THE USE OF HWSSSP
C
      COMMON/UNIT/LUN
      COMMON/MSG/ICNT,ITEST(38)
      COMMON       F(19,73)   ,BDTF(73)   ,SINT(19)   ,SINP(73)   ,
     1                W(1200)
C
C     THE VALUE OF IDIMF IS THE FIRST DIMENSION OF F.  W IS
C     DIMENSIONED 11*(M+1)+6*(N+1)=647 SINCE M=18 AND N=72.
C
      PI = ACOS(-1.)
      ERMAX=5.E-3
      TS = 0
      TF = PI/2.
      M = 18
      MBDCND = 6
      PS = 0
      PF = PI+PI
      N = 72
      NBDCND = 0
      ELMBDA = 0.
      IDIMF = 19
C
C     GENERATE SINES FOR USE IN SUBSEQUENT COMPUTATIONS
C
      DTHETA = TF/FLOAT(M)
      MP1 = M+1
      DO 101 I=1,MP1
         SINT(I) = SIN(FLOAT(I-1)*DTHETA)
  101 CONTINUE
      DPHI = (PI+PI)/FLOAT(N)
      NP1 = N+1
      DO 102 J=1,NP1
         SINP(J) = SIN(FLOAT(J-1)*DPHI)
  102 CONTINUE
C
C     COMPUTE RIGHT SIDE OF EQUATION AND STORE IN F
C
      DO 104 J=1,NP1
         DO 103 I=1,MP1
            F(I,J) = 2.-6.*(SINT(I)*SINP(J))**2
  103    CONTINUE
  104 CONTINUE
C
C     STORE DERIVATIVE DATA AT THE EQUATOR
C
      DO 105 J=1,NP1
         BDTF(J) = 0.
  105 CONTINUE
C
      CALL HWSSSP(TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,N,NBDCND,BDPS,BDPF,
     1             ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR. SINCE PROBLEM IS SINGULAR, THE
C     SOLUTION MUST BE NORMALIZED.
C
      ERR = 0
      DO 107 J=1,NP1
         DO 106 I=1,MP1
            Z = ABS(F(I,J)-(SINT(I)*SINP(J))**2-F(1,1))
            IF (Z .GT. ERR) ERR = Z
  106    CONTINUE
  107 CONTINUE
C
      IW = INT(W(1))
      IPASS=1
      IF(ERR.GT.ERMAX)IPASS=0
      IF(KPRINT.EQ.0)RETURN
      IF(KPRINT.GE.3 .OR. (KPRINT.EQ.2.AND.IPASS.EQ.0))
     +  WRITE(LUN,1001) IERROR,ERR,IW
      RETURN
C
C
 1001 FORMAT ('1',20X,'SUBROUTINE HWSSSP EXAMPLE'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/
     3        18X,'DISCRETIZATION ERROR = 3.38107E-03'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 600'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/18X,'DISCRETIZATION ERROR =',E12.5 /
     7        12X,'REQUIRED LENGTH OF W ARRAY =',I4)
C
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
