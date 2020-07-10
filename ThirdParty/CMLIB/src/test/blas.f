C
C   DRIVER FOR TESTING CMLIB ROUTINES
C      "BLAS SUBPROGRAMS"
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
C   TEST BLAS
      CALL BLACHK(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C
      IF(KPRINT.GE.1.AND.ITEST.NE.1) WRITE(LUN,2)
2     FORMAT(/'   ********** WARNING -- AT LEAST ONE TEST FOR THE BLAS,
     1 SUBLIBRARY  HAS FAILED ****************** ')
      IF(KPRINT.GE.1.AND.ITEST.EQ.1) WRITE(LUN,3)
3     FORMAT(/' ----- THE BLAS SUBLIBRARY PASSED ALL TESTS ----- ')
      END
      SUBROUTINE BLACHK (KPRINT,IPASS)
C1    ********************************* TBLA ***************************
C     TEST DRIVER FOR BASIC LINEAR ALGEBRA SUBPROGRAMS.
C     C. L. LAWSON,JPL, 1974 DEC 10, 1975 MAY 28
C2
C     MODIFIED FOR SANDIA MATH LIBRARY USE BY K. HASKELL, 6/23/77
C     UPDATED BY K. HASKELL - JUNE 23,1980
C
      COMMON /UNIT/ LUN
      COMMON/COMBLA/NPRINT,ICASE,N,INCX,INCY,MODE,PASS
      COMMON /MSG/ ICNT,JTEST
      DIMENSION JTEST(38)
      LOGICAL          PASS
      INTEGER          ITEST(38)
      DOUBLE PRECISION DFAC,DQFAC
      DATA SFAC,SDFAC,DFAC,DQFAC / .625E-1, .50, .625D-1, 0.625D-1/
      DATA ITEST /1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,
     1            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/
C     THE ZEROS IN THE ABOVE DATA STATEMENT ARE TO SUPPRESS
C     TESTING OF DQDOTI AND DQDOTA, WHICH DO NOT EXIST IN
C     NONTRIVIAL SUBROUTINES ON THE SANDIA MATH. LIBRARY.
      NPRINT = LUN
      ICNT=0
C
C
    5 CONTINUE
      IF (KPRINT.GE.2) WRITE (NPRINT,1005)
 1005 FORMAT(/' QUICK CHECK OF 38 BASIC LINEAR ALGEBRA SUBROUTINES  '/)
          DO 60 ICASE=1,38
          IF(ITEST(ICASE) .EQ. 0) GO TO 60
          ICNT = ICNT+1
          CALL HEADER (KPRINT)
C
C         INITIALIZE  PASS, INCX, INCY, AND MODE FOR A NEW CASE.
C         THE VALUE 9999 FOR INCX, INCY OR MODE WILL APPEAR IN THE
C         DETAILED  OUTPUT, IF ANY, FOR CASES THAT DO NOT INVOLVE
C         THESE PARAMETERS.
C
          PASS=.TRUE.
          INCX=9999
          INCY=9999
          MODE=9999
              GO TO (12,12,12,12,12,12,12,12,12,12,
     A               12,10,10,12,12,10,10,12,12,12,
     B               12,12,12,12,12,11,11,11,11,11,
     C               11,11,11,11,11,11,11,11),  ICASE
C                                       ICASE = 12-13 OR 16-17
   10         CALL CHECK0(SFAC,DFAC,KPRINT)
              GO TO 50
C                                       ICASE = 26-38
   11         CALL CHECK1(SFAC,DFAC,KPRINT)
              GO TO 50
C                                       ICASE =  1-11, 14-15, OR 18-25
   12         CALL CHECK2(SFAC,SDFAC,DFAC,DQFAC,KPRINT)
   50         CONTINUE
C                                                  PRINT
          IF (KPRINT.GE.2 .AND. PASS) WRITE (NPRINT,1001)
      JTEST(ICNT) = 1
      IF (.NOT.PASS) JTEST(ICNT) = 0
   60     CONTINUE
      IPASS=1
      DO 70 I = 1, ICNT
      IPASS = IPASS*JTEST(I)
   70 CONTINUE
      IF (KPRINT.GE.2 .AND. IPASS.EQ.1) WRITE (NPRINT,1006)
      IF (KPRINT.GE.1 .AND. IPASS.EQ.0) WRITE (NPRINT,1007)
      RETURN
 1001 FORMAT('+                                       PASS')
 1006 FORMAT(/'   ****************BLAS PASSED ALL TESTS***************')
 1007 FORMAT(/'   ****************BLAS FAILED SOME TESTS**************')
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
      SUBROUTINE CHECK0(SFAC,DFAC,KPRINT)
C1    ********************************* CHECK0 *************************
C     THIS SUBROUTINE TESTS SUBPROGRAMS 12-13 AND 16-17.
C     THESE SUBPROGRAMS HAVE NO ARRAY ARGUMENTS.
C
C     C. L. LAWSON, JPL, 1975 MAR 07, MAY 28
C     R. J. HANSON, J. A. WISNIEWSKI, SANDIA LABS, APRIL 25,1977.
C2
      COMMON/COMBLA/NPRINT,ICASE,N,INCX,INCY,MODE,PASS
      LOGICAL          PASS
      REAL             STRUE(9),STEMP(9)
      DOUBLE PRECISION DC,DS,DA1(8),DB1(8),DC1(8),DS1(8)
      DOUBLE PRECISION DA,DATRUE(8),DBTRUE(8),DZERO,DFAC,DB
      DOUBLE PRECISION DAB(4,9),DTEMP(9),DTRUE(9,9),D12
      DATA ZERO, DZERO / 0., 0.D0 /
      DATA DA1/ .3D0,  .4D0, -.3D0, -.4D0, -.3D0,  0.D0,  0.D0,  1.D0/
      DATA DB1/ .4D0,  .3D0,  .4D0,  .3D0, -.4D0,  0.D0,  1.D0,  0.D0/
      DATA DC1/ .6D0,  .8D0, -.6D0,  .8D0,  .6D0,  1.D0,  0.D0,  1.D0/
      DATA DS1/ .8D0,  .6D0,  .8D0, -.6D0,  .8D0,  0.D0,  1.D0,  0.D0/
      DATA DATRUE/ .5D0,  .5D0,  .5D0, -.5D0, -.5D0, 0.D0, 1.D0, 1.D0/
      DATA DBTRUE/ 0.D0,  .6D0,  0.D0, -.6D0,  0.D0, 0.D0, 1.D0, 0.D0/
C                                              INPUT FOR MODIFIED GIVENS
      DATA DAB/ .1D0,.3D0,1.2D0,.2D0,
     A          .7D0, .2D0, .6D0, 4.2D0,
     B          0.D0,0.D0,0.D0,0.D0,
     C          4.D0, -1.D0, 2.D0, 4.D0,
     D          6.D-10, 2.D-2, 1.D5, 10.D0,
     E          4.D10, 2.D-2, 1.D-5, 10.D0,
     F          2.D-10, 4.D-2, 1.D5, 10.D0,
     G          2.D10, 4.D-2, 1.D-5, 10.D0,
     H          4.D0, -2.D0, 8.D0, 4.D0    /
C                                       TRUE RESULTS FOR MODIFIED GIVENS
      DATA DTRUE/0.D0,0.D0, 1.3D0, .2D0, 0.D0,0.D0,0.D0, .5D0, 0.D0,
     A           0.D0,0.D0, 4.5D0, 4.2D0, 1.D0, .5D0, 0.D0,0.D0,0.D0,
     B           0.D0,0.D0,0.D0,0.D0, -2.D0, 0.D0,0.D0,0.D0,0.D0,
     C           0.D0,0.D0,0.D0, 4.D0, -1.D0, 0.D0,0.D0,0.D0,0.D0,
     D           0.D0, 15.D-3, 0.D0, 10.D0, -1.D0, 0.D0, -1.D-4,
     E           0.D0, 1.D0,
     F           0.D0,0.D0, 6144.D-5, 10.D0, -1.D0, 4096.D0, -1.D6,
     G           0.D0, 1.D0,
     H           0.D0,0.D0,15.D0,10.D0,-1.D0, 5.D-5, 0.D0,1.D0,0.D0,
     I           0.D0,0.D0, 15.D0, 10.D0, -1. D0, 5.D5, -4096.D0,
     J           1.D0, 4096.D-6,
     K           0.D0,0.D0, 7.D0, 4.D0, 0.D0,0.D0, -.5D0, -.25D0, 0.D0/
C                   4096 = 2 ** 12
      DATA D12  /4096.D0/
C
C                   COMPUTE TRUE VALUES WHICH CANNOT BE PRESTORED
C                   IN DECIMAL NOTATION.
      DTRUE(1,1) = 12.D0 / 130.D0
      DTRUE(2,1) = 36.D0 / 130.D0
      DTRUE(7,1) = -1.D0 / 6.D0
      DTRUE(1,2) = 14.D0 / 75.D0
      DTRUE(2,2) = 49.D0 / 75.D0
      DTRUE(9,2) = 1.D0 / 7.D0
      DTRUE(1,5) = 45.D-11 * (D12 * D12)
      DTRUE(3,5) = 4.D5 / (3.D0 * D12)
      DTRUE(6,5) = 1.D0 / D12
      DTRUE(8,5) = 1.D4 / (3.D0 * D12)
      DTRUE(1,6) = 4.D10 / (1.5D0 * D12 * D12)
      DTRUE(2,6) = 2.D-2 / 1.5D0
      DTRUE(8,6) = 5.D-7 * D12
      DTRUE(1,7) = 4.D0 / 150.D0
      DTRUE(2,7) = (2.D-10 / 1.5D0) * (D12 * D12)
      DTRUE(7,7) = -DTRUE(6,5)
      DTRUE(9,7) = 1.D4 / D12
      DTRUE(1,8) = DTRUE(1,7)
      DTRUE(2,8) = 2.D10 / (1.5D0 * D12 * D12)
      DTRUE(1,9) = 32.D0 / 7.D0
      DTRUE(2,9) = -16.D0 / 7.D0
      DBTRUE(1) = 1.D0/.6D0
      DBTRUE(3) = -1.D0/.6D0
      DBTRUE(5) = 1.D0/.6D0
C
      JUMP= ICASE-11
          DO 500 K = 1, 9
C                        SET N=K FOR IDENTIFICATION IN OUTPUT IF ANY.
          N=K
C                             BRANCH TO SELECT SUBPROGRAM TO BE TESTED.
C
          GO TO (120,130,999,999,160,170), JUMP
C                                                             12. SROTG
  120 IF(K.GT.8) GO TO 600
          SA=SNGL(DA1(K))
          SB = SNGL(DB1(K))
          CALL SROTG(SA,SB,SC,SS)
          CALL STEST(1,SA,SNGL(DATRUE(K)),SNGL(DATRUE(K)),SFAC,KPRINT)
          CALL STEST(1,SB,SNGL(DBTRUE(K)),SNGL(DBTRUE(K)),SFAC,KPRINT)
          CALL STEST(1,SC,SNGL(DC1(K)),SNGL(DC1(K)),SFAC,KPRINT)
          CALL STEST(1,SS,SNGL(DS1(K)),SNGL(DS1(K)),SFAC,KPRINT)
          GO TO 500
C                                                             13. DROTG
  130 IF(K.GT.8) GO TO 600
          DA=DA1(K)
          DB = DB1(K)
          CALL DROTG(DA,DB,DC,DS)
          CALL DTEST(1,DA,DATRUE(K),DATRUE(K),DFAC,KPRINT)
          CALL DTEST(1,DB,DBTRUE(K),DBTRUE(K),DFAC,KPRINT)
          CALL DTEST(1,DC,DC1(K),DC1(K),DFAC,KPRINT)
          CALL DTEST(1,DS,DS1(K),DS1(K),DFAC,KPRINT)
          GO TO 500
C                                                             16. SROTMG
  160     CONTINUE
               DO 162 I = 1, 4
               STEMP(I) = SNGL(DAB(I,K))
               STEMP(I+4) = ZERO
  162          CONTINUE
           STEMP(9) = ZERO
           CALL SROTMG(STEMP(1),STEMP(2),STEMP(3),STEMP(4),STEMP(5))
C
               DO 166 I = 1, 9
  166          STRUE(I) = SNGL(DTRUE(I,K))
          CALL STEST(9,STEMP,STRUE,STRUE,SFAC,KPRINT)
          GO TO 500
C                                                             17. DROTMG
  170     CONTINUE
               DO 172 I = 1, 4
               DTEMP(I) = DAB(I,K)
               DTEMP(I+4) = DZERO
  172          CONTINUE
          DTEMP(9) = DZERO
          CALL DROTMG(DTEMP(1),DTEMP(2),DTEMP(3),DTEMP(4),DTEMP(5))
          CALL DTEST(9,DTEMP,DTRUE(1,K),DTRUE(1,K),DFAC,KPRINT)
  500     CONTINUE
  600 RETURN
C                     THE FOLLOWING STOP SHOULD NEVER BE REACHED.
  999 STOP
      END
      SUBROUTINE CHECK1(SFAC,DFAC,KPRINT)
C1    ********************************* CHECK1 *************************
C     THIS SUBPROGRAM TESTS THE INCREMENTING AND ACCURACY OF THE LINEAR
C     ALGEBRA SUBPROGRAMS 26 - 38 (SNRM2 TO ICAMAX). STORED RESULTS ARE
C     COMPARED WITH THE RESULT RETURNED BY THE SUBPROGRAM.
C
C     THESE SUBPROGRAMS REQUIRE A SINGLE VECTOR ARGUMENT.
C
C     ICASE            DESIGNATES WHICH SUBPROGRAM TO TEST.
C                      26 .LE. ICASE .LE. 38
C     C. L. LAWSON, JPL, 1974 DEC 10, MAY 28
C2
      COMMON/COMBLA/NPRINT,ICASE,N,INCX,INCY,MODE,PASS
      LOGICAL          PASS
      INTEGER          ITRUE2(5),ITRUE3(5)
      DOUBLE PRECISION DA,DX(8)
      DOUBLE PRECISION  DV(8,5,2)
      DOUBLE PRECISION DFAC
      DOUBLE PRECISION DNRM2,DASUM
      DOUBLE PRECISION DTRUE1(5),DTRUE3(5),DTRUE5(8,5,2)
      REAL             STRUE2(5),STRUE4(5),STRUE(8),SX(8)
      COMPLEX          CA,CV(8,5,2),CTRUE5(8,5,2),CTRUE6(8,5,2),CX(8)
C
      DATA SA, DA, CA        / .3, .3D0, (.4,-.7)    /
      DATA DV/.1D0,2.D0,2.D0,2.D0,2.D0,2.D0,2.D0,2.D0,
     1        .3D0,3.D0,3.D0,3.D0,3.D0,3.D0,3.D0,3.D0,
     2        .3D0,-.4D0,4.D0,4.D0,4.D0,4.D0,4.D0,4.D0,
     3        .2D0,-.6D0,.3D0,5.D0,5.D0,5.D0,5.D0,5.D0,
     4        .1D0,-.3D0,.5D0,-.1D0,6.D0,6.D0,6.D0,6.D0,
     5        .1D0,8.D0,8.D0,8.D0,8.D0,8.D0,8.D0,8.D0,
     6        .3D0,9.D0,9.D0,9.D0,9.D0,9.D0,9.D0,9.D0,
     7        .3D0,2.D0,-.4D0,2.D0,2.D0,2.D0,2.D0,2.D0,
     8        .2D0,3.D0,-.6D0,5.D0,.3D0,2.D0,2.D0,2.D0,
     9         .1D0,4.D0,-.3D0,6.D0,-.5D0,7.D0,-.1D0,              3.D0/
C     COMPLEX TEST VECTORS
      DATA CV/
     1(.1,.1),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),
     2(.3,-.4),(3.,4.),(3.,4.),(3.,4.),(3.,4.),(3.,4.),(3.,4.),(3.,4.),
     3(.1,-.3),(.5,-.1),(5.,6.),(5.,6.),(5.,6.),(5.,6.),(5.,6.),(5.,6.),
     4(.1,.1),(-.6,.1),(.1,-.3),(7.,8.),(7.,8.),(7.,8.),(7.,8.),(7.,8.),
     5(.3,.1),(.1,.4),(.4,.1),(.1,.2),(2.,3.),(2.,3.),(2.,3.),(2.,3.),
     6(.1,.1),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),
     7(.3,-.4),(6.,7.),(6.,7.),(6.,7.),(6.,7.),(6.,7.),(6.,7.),(6.,7.),
     8(.1,-.3),(8.,9.),(.5,-.1),(2.,5.),(2.,5.),(2.,5.),(2.,5.),(2.,5.),
     9(.1,.1),(3.,6.),(-.6,.1),(4.,7.),(.1,-.3),(7.,2.),(7.,2.),(7.,2.),
     T(.3,.1),(5.,8.),(.1,.4),(6.,9.),(.4,.1),(8.,3.),(.1,.2),(9.,4.) /
C
      DATA STRUE2/.0,.5,.6,.7,.7/
      DATA STRUE4/.0,.7,1.,1.3,1.7/
      DATA DTRUE1/.0D0,.3D0,.5D0,.7D0,.6D0/
      DATA DTRUE3/.0D0,.3D0,.7D0,1.1D0,1.D0/
      DATA DTRUE5/.10D0,2.D0,2.D0,2.D0,2.D0,2.D0,2.D0,2.D0,
     1            .09D0,3.D0,3.D0,3.D0,3.D0,3.D0,3.D0,3.D0,
     2            .09D0,-.12D0,4.D0,4.D0,4.D0,4.D0,4.D0,4.D0,
     3            .06D0,-.18D0,.09D0,5.D0,5.D0,5.D0,5.D0,5.D0,
     4            .03D0,-.09D0,.15D0,-.03D0,6.D0,6.D0,6.D0,6.D0,
     5            .10D0,8.D0,8.D0,8.D0,8.D0,8.D0,8.D0,8.D0,
     6            .09D0,9.D0,9.D0,9.D0,9.D0,9.D0,9.D0,9.D0,
     7            .09D0,2.D0,-.12D0,2.D0,2.D0,2.D0,2.D0,2.D0,
     8            .06D0,3.D0,-.18D0,5.D0,.09D0,2.D0,2.D0,2.D0,
     9            .03D0,4.D0, -.09D0,6.D0, -.15D0,7.D0, -.03D0,  3.D0/
C
      DATA CTRUE5/
     A(.1,.1),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),
     B(-.16,-.37),(3.,4.),(3.,4.),(3.,4.),(3.,4.),(3.,4.),(3.,4.),
     C                                                         (3.,4.),
     D(-.17,-.19),(.13,-.39),(5.,6.),(5.,6.),(5.,6.),(5.,6.),(5.,6.),
     E                                                         (5.,6.),
     F(.11,-.03),(-.17,.46),(-.17,-.19),(7.,8.),(7.,8.),(7.,8.),(7.,8.),
     G                                                         (7.,8.),
     H(.19,-.17),(.32,.09),(.23,-.24),(.18,.01),(2.,3.),(2.,3.),(2.,3.),
     I                                                         (2.,3.),
     J(.1,.1),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),
     K(-.16,-.37),(6.,7.),(6.,7.),(6.,7.),(6.,7.),(6.,7.),(6.,7.),
     L                                                         (6.,7.),
     M(-.17,-.19),(8.,9.),(.13,-.39),(2.,5.),(2.,5.),(2.,5.),(2.,5.),
     N                                                         (2.,5.),
     O(.11,-.03),(3.,6.),(-.17,.46),(4.,7.),(-.17,-.19),(7.,2.),(7.,2.),
     P                                                         (7.,2.),
     Q(.19,-.17),(5.,8.),(.32,.09),(6.,9.),(.23,-.24),(8.,3.),(.18,.01),
     R                                                         (9.,4.) /
C
      DATA CTRUE6/
     A(.1,.1),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),
     B(.09,-.12),(3.,4.),(3.,4.),(3.,4.),(3.,4.),(3.,4.),(3.,4.),
     C                                                         (3.,4.),
     D(.03,-.09),(.15,-.03),(5.,6.),(5.,6.),(5.,6.),(5.,6.),(5.,6.),
     E                                                         (5.,6.),
     F(.03,.03),(-.18,.03),(.03,-.09),(7.,8.),(7.,8.),(7.,8.),(7.,8.),
     G                                                         (7.,8.),
     H(.09,.03),(.03,.12),(.12,.03),(.03,.06),(2.,3.),(2.,3.),(2.,3.),
     I                                                         (2.,3.),
     J(.1,.1),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),
     K(.09,-.12),(6.,7.),(6.,7.),(6.,7.),(6.,7.),(6.,7.),(6.,7.),
     L                                                         (6.,7.),
     M(.03,-.09),(8.,9.),(.15,-.03),(2.,5.),(2.,5.),(2.,5.),(2.,5.),
     N                                                         (2.,5.),
     O(.03,.03),(3.,6.),(-.18,.03),(4.,7.),(.03,-.09),(7.,2.),(7.,2.),
     P                                                         (7.,2.),
     Q(.09,.03),(5.,8.),(.03,.12),(6.,9.),(.12,.03),(8.,3.),(.03,.06),
     R                                                         (9.,4.) /
C
C
      DATA ITRUE2/ 0, 1, 2, 2, 3/
      DATA ITRUE3/ 0, 1, 2, 2, 2/
C
      JUMP=ICASE-25
         DO 520 INCX=1,2
            DO 500 NP1=1,5
            N=NP1-1
            LEN= 2*MAX0(N,1)
C                                                  SET VECTOR ARGUMENTS.
                    DO 22 I = 1, LEN
                    SX(I) = SNGL(DV(I,NP1,INCX))
                    DX(I) = DV(I,NP1,INCX)
   22               CX(I) = CV(I,NP1,INCX)
C
C                        BRANCH TO INVOKE SUBPROGRAM TO BE TESTED.
C
               GO TO (260,270,280,290,300,310,320,
     *                330,340,350,360,370,380),JUMP
C                                                             26. SNRM2
  260       STEMP=SNGL(DTRUE1(NP1))
            CALL STEST(1,SNRM2(N,SX,INCX),STEMP,STEMP,SFAC,KPRINT)
            GO TO 500
C                                                             27. DNRM2
  270       CALL DTEST(1,DNRM2(N,DX,INCX),DTRUE1(NP1),DTRUE1(NP1),DFAC,
     1                 KPRINT)
            GO TO 500
C                                                             28. SCNRM2
  280       CALL STEST(1,SCNRM2(N,CX,INCX),STRUE2(NP1),STRUE2(NP1),
     1                 SFAC,KPRINT)
            GO TO 500
C                                                             29. SASUM
  290       STEMP=SNGL(DTRUE3(NP1))
            CALL STEST(1,SASUM(N,SX,INCX),STEMP,STEMP,SFAC,KPRINT)
            GO TO 500
C                                                             30. DASUM
  300       CALL DTEST(1,DASUM(N,DX,INCX),DTRUE3(NP1),DTRUE3(NP1),DFAC,
     1                 KPRINT)
            GO TO 500
C                                                             31. SCASUM
  310       CALL STEST(1,SCASUM(N,CX,INCX),STRUE4(NP1),STRUE4(NP1),SFAC,
     1                 KPRINT)
            GO TO 500
C                                                             32. SSCALE
  320       CALL SSCAL(N,SA,SX,INCX)
               DO 322 I = 1, LEN
  322          STRUE(I) = SNGL(DTRUE5(I,NP1,INCX))
            CALL STEST(LEN,SX,STRUE,STRUE,SFAC,KPRINT)
            GO TO 500
C                                                             33. DSCALE
  330       CALL DSCAL(N,DA,DX,INCX)
           CALL DTEST(LEN,DX,DTRUE5(1,NP1,INCX),DTRUE5(1,NP1,INCX),
     1                 DFAC,KPRINT)
            GO TO 500
C                                                             34. CSCALE
  340       CALL CSCAL(N,CA,CX,INCX)
        CALL STEST(2*LEN,CX,CTRUE5(1,NP1,INCX),CTRUE5(1,NP1,INCX),
     1                 SFAC,KPRINT)
            GO TO 500
C                                                             35. CSSCAL
  350       CALL CSSCAL(N,SA,CX,INCX)
         CALL STEST(2*LEN,CX,CTRUE6(1,NP1,INCX),CTRUE6(1,NP1,INCX),
     1                 SFAC,KPRINT)
            GO TO 500
C                                                             36. ISAMAX
  360       CALL ITEST(1,ISAMAX(N,SX,INCX),ITRUE2(NP1),KPRINT)
            GO TO 500
C                                                             37. IDAMAX
  370       CALL ITEST(1,IDAMAX(N,DX,INCX),ITRUE2(NP1),KPRINT)
            GO TO 500
C                                                             38. ICAMAX
  380       CALL ITEST(1,ICAMAX(N,CX,INCX),ITRUE3(NP1),KPRINT)
C
  500       CONTINUE
  520    CONTINUE
      RETURN
      END
      SUBROUTINE CHECK2(SFAC,SDFAC,DFAC,DQFAC,KPRINT)
C1    ********************************* CHECK2 *************************
C     THIS SUBPROGRAM TESTS THE BASIC LINEAR ALGEBRA SUBPROGRAMS 1-11,
C     14-15, AND 18-25. SUBPROGRAMS IN THIS SET EACH REQUIRE TWO ARRAYS
C     IN THE PARAMETER LIST.
C
C     C. L. LAWSON, JPL, 1975 FEB 26, APR 29, MAY 8, MAY 28
C2
      COMMON/COMBLA/NPRINT,ICASE,N,INCX,INCY,MODE,PASS
C
      LOGICAL          PASS
      INTEGER          INCXS(4),INCYS(4),LENS(4,2),NS(4)
      REAL             SX(7),SY(7),STX(7),STY(7),SSIZE1(4),SSIZE2(14,2)
      REAL             SSIZE(7),QC(10),SPARAM(5),ST7B(4,4),SSIZE3(4)
      DOUBLE PRECISION DX(7),DA,DX1(7),DY1(7),DY(7),DT7(4,4),DT8(7,4,4)
      DOUBLE PRECISION DX2(7), DY2(7), DT2(4,4,2), DPARAM(5), DPAR(5,4)
      DOUBLE PRECISION DSDOT,DDOT,DQDOTI,DQDOTA,DFAC,DQFAC
      DOUBLE PRECISION DT10X(7,4,4),DT10Y(7,4,4),DB
      DOUBLE PRECISION DSIZE1(4),DSIZE2(7,2),DSIZE(7)
      DOUBLE PRECISION DC,DS,DT9X(7,4,4),DT9Y(7,4,4),DTX(7),DTY(7)
      DOUBLE PRECISION DT19X(7,4,16),DT19XA(7,4,4),DT19XB(7,4,4)
      DOUBLE PRECISION DT19XC(7,4,4),DT19XD(7,4,4),DT19Y(7,4,16)
      DOUBLE PRECISION DT19YA(7,4,4),DT19YB(7,4,4),DT19YC(7,4,4)
      DOUBLE PRECISION DT19YD(7,4,4)
C
      EQUIVALENCE (DT19X(1,1,1),DT19XA(1,1,1)),(DT19X(1,1,5),
     A   DT19XB(1,1,1)),(DT19X(1,1,9),DT19XC(1,1,1)),
     B   (DT19X(1,1,13),DT19XD(1,1,1))
      EQUIVALENCE (DT19Y(1,1,1),DT19YA(1,1,1)),(DT19Y(1,1,5),
     A   DT19YB(1,1,1)),(DT19Y(1,1,9),DT19YC(1,1,1)),
     B   (DT19Y(1,1,13),DT19YD(1,1,1))
      COMPLEX          CX(7),CA,CX1(7),CY1(7),CY(7),CT6(4,4),CT7(4,4)
      COMPLEX          CT8(7,4,4),CSIZE1(4),CSIZE2(7,2)
      COMPLEX          CT10X(7,4,4), CT10Y(7,4,4)
      COMPLEX          CDOTC,CDOTU
      DATA SA,DA,CA,DB,SB/.3,.3D0,(.4,-.7),.25D0,.1/
      DATA INCXS/   1,   2,  -2,  -1 /
      DATA INCYS/   1,  -2,   1,  -2 /
      DATA LENS/1, 1, 2, 4,   1, 1, 3, 7/
      DATA NS   /   0,   1,   2,   4 /
      DATA SC,SS,DC,DS/ .8,.6,.8D0,.6D0/
      DATA DX1/ .6D0, .1D0,-.5D0, .8D0, .9D0,-.3D0,-.4D0/
      DATA DY1/ .5D0,-.9D0, .3D0, .7D0,-.6D0, .2D0, .8D0/
      DATA DX2/ 1.D0,.01D0, .02D0,1.D0,.06D0, 2.D0, 1.D0/
      DATA DY2/ 1.D0,.04D0,-.03D0,-1.D0,.05D0,3.D0,-1.D0/
C            THE TERMS D11(3,2) AND D11(4,2) WILL BE SET BY
C            COMPUTATION AT RUN TIME.
      DATA CX1/(.7,-.8),(-.4,-.7),(-.1,-.9),(.2,-.8),(-.9,-.4),(.1,.4),
     *                                                        (-.6,.6)/
      DATA CY1/(.6,-.6),(-.9,.5),(.7,-.6),(.1,-.5),(-.1,-.2),(-.5,-.3),
     *                                                       (.8,-.7) /
C
C                             FOR DQDOTI AND DQDOTA
C
      DATA DT2/0.25D0,1.25D0,1.2504D0,0.2498D0,
     A         0.25D0,1.25D0,0.24D0,0.2492D0,
     B         0.25D0,1.25D0,0.31D0,0.2518D0,
     C         0.25D0,1.25D0,1.2497D0,0.2507D0,
     D         0.D0,2.D0,2.0008D0,-.0004D0,
     E         0.D0,2.D0,-.02D0,-.0016D0,
     F         0.D0,2.D0,.12D0,.0036D0,
     G         0.D0,2.D0,1.9994D0,.0014D0/
      DATA DT7/ 0.D0,.30D0,.21D0,.62D0,      0.D0,.30D0,-.07D0,.85D0,
     *          0.D0,.30D0,-.79D0,-.74D0,    0.D0,.30D0,.33D0,1.27D0/
      DATA ST7B/ .1, .4, .31, .72,     .1, .4, .03, .95,
     *           .1, .4, -.69, -.64,   .1, .4, .43, 1.37/
C
C                       FOR CDOTU
C
      DATA CT7/(0.,0.),(-.06,-.90),(.65,-.47),(-.34,-1.22),
     1         (0.,0.),(-.06,-.90),(-.59,-1.46),(-1.04,-.04),
     2         (0.,0.),(-.06,-.90),(-.83,.59),  (  .07,-.37),
     3         (0.,0.),(-.06,-.90),(-.76,-1.15),(-1.33,-1.82)/
C
C                       FOR CDOTC
C
      DATA CT6/(0.,0.),(.90,0.06), (.91,-.77),    (1.80,-.10),
     A         (0.,0.),(.90,0.06), (1.45,.74),    (.20,.90),
     B         (0.,0.),(.90,0.06), (-.55,.23),    (.83,-.39),
     C         (0.,0.),(.90,0.06), (1.04,0.79),    (1.95,1.22)/
C
      DATA DT8/.5D0,                     0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     1         .68D0,                    0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     2         .68D0,-.87D0,                 0.D0,0.D0,0.D0,0.D0,0.D0,
     3         .68D0,-.87D0,.15D0,.94D0,          0.D0,0.D0,0.D0,
     4         .5D0,                     0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     5         .68D0,                    0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     6         .35D0,-.9D0,.48D0,                   0.D0,0.D0,0.D0,0.D0,
     7         .38D0,-.9D0,.57D0,.7D0,-.75D0,.2D0,.98D0,
     8         .5D0,                      0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     9         .68D0,                     0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A         .35D0,-.72D0,                0.D0,0.D0,0.D0,0.D0,0.D0,
     B         .38D0,-.63D0,.15D0,.88D0,                 0.D0,0.D0,0.D0,
     C         .5D0,                      0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D         .68D0,                     0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E         .68D0,-.9D0,.33D0,                0.D0,0.D0,0.D0,0.D0,
     F         .68D0,-.9D0,.33D0,.7D0,-.75D0,.2D0,1.04D0/
C
      DATA CT8/
     A(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     B(.32,-1.41),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     C(.32,-1.41),(-1.55,.5),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     D(.32,-1.41),(-1.55,.5),(.03,-.89),(-.38,-.96),(0.,0.),(0.,0.),
     E                                                         (0.,0.),
     F(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     G(.32,-1.41),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     H(-.07,-.89),(-.9,.5),(.42,-1.41),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     I(.78,.06),(-.9,.5),(.06,-.13),(.1,-.5),(-.77,-.49),(-.5,-.3),
     J                                                     (.52,-1.51),
     K(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     L(.32,-1.41),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     M(-.07,-.89),(-1.18,-.31),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     N(.78,.06),(-1.54,.97),(.03,-.89),(-.18,-1.31),(0.,0.),(0.,0.),
     O(0.,0.),(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     P(.32,-1.41),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     Q(.32,-1.41),(-.9,.5),(.05,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     R(.32,-1.41),(-.9,.5),(.05,-.6),(.1,-.5),(-.77,-.49),(-.5,-.3),
     S                                                     (.32,-1.16) /
C
C
C                TRUE X VALUES AFTER ROTATION USING SROT OR DROT.
      DATA DT9X/.6D0,                    0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A          .78D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B          .78D0,-.46D0,               0.D0,0.D0,0.D0,0.D0,0.D0,
     C          .78D0,-.46D0,-.22D0,1.06D0,              0.D0,0.D0,0.D0,
     D          .6D0,                    0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E          .78D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F          .66D0,.1D0,-.1D0,                   0.D0,0.D0,0.D0,0.D0,
     G          .96D0,.1D0,-.76D0,.8D0,.90D0,-.3D0,-.02D0,
     H          .6D0,                    0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     I          .78D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     J          -.06D0,.1D0,-.1D0,                  0.D0,0.D0,0.D0,0.D0,
     K          .90D0,.1D0,-.22D0,.8D0,.18D0,-.3D0,-.02D0,
     L          .6D0,                    0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     M          .78D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     N          .78D0,.26D0,                0.D0,0.D0,0.D0,0.D0,0.D0,
     O          .78D0,.26D0,-.76D0,1.12D0,               0.D0,0.D0,0.D0/
C
C                TRUE Y VALUES AFTER ROTATION USING SROT OR DROT.
C
      DATA DT9Y/ .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A           .04D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B           .04D0,-.78D0,              0.D0,0.D0,0.D0,0.D0,0.D0,
     C           .04D0,-.78D0, .54D0, .08D0,             0.D0,0.D0,0.D0,
     D           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E           .04D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           .7D0,-.9D0,-.12D0,                 0.D0,0.D0,0.D0,0.D0,
     G           .64D0,-.9D0,-.30D0, .7D0,-.18D0, .2D0, .28D0,
     H           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     I           .04D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     J           .7D0,-1.08D0,              0.D0,0.D0,0.D0,0.D0,0.D0,
     K           .64D0,-1.26D0,.54D0, .20D0,             0.D0,0.D0,0.D0,
     L           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     M          .04D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     N           .04D0,-.9D0, .18D0,                0.D0,0.D0,0.D0,0.D0,
     O           .04D0,-.9D0, .18D0, .7D0,-.18D0, .2D0, .16D0/
C
      DATA DT10X/.6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B           .5D0,-.9D0,                0.D0,0.D0,0.D0,0.D0,0.D0,
     C           .5D0,-.9D0,.3D0,.7D0,                   0.D0,0.D0,0.D0,
     D           .6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           .3D0,.1D0 ,.5D0,                   0.D0,0.D0,0.D0,0.D0,
     G           .8D0,.1D0 ,-.6D0,.8D0 ,.3D0,-.3D0,.5D0,
     H           .6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     I           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     J           -.9D0,.1D0,.5D0,                   0.D0,0.D0,0.D0,0.D0,
     K           .7D0, .1D0,.3D0, .8D0,-.9D0,-.3D0,.5D0,
     L           .6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     M           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     N           .5D0,.3D0,                 0.D0,0.D0,0.D0,0.D0,0.D0,
     O           .5D0,.3D0,-.6D0,.8D0,                   0.D0,0.D0,0.D0/
C
      DATA DT10Y/.5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A           .6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B           .6D0,.1D0,                 0.D0,0.D0,0.D0,0.D0,0.D0,
     C           .6D0,.1D0,-.5D0,.8D0,                   0.D0,0.D0,0.D0,
     D           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E           .6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           -.5D0,-.9D0,.6D0,                  0.D0,0.D0,0.D0,0.D0,
     G           -.4D0,-.9D0,.9D0, .7D0,-.5D0, .2D0,.6D0,
     H           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     I           .6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     J           -.5D0,.6D0,                0.D0,0.D0,0.D0,0.D0,0.D0,
     K           -.4D0,.9D0,-.5D0,.6D0,                  0.D0,0.D0,0.D0,
     L           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     M           .6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     N           .6D0,-.9D0,.1D0,                   0.D0,0.D0,0.D0,0.D0,
     O           .6D0,-.9D0,.1D0, .7D0,-.5D0, .2D0, .8D0/
C
      DATA CT10X/
     A(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     B(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     C(.6,-.6),(-.9,.5),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     D(.6,-.6),(-.9,.5),(.7,-.6),(.1,-.5),(0.,0.),(0.,0.),(0.,0.),
     E(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     F(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     G(.7,-.6),(-.4,-.7),(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     H(.8,-.7),(-.4,-.7),(-.1,-.2),(.2,-.8),(.7,-.6),(.1,.4),(.6,-.6),
     I(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     J(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     K(-.9,.5),(-.4,-.7),(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     L(.1,-.5),(-.4,-.7),(.7,-.6),(.2,-.8),(-.9,.5),(.1,.4),(.6,-.6),
     M(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     N(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     O(.6,-.6),(.7,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     P(.6,-.6),(.7,-.6),(-.1,-.2),(.8,-.7),(0.,0.),(0.,0.),(0.,0.)   /
C
      DATA CT10Y/
     A(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     B(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     C(.7,-.8),(-.4,-.7),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     D(.7,-.8),(-.4,-.7),(-.1,-.9),(.2,-.8),(0.,0.),(0.,0.),(0.,0.),
     E(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     F(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     G(-.1,-.9),(-.9,.5),(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     H(-.6,.6),(-.9,.5),(-.9,-.4),(.1,-.5),(-.1,-.9),(-.5,-.3),(.7,-.8),
     I(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     J(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     K(-.1,-.9),(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     L(-.6,.6),(-.9,-.4),(-.1,-.9),(.7,-.8),(0.,0.),(0.,0.),(0.,0.),
     M(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     N(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     O(.7,-.8),(-.9,.5),(-.4,-.7),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     P(.7,-.8),(-.9,.5),(-.4,-.7),(.1,-.5),(-.1,-.9),(-.5,-.3),(.2,-.8)/
C                        TRUE X RESULTS F0R ROTATIONS SROTM AND DROTM
      DATA DT19XA/.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E           -.8D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           -.9D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G           3.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .6D0,   .1D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     I           -.8D0,  3.8D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     J           -.9D0,  2.8D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     K           3.5D0,  -.4D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     L            .6D0,   .1D0,  -.5D0,   .8D0,          0.D0,0.D0,0.D0,
     M           -.8D0,  3.8D0, -2.2D0, -1.2D0,          0.D0,0.D0,0.D0,
     N           -.9D0,  2.8D0, -1.4D0, -1.3D0,          0.D0,0.D0,0.D0,
     O           3.5D0,  -.4D0, -2.2D0,  4.7D0,          0.D0,0.D0,0.D0/
C
      DATA DT19XB/.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E           -.8D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           -.9D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G           3.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .6D0,   .1D0,  -.5D0,             0.D0,0.D0,0.D0,0.D0,
     I           0.D0,    .1D0, -3.0D0,             0.D0,0.D0,0.D0,0.D0,
     J           -.3D0,   .1D0, -2.0D0,             0.D0,0.D0,0.D0,0.D0,
     K           3.3D0,   .1D0, -2.0D0,             0.D0,0.D0,0.D0,0.D0,
     L            .6D0,   .1D0,  -.5D0,   .8D0,   .9D0,  -.3D0,  -.4D0,
     M          -2.0D0,   .1D0,  1.4D0,   .8D0,   .6D0,  -.3D0, -2.8D0,
     N          -1.8D0,   .1D0,  1.3D0,   .8D0,  0.D0,   -.3D0, -1.9D0,
     O           3.8D0,   .1D0, -3.1D0,   .8D0,  4.8D0,  -.3D0, -1.5D0 /
C
      DATA DT19XC/.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E           -.8D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           -.9D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G           3.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .6D0,   .1D0,  -.5D0,             0.D0,0.D0,0.D0,0.D0,
     I           4.8D0,   .1D0, -3.0D0,             0.D0,0.D0,0.D0,0.D0,
     J           3.3D0,   .1D0, -2.0D0,             0.D0,0.D0,0.D0,0.D0,
     K           2.1D0,   .1D0, -2.0D0,             0.D0,0.D0,0.D0,0.D0,
     L            .6D0,   .1D0,  -.5D0,   .8D0,   .9D0,  -.3D0,  -.4D0,
     M          -1.6D0,   .1D0, -2.2D0,   .8D0,  5.4D0,  -.3D0, -2.8D0,
     N          -1.5D0,   .1D0, -1.4D0,   .8D0,  3.6D0,  -.3D0, -1.9D0,
     O           3.7D0,   .1D0, -2.2D0,   .8D0,  3.6D0,  -.3D0, -1.5D0 /
C
      DATA DT19XD/.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E           -.8D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           -.9D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G           3.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .6D0,   .1D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     I           -.8D0, -1.0D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     J           -.9D0,  -.8D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     K           3.5D0,   .8D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     L            .6D0,   .1D0,  -.5D0,   .8D0,          0.D0,0.D0,0.D0,
     M           -.8D0, -1.0D0,  1.4D0, -1.6D0,          0.D0,0.D0,0.D0,
     N           -.9D0,  -.8D0,  1.3D0, -1.6D0,          0.D0,0.D0,0.D0,
     O           3.5D0,   .8D0, -3.1D0,  4.8D0,          0.D0,0.D0,0.D0/
C                        TRUE Y RESULTS FOR ROTATIONS SROTM AND DROTM
      DATA DT19YA/.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E            .7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           1.7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G          -2.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .5D0,  -.9D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     I            .7D0, -4.8D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     J           1.7D0,  -.7D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     K          -2.6D0,  3.5D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     L            .5D0,  -.9D0,   .3D0,   .7D0,          0.D0,0.D0,0.D0,
     M            .7D0, -4.8D0,  3.0D0,  1.1D0,          0.D0,0.D0,0.D0,
     N           1.7D0,  -.7D0,  -.7D0,  2.3D0,          0.D0,0.D0,0.D0,
     O          -2.6D0,  3.5D0,  -.7D0, -3.6D0,          0.D0,0.D0,0.D0/
C
      DATA DT19YB/.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E            .7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           1.7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G          -2.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .5D0,  -.9D0,   .3D0,             0.D0,0.D0,0.D0,0.D0,
     I           4.0D0,  -.9D0,  -.3D0,             0.D0,0.D0,0.D0,0.D0,
     J           -.5D0,  -.9D0,  1.5D0,             0.D0,0.D0,0.D0,0.D0,
     K          -1.5D0,  -.9D0, -1.8D0,             0.D0,0.D0,0.D0,0.D0,
     L            .5D0,  -.9D0,   .3D0,   .7D0,  -.6D0,   .2D0,   .8D0,
     M           3.7D0,  -.9D0, -1.2D0,   .7D0, -1.5D0,   .2D0,  2.2D0,
     N           -.3D0,  -.9D0,  2.1D0,   .7D0, -1.6D0,   .2D0,  2.0D0,
     O          -1.6D0,  -.9D0, -2.1D0,   .7D0,  2.9D0,   .2D0, -3.8D0 /
C
      DATA DT19YC/.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E            .7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           1.7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G          -2.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .5D0,  -.9D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     I           4.0D0, -6.3D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     J           -.5D0,   .3D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     K          -1.5D0,  3.0D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     L            .5D0,  -.9D0,   .3D0,   .7D0,          0.D0,0.D0,0.D0,
     M           3.7D0, -7.2D0,  3.0D0,  1.7D0,          0.D0,0.D0,0.D0,
     N           -.3D0,   .9D0,  -.7D0,  1.9D0,          0.D0,0.D0,0.D0,
     O          -1.6D0,  2.7D0,  -.7D0, -3.4D0,          0.D0,0.D0,0.D0/
C
      DATA DT19YD/.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E            .7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           1.7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G          -2.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .5D0,  -.9D0,   .3D0,             0.D0,0.D0,0.D0,0.D0,
     I            .7D0,  -.9D0,  1.2D0,             0.D0,0.D0,0.D0,0.D0,
     J           1.7D0,  -.9D0,   .5D0,             0.D0,0.D0,0.D0,0.D0,
     K          -2.6D0,  -.9D0, -1.3D0,             0.D0,0.D0,0.D0,0.D0,
     L            .5D0,  -.9D0,   .3D0,   .7D0,  -.6D0,   .2D0,   .8D0,
     M            .7D0,  -.9D0,  1.2D0,   .7D0, -1.5D0,   .2D0,  1.6D0,
     N           1.7D0,  -.9D0,   .5D0,   .7D0, -1.6D0,   .2D0,  2.4D0,
     O          -2.6D0,  -.9D0, -1.3D0,   .7D0,  2.9D0,   .2D0, -4.0D0 /
C
      DATA SSIZE1/ 0.  , .3  , 1.6  , 3.2   /
      DATA DSIZE1/ 0.D0, .3D0, 1.6D0, 3.2D0 /
      DATA SSIZE3/ .1, .4, 1.7, 3.3 /
C
C                         FOR CDOTC AND CDOTU
C
      DATA CSIZE1/ (0.,0.), (.9,.9), (1.63,1.73), (2.90,2.78) /
      DATA SSIZE2/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     A  1.17,1.17,1.17,1.17,1.17,1.17,1.17,
     B  1.17,1.17,1.17,1.17,1.17,1.17,1.17/
      DATA DSIZE2/0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A  1.17D0,1.17D0,1.17D0,1.17D0,1.17D0,1.17D0,1.17D0/
C
C                         FOR CAXPY
C
      DATA CSIZE2/
     A (0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     B (1.54,1.54),(1.54,1.54),(1.54,1.54),(1.54,1.54),(1.54,1.54),
     C                                     (1.54,1.54),(1.54,1.54) /
C
C                         FOR SROTM AND DROTM
C
      DATA DPAR/-2.D0,  0.D0,0.D0,0.D0,0.D0,
     A          -1.D0,  2.D0, -3.D0, -4.D0,  5.D0,
     B           0.D0,  0.D0,  2.D0, -3.D0,  0.D0,
     C           1.D0,  5.D0,  2.D0,  0.D0, -4.D0/
C
        DO 520 KI = 1, 4
        INCX = INCXS(KI)
        INCY = INCYS(KI)
        MX   = IABS(INCX)
        MY   = IABS(INCY)
C
          DO 500 KN=1,4
          N= NS(KN)
          KSIZE=MIN0(2,KN)
          LENX = LENS(KN,MX)
          LENY = LENS(KN,MY)
C                                       INITIALIZE ALL ARGUMENT ARRAYS.
               DO 5 I = 1, 7
               SX(I) = SNGL(DX1(I))
               SY(I) = SNGL(DY1(I))
               DX(I) = DX1(I)
               DY(I) = DY1(I)
               CX(I) = CX1(I)
    5          CY(I) = CY1(I)
C
C                             BRANCH TO SELECT SUBPROGRAM TO BE TESTED.
C
          GO TO ( 10, 20, 30, 40, 50, 60, 70, 80, 90,100,
     A           110,999,999,140,150,999,999,180,190,200,
     B           210,220,230,240,250), ICASE
C                                                              1. SDOT
   10     CALL STEST(1,SDOT(N,SX,INCX,SY,INCY),SNGL(DT7(KN,KI)),
     *                                         SSIZE1(KN),SFAC,KPRINT)
          GO TO 500
C                                                              2. DSDOT
   20     CALL STEST(1,SNGL(DSDOT(N,SX,INCX,SY,INCY)),
     *               SNGL(DT7(KN,KI)),SSIZE1(KN),SFAC,KPRINT)
          GO TO 500
C                                                              3. SDSDOT
   30     CALL STEST(1,SDSDOT(N,SB,SX,INCX,SY,INCY),
     *               ST7B(KN,KI),SSIZE3(KN),SFAC,KPRINT)
          GO TO 500
C                                                              4. DDOT
   40     CALL DTEST(1,DDOT(N,DX,INCX,DY,INCY),DT7(KN,KI),
     *               DSIZE1(KN),DFAC,KPRINT)
          GO TO 500
C                                                              5. DQDOTI
   50 CONTINUE
C                        DQDOTI AND DQDOTA ARE SUPPOSED TO USE EXTENDED
C                        PRECISION ARITHMETIC INTERNALLY.
C     SET MODE = 1 OR 2 TO DISTINGUISH TESTS OF DQDOTI OR DQDOTA
C     IN THE DIAGNOSTIC OUTPUT.
C
C         MODE = 1
C         CALL DTEST(1,DQDOTI(N,DB,QC,DX2,INCX,DY2,INCY),
C    *               DT2(KN,KI,1),DT2(KN,KI,1),DQFAC,KPRINT)
C     GO TO 500
C                                                              6. DQDOTA
   60 CONTINUE
C     TO TEST DQDOTA WE ACTUALLY TEST BOTH DQDOTI AND DQDOTA.
C     THE OUTPUT VALUE OF QX FROM DQDOTI WILL BE USED AS INPUT
C     TO DQDOTA.  QX IS SUPPOSED TO BE IN A MACHINE-DEPENDENT
C     EXTENDED PRECISION FORM.
C     MODE IS SET TO 1 OR 2 TO DISTINGUISH TESTS OF
C     DQDOTI OR DQDOTA IN THE DIAGNOSTIC OUTPUT.
C
C         MODE = 1
C         CALL DTEST(1,DQDOTI(N,DB,QC,DX2,INCX,DY2,INCY),
C    *               DT2(KN,KI,1),DT2(KN,KI,1),DQFAC,KPRINT)
C         MODE = 2
C         CALL DTEST(1,DQDOTA(N,-DB,QC,DX2,INCX,DY2,INCY),
C    *               DT2(KN,KI,2),DT2(KN,KI,2),DQFAC,KPRINT)
C         GO TO 500
C                                                              7. CDOTC
   70     CALL STEST(2, CDOTC(N,CX,INCX,CY,INCY),
     *               CT6(KN,KI),CSIZE1(KN),SFAC,KPRINT)
          GO TO 500
C                                                              8. CDOTU
   80     CALL STEST(2,CDOTU(N,CX,INCX,CY,INCY),
     *               CT7(KN,KI),CSIZE1(KN),SFAC,KPRINT)
          GO TO 500
C                                                              9. SAXPY
   90     CALL SAXPY(N,SA,SX,INCX,SY,INCY)
               DO 95 J = 1, LENY
   95          STY(J) = SNGL(DT8(J,KN,KI))
          CALL STEST(LENY,SY,STY,SSIZE2(1,KSIZE),SFAC,KPRINT)
          GO TO 500
C                                                              10. DAXPY
  100      CALL DAXPY(N,DA,DX,INCX,DY,INCY)
          CALL DTEST(LENY,DY,DT8(1,KN,KI),DSIZE2(1,KSIZE),DFAC,KPRINT)
          GO TO 500
C                                                              11. CAXPY
  110     CALL CAXPY(N,CA,CX,INCX,CY,INCY)
          CALL STEST(2*LENY,CY,CT8(1,KN,KI),CSIZE2(1,KSIZE),SFAC,KPRINT)
          GO TO 500
C                                                              14. SROT
  140     CONTINUE
               DO 144 I = 1, 7
               SX(I) = SNGL(DX1(I))
               SY(I) = SNGL(DY1(I))
               STX(I) = SNGL(DT9X(I,KN,KI))
               STY(I) = SNGL(DT9Y(I,KN,KI))
  144         CONTINUE
          CALL SROT   (N,SX,INCX,SY,INCY,SC,SS)
          CALL STEST(LENX,SX,STX,SSIZE2(1,KSIZE),SFAC,KPRINT)
          CALL STEST(LENY,SY,STY,SSIZE2(1,KSIZE),SFAC,KPRINT)
          GO TO 500
C                                                             15. DROT
  150     CONTINUE
               DO 154 I = 1, 7
               DX(I) = DX1(I)
               DY(I) = DY1(I)
  154          CONTINUE
          CALL DROT   (N,DX,INCX,DY,INCY,DC,DS)
          CALL DTEST(LENX,DX,DT9X(1,KN,KI),DSIZE2(1,KSIZE),DFAC,KPRINT)
          CALL DTEST(LENY,DY,DT9Y(1,KN,KI),DSIZE2(1,KSIZE),DFAC,KPRINT)
          GO TO 500
C                                                             18. SROTM
  180     KNI = KN + 4*(KI-1)
          DO 189 KPAR=1,4
          DO 182 I = 1, 7
          SX(I) = SNGL(DX1(I))
          SY(I) = SNGL(DY1(I))
          STX(I) = SNGL(DT19X(I,KPAR,KNI))
  182     STY(I) = SNGL(DT19Y(I,KPAR,KNI))
C
          DO 186 I = 1, 5
  186     SPARAM(I) = SNGL(DPAR(I,KPAR))
C                          SET MODE TO IDENTIFY DIAGNOSTIC OUTPUT,
C                          IF ANY
          MODE = INT(SPARAM(1))
C
          DO 187 I = 1, LENX
  187     SSIZE(I) = STX(I)
C                         THE TRUE RESULTS DT19X(1,2,7) AND
C                         DT19X(5,3,8) ARE ZERO DUE TO CANCELLATION.
C                         DT19X(1,2,7) = 2.*.6 - 4.*.3 = 0
C                         DT19X(5,3,8) = .9 - 3.*.3 = 0
C                         FOR THESE CASES RESPECTIVELY SET SIZE( )
C                         EQUAL TO 2.4 AND 1.8
          IF ((KPAR .EQ. 2) .AND. (KNI .EQ. 7))
     1           SSIZE(1) = 2.4E0
          IF ((KPAR .EQ. 3) .AND. (KNI .EQ. 8))
     1           SSIZE(5) = 1.8E0
C
          CALL SROTM(N,SX,INCX,SY,INCY,SPARAM)
          CALL STEST(LENX,SX,STX,SSIZE,SFAC,KPRINT)
          CALL STEST(LENY,SY,STY,STY,SFAC,KPRINT)
  189     CONTINUE
          GO TO 500
C                                                             19. DROTM
  190     KNI = KN + 4*(KI-1)
          DO 199 KPAR=1,4
            DO 192 I = 1, 7
            DX(I) = DX1(I)
            DY(I) = DY1(I)
            DTX(I) = DT19X(I,KPAR,KNI)
  192       DTY(I) = DT19Y(I,KPAR,KNI)
C
            DO 196 I = 1, 5
  196       DPARAM(I) = DPAR(I,KPAR)
C                            SET MODE TO IDENTIFY DIAGNOSTIC OUTPUT,
C                            IF ANY
*-----------------------------------------------------------------------
*  Changed by RFB 2-Apr-98 ... f77 compiler on SGI Origin fails on this
*          MODE = IDINT(DPARAM(1))
          MODE = INT(DPARAM(1))
*-----------------------------------------------------------------------
C
            DO 197 I = 1, LENX
  197       DSIZE(I) = DTX(I)
C                             SEE REMARK ABOVE ABOUT DT11X(1,2,7)
C                             AND DT11X(5,3,8).
          IF ((KPAR .EQ. 2) .AND. (KNI .EQ. 7))
     1               DSIZE(1) = 2.4D0
          IF ((KPAR .EQ. 3) .AND. (KNI .EQ. 8))
     1               DSIZE(5) = 1.8D0
C
          CALL   DROTM(N,DX,INCX,DY,INCY,DPARAM)
          CALL DTEST(LENX,DX,DTX,DSIZE,DFAC,KPRINT)
          CALL DTEST(LENY,DY,DTY,DTY,DFAC,KPRINT)
  199     CONTINUE
          GO TO 500
C                                                             20. SCOPY
  200     DO 205 I = 1, 7
  205     STY(I) = SNGL(DT10Y( I,KN,KI))
          CALL SCOPY(N,SX,INCX,SY,INCY)
          CALL STEST(LENY,SY,STY,SSIZE2(1,1),1.,KPRINT)
          GO TO 500
C                                                             21. DCOPY
  210     CALL DCOPY(N,DX,INCX,DY,INCY)
          CALL DTEST(LENY,DY,DT10Y(1,KN,KI),DSIZE2(1,1),1.D0,KPRINT)
          GO TO 500
C                                                             22. CCOPY
  220     CALL CCOPY(N,CX,INCX,CY,INCY)
          CALL STEST(2*LENY,CY,CT10Y(1,KN,KI),SSIZE2(1,1),1.,KPRINT)
          GO TO 500
C                                                             23. SSWAP
  230     CALL SSWAP(N,SX,INCX,SY,INCY)
               DO 235 I = 1, 7
               STX(I) = SNGL(DT10X(I,KN,KI))
  235          STY(I) = SNGL(DT10Y(I,KN,KI))
          CALL STEST(LENX,SX,STX,SSIZE2(1,1),1.,KPRINT)
          CALL STEST(LENY,SY,STY,SSIZE2(1,1),1.,KPRINT)
          GO TO 500
C                                                             24. DSWAP
  240     CALL DSWAP(N,DX,INCX,DY,INCY)
          CALL DTEST(LENX,DX,DT10X(1,KN,KI),DSIZE2(1,1),1.D0,KPRINT)
          CALL DTEST(LENY,DY,DT10Y(1,KN,KI),DSIZE2(1,1),1.D0,KPRINT)
          GO TO 500
C                                                             25. CSWAP
  250     CALL CSWAP(N,CX,INCX,CY,INCY)
          CALL STEST(2*LENX,CX,CT10X(1,KN,KI),SSIZE2(1,1),1.,KPRINT)
          CALL STEST(2*LENY,CY,CT10Y(1,KN,KI),SSIZE2(1,1),1.,KPRINT)
C
C
C
  500     CONTINUE
  520   CONTINUE
      RETURN
C                 THE FOLLOWING STOP SHOULD NEVER BE REACHED.
  999 STOP
      END
      SUBROUTINE HEADER (KPRINT)
C1    ********************************* HEADER *************************
C     PRINT HEADER FOR CASE
C     C. L. LAWSON, JPL, 1974 DEC 12
C2
      COMMON/COMBLA/NPRINT,ICASE,N,INCX,INCY,MODE,PASS
      LOGICAL          PASS
      DIMENSION        L(3,38)
C
      DATA L(1, 1),L(2, 1),L(3, 1)/2H  ,2HSD,2HOT/
      DATA L(1, 2),L(2, 2),L(3, 2)/2H D,2HSD,2HOT/
      DATA L(1, 3),L(2, 3),L(3, 3)/2HSD,2HSD,2HOT/
      DATA L(1, 4),L(2, 4),L(3, 4)/2H  ,2HDD,2HOT/
      DATA L(1, 5),L(2, 5),L(3, 5)/2HDQ,2HDO,2HTI/
      DATA L(1, 6),L(2, 6),L(3, 6)/2HDQ,2HDO,2HTA/
      DATA L(1,7),L(2,7),L(3,7)/2H C,2HDO,2HTC/
      DATA L(1, 8),L(2, 8),L(3, 8)/2H C,2HDO,2HTU/
      DATA L(1, 9),L(2, 9),L(3, 9)/2H S,2HAX,2HPY/
      DATA L(1,10),L(2,10),L(3,10)/2H D,2HAX,2HPY/
      DATA L(1,11),L(2,11),L(3,11)/2H C,2HAX,2HPY/
      DATA L(1,12),L(2,12),L(3,12)/2H S,2HRO,2HTG/
      DATA L(1,13),L(2,13),L(3,13)/2H D,2HRO,2HTG/
      DATA L(1,14),L(2,14),L(3,14)/2H  ,2HSR,2HOT/
      DATA L(1,15),L(2,15),L(3,15)/2H  ,2HDR,2HOT/
      DATA L(1,16),L(2,16),L(3,16)/2HSR,2HOT,2HMG/
      DATA L(1,17),L(2,17),L(3,17)/2HDR,2HOT,2HMG/
      DATA L(1,18),L(2,18),L(3,18)/2H S,2HRO,2HTM/
      DATA L(1,19),L(2,19),L(3,19)/2H D,2HRO,2HTM/
      DATA L(1,20),L(2,20),L(3,20)/2H S,2HCO,2HPY/
      DATA L(1,21),L(2,21),L(3,21)/2H D,2HCO,2HPY/
      DATA L(1,22),L(2,22),L(3,22)/2H C,2HCO,2HPY/
      DATA L(1,23),L(2,23),L(3,23)/2H S,2HSW,2HAP/
      DATA L(1,24),L(2,24),L(3,24)/2H D,2HSW,2HAP/
      DATA L(1,25),L(2,25),L(3,25)/2H C,2HSW,2HAP/
      DATA L(1,26),L(2,26),L(3,26)/2H S,2HNR,2HM2/
      DATA L(1,27),L(2,27),L(3,27)/2H D,2HNR,2HM2/
      DATA L(1,28),L(2,28),L(3,28)/2HSC,2HNR,2HM2/
      DATA L(1,29),L(2,29),L(3,29)/2H S,2HAS,2HUM/
      DATA L(1,30),L(2,30),L(3,30)/2H D,2HAS,2HUM/
      DATA L(1,31),L(2,31),L(3,31)/2HSC,2HAS,2HUM/
      DATA L(1,32),L(2,32),L(3,32)/2H S,2HSC,2HAL/
      DATA L(1,33),L(2,33),L(3,33)/2H D,2HSC,2HAL/
      DATA L(1,34),L(2,34),L(3,34)/2H C,2HSC,2HAL/
      DATA L(1,35),L(2,35),L(3,35)/2HCS,2HSC,2HAL/
      DATA L(1,36),L(2,36),L(3,36)/2HIS,2HAM,2HAX/
      DATA L(1,37),L(2,37),L(3,37)/2HID,2HAM,2HAX/
      DATA L(1,38),L(2,38),L(3,38)/2HIC,2HAM,2HAX/
C
      IF (KPRINT.GE.2) WRITE(NPRINT,1000)ICASE,(L(I,ICASE),I = 1, 3)
      RETURN
C
 1000 FORMAT('0TEST OF SUBPROGRAM NO.',I3,2X,3A2)
      END
      SUBROUTINE DTEST(LEN,DCOMP,DTRUE,DSIZE,DFAC,KPRINT)
C1    ********************************* DTEST **************************
C
C     THIS SUBR COMPARES ARRAYS  DCOMP() AND DTRUE() OF LENGTH LEN TO
C     SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY DFAC, ARE
C     NEGLIGIBLE.
C
C     C. L. LAWSON, JPL, 1974 DEC 10
C2
      COMMON/COMBLA/NPRINT,ICASE,N,INCX,INCY,MODE,PASS
      LOGICAL          PASS
      DOUBLE PRECISION DCOMP(LEN),DTRUE(LEN),DSIZE(LEN),DFAC,DDIFF,DD
C
        DO 10 I = 1, LEN
        DD = DCOMP(I)-DTRUE(I)
        IF(DDIFF(DABS(DSIZE(I))+DABS(DFAC*DD),DABS(DSIZE(I))) .EQ. 0.D0)
     *      GO TO 10
C
C                             HERE DCOMP(I) IS NOT CLOSE TO DTRUE(I).
C
        IF(.NOT. PASS) GO TO 5
C                             PRINT FAIL MESSAGE AND HEADER.
        PASS = .FALSE.
        IF (KPRINT.LT.2) GO TO 10
        WRITE(NPRINT,1000)
        WRITE(NPRINT,1001)
    5   IF (KPRINT.GE.2) WRITE(NPRINT,1002)ICASE,N,INCX,INCY,MODE,I,
     *                      DCOMP(I),DTRUE(I),DD,DSIZE(I)
   10   CONTINUE
      RETURN
 1000 FORMAT('+',39X,'FAIL')
 1001 FORMAT('0CASE  N INCX INCY MODE  I',
     1       29X,'COMP(I)',29X,'TRUE(I)',2X,'DIFFERENCE',
     2       5X,'SIZE(I)'/1X)
 1002 FORMAT(1X,I4,I3,3I5,I3,2D36.18,2D12.4)
      END
      SUBROUTINE ITEST(LEN,ICOMP,ITRUE,KPRINT)
C1    ********************************* ITEST **************************
C
C     THIS SUBROUTINE COMPARES THE ARRAYS ICOMP() AND ITRUE() OF
C     LENGTH LEN FOR EQUALITY.
C     C. L. LAWSON, JPL, 1974 DEC 10
C2
      COMMON/COMBLA/NPRINT,ICASE,N,INCX,INCY,MODE,PASS
      LOGICAL          PASS
      INTEGER          ICOMP(LEN), ITRUE(LEN)
C
        DO 10 I = 1, LEN
        IF(ICOMP(I) .EQ. ITRUE(I)) GO TO 10
C
C                            HERE ICOMP(I) IS NOT EQUAL TO ITRUE(I).
C
        IF(.NOT. PASS) GO TO 5
C                             PRINT FAIL MESSAGE AND HEADER.
        PASS = .FALSE.
      IF (KPRINT.LT.2) GO TO 2
         WRITE(NPRINT,1000)
        WRITE(NPRINT,1001)
    2 CONTINUE
    5   ID=ICOMP(I)-ITRUE(I)
      IF (KPRINT.LT.2) GO TO 10
        WRITE(NPRINT,1002) ICASE,N,INCX,INCY,MODE,I,ICOMP(I),ITRUE(I),ID
  10    CONTINUE
      RETURN
 1000 FORMAT('+',39X,'FAIL')
 1001 FORMAT('0CASE  N INCX INCY MODE  I',
     1       29X,'COMP(I)',29X,'TRUE(I)',2X,'DIFFERENCE'/1X)
 1002 FORMAT(1X,I4,I3,3I5,I3,2I36,I12)
      END
      SUBROUTINE STEST(LEN,SCOMP,STRUE,SSIZE,SFAC,KPRINT)
C1    ********************************* STEST **************************
C
C     THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
C     SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY SFAC, ARE
C     NEGLIGIBLE.
C
C     C. L. LAWSON, JPL, 1974 DEC 10
C2
      REAL             SCOMP(LEN),STRUE(LEN),SSIZE(LEN),SFAC,SDIFF,SD
      LOGICAL          PASS
      COMMON/COMBLA/NPRINT,ICASE,N,INCX,INCY,MODE,PASS
C
         DO 10 I = 1, LEN
         SD = SCOMP(I)-STRUE(I)
         IF( SDIFF(ABS(SSIZE(I))+ABS(SFAC*SD), ABS(SSIZE(I))) .EQ. 0.)
     *      GO TO 10
C
C                             HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I).
C
         IF(.NOT. PASS) GO TO 5
C                             PRINT FAIL MESSAGE AND HEADER.
         PASS = .FALSE.
         IF (KPRINT.LT.2) GO TO 10
         WRITE(NPRINT,1000)
         WRITE(NPRINT,1001)
         PASS = .FALSE.
    5    IF (KPRINT.GE.2)WRITE(NPRINT,1002)ICASE,N,INCX,INCY,MODE,I,
     *                      SCOMP(I),STRUE(I),SD,SSIZE(I)
   10    CONTINUE
      RETURN
 1000 FORMAT('+',39X,'FAIL')
 1001 FORMAT('0CASE  N INCX INCY MODE  I',
     1       29X,'COMP(I)',29X,'TRUE(I)',2X,'DIFFERENCE',
     2       5X,'SIZE(I)'/1X)
 1002 FORMAT(1X,I4,I3,3I5,I3,2E36.8,2E12.4)
      END
      DOUBLE PRECISION FUNCTION DDIFF(DA,DB)
C1    ********************************* DDIFF **************************
C     COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15
C2
      DOUBLE PRECISION DA,DB
      DDIFF=DA-DB
      RETURN
      END
      FUNCTION SDIFF(SA,SB)
C1    ********************************* SDIFF **************************
C     COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15
C2
      SDIFF=SA-SB
      RETURN
      END
