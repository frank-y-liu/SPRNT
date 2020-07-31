 
      REAL FUNCTION RNOR(JD)
C***BEGIN PROLOGUE  RNOR
C***DATE WRITTEN   810915
C***REVISION DATE  830805
C***CATEGORY NO.  L6A14
C***KEYWORDS  RANDOM NUMBERS, UNIFORM RANDOM NUMBERS
C***AUTHOR    KAHANER, DAVID, SCIENTIFIC COMPUTING DIVISION, NBS
C             MARSAGLIA, GEORGE, COMPUTER SCIENCE DEPT., WASH STATE UNIV
C
C***PURPOSE  GENERATES QUASI NORMAL RANDOM NUMBERS, WITH MEAN ZERO AND
C             UNIT STANDARD DEVIATION, AND CAN BE USED WITH ANY COMPUTER
C             WITH INTEGERS AT LEAST AS LARGE AS 32767.
C***DESCRIPTION
C
C       RNOR generates quasi normal random numbers with zero mean and
C       unit standard deviation.
C       It can be used with any computer with integers at least as
C       large as 32767.
C
C
C   Use
C       First time....
C                   Z = RNOR(JD)
C                     Here JD is any  n o n - z e r o  integer.
C                     This causes initialization of the program
C                     and the first random number to be returned as Z.
C       Subsequent times...
C                   Z = RNOR(0)
C                     Causes the next random number to be returned as Z.
C
C.....................................................................
C
C    Note: Users who wish to transport this program to other
C           computers should read the following ....
C
C   Machine dependencies...
C      MDIG = A lower bound on the number of binary digits available
C              for representing integers, including the sign bit.
C              This must be at least 16, but can be increased in
C              line with remark A below.
C
C   Remarks...
C     A. This program can be used in two ways:
C        (1) To obtain repeatable results on different computers,
C            set 'MDIG' to the smallest of its values on each, or,
C        (2) To allow the longest sequence of random numbers to be
C            generated without cycling (repeating) set 'MDIG' to the
C            largest possible value.
C     B. The sequence of numbers generated depends on the initial
C          input 'JD' as well as the value of 'MDIG'.
C          If MDIG=16 one should find that
C            the first evaluation
C              Z=RNOR(87) gives  Z=-.40079207...
C            The second evaluation
C              Z=RNOR(0) gives   Z=-1.8728870...
C            The third evaluation
C              Z=RNOR(0) gives   Z=1.8216004...
C            The fourth evaluation
C              Z=RNOR(0) gives   Z=.69410355...
C            The thousandth evaluation
C              Z=RNOR(0) gives   Z=.96782424...
C
C***REFERENCES  MARSAGLIA & TSANG, "A FAST, EASILY IMPLEMENTED
C                 METHOD FOR SAMPLING FROM DECREASING OR
C                 SYMMETRIC UNIMODAL DENSITY FUNCTIONS", TO BE
C                 PUBLISHED IN SIAM J SISC 1983.
C***ROUTINES CALLED  I1MACH,XERROR
C***END PROLOGUE  RNOR
      REAL V(65),W(65)
      INTEGER M(17)
      SAVE I1,J1,M,M1,M2,RMAX
      DATA AA,B,C,RMAX/12.37586,.4878992,12.67706,3.0518509E-5/
      DATA C1,C2,PC,XN/.9689279,1.301198,.1958303E-1,2.776994/
      DATA V/ .3409450, .4573146, .5397793, .6062427, .6631691
     +, .7136975, .7596125, .8020356, .8417227, .8792102, .9148948
     +, .9490791, .9820005, 1.0138492, 1.0447810, 1.0749254, 1.1043917
     +,1.1332738, 1.1616530, 1.1896010, 1.2171815, 1.2444516, 1.2714635
     +,1.2982650, 1.3249008, 1.3514125, 1.3778399, 1.4042211, 1.4305929
     +,1.4569915, 1.4834526, 1.5100121, 1.5367061, 1.5635712, 1.5906454
     +,1.6179680, 1.6455802, 1.6735255, 1.7018503, 1.7306045, 1.7598422
     +,1.7896223, 1.8200099, 1.8510770, 1.8829044, 1.9155830, 1.9492166
     +,1.9839239, 2.0198430, 2.0571356, 2.0959930, 2.1366450, 2.1793713
     +,2.2245175, 2.2725185, 2.3239338, 2.3795007, 2.4402218, 2.5075117
     +,2.5834658, 2.6713916, 2.7769943, 2.7769943, 2.7769943, 2.7769943/
      DATA W/   .10405134E-04, .13956560E-04, .16473259E-04,
     + .18501623E-04, .20238931E-04, .21780983E-04, .23182241E-04,
     + .24476931E-04, .25688121E-04, .26832186E-04, .27921226E-04,
     + .28964480E-04, .29969191E-04, .30941168E-04, .31885160E-04,
     + .32805121E-04, .33704388E-04, .34585827E-04, .35451919E-04,
     + .36304851E-04, .37146564E-04, .37978808E-04, .38803170E-04,
     + .39621114E-04, .40433997E-04, .41243096E-04, .42049621E-04,
     + .42854734E-04, .43659562E-04, .44465208E-04, .45272764E-04,
     + .46083321E-04, .46897980E-04, .47717864E-04, .48544128E-04,
     + .49377973E-04, .50220656E-04, .51073504E-04, .51937936E-04,
     + .52815471E-04, .53707761E-04, .54616606E-04, .55543990E-04,
     + .56492112E-04, .57463436E-04, .58460740E-04, .59487185E-04,
     + .60546402E-04, .61642600E-04, .62780711E-04, .63966581E-04,
     + .65207221E-04, .66511165E-04, .67888959E-04, .69353880E-04,
     + .70922996E-04, .72618816E-04, .74471933E-04, .76525519E-04,
     + .78843526E-04, .81526890E-04, .84749727E-04,
     + .84749727E-04, .84749727E-04, .84749727E-04/
      DATA M(1),M(2),M(3),M(4),M(5),M(6),M(7),M(8),M(9),M(10),M(11),
     1     M(12),M(13),M(14),M(15),M(16),M(17)
     2                   / 30788,23052,2053,19346,10646,19427,23975,
     3                     19049,10949,19693,29746,26748,2796,23890,
     4                     29168,31924,16499 /
      DATA M1,M2,I1,J1 / 32767,256,5,17 /
C Fast part...
C
C
C***FIRST EXECUTABLE STATEMENT  RNOR
      IF(JD.NE.0)GO TO 27
   10 CONTINUE
      I=M(I1)-M(J1)
      IF(I .LT. 0) I=I+M1
      M(J1)=I
      I1=I1-1
      IF(I1 .EQ. 0) I1=17
      J1=J1-1
      IF(J1 .EQ. 0) J1=17
      J=MOD(I,64)+1
      RNOR=I*W(J+1)
      IF( ( (I/M2)/2 )*2.EQ.(I/M2))RNOR=-RNOR
      IF(ABS(RNOR).LE.V(J))RETURN
C Slow part; AA is a*f(0)
      X=(ABS(RNOR)-V(J))/(V(J+1)-V(J))
      Y=UNI(0)
      S=X+Y
      IF(S.GT.C2)GO TO 11
      IF(S.LE.C1)RETURN
      IF(Y.GT.C-AA*EXP(-.5*(B-B*X)**2))GO TO 11
      IF(EXP(-.5*V(J+1)**2)+Y*PC/V(J+1).LE.EXP(-.5*RNOR**2))RETURN
C Tail part; 3.855849 is .5*XN**2
   22 S=XN-ALOG(UNI(0))/XN
      IF(3.855849+ALOG(UNI(0))-XN*S.GT.-.5*S**2)GO TO 22
      RNOR=SIGN(S,RNOR)
      RETURN
   11 RNOR=SIGN(B-B*X,RNOR)
      RETURN
C  FILL
   27 CONTINUE
      MDIG=I1MACH(8)+1
C          BE SURE THAT MDIG AT LEAST 16...
      IF(MDIG.LT.16)CALL XERROR('RNOR--MDIG LESS THAN 16',23,1,2)
      M1 = 2**(MDIG-2) + (2**(MDIG-2)-1)
      M2 = 2**(MDIG/2)
      JSEED = MIN0(IABS(JD),M1)
      IF( MOD(JSEED,2).EQ.0 ) JSEED=JSEED-1
      K0 =MOD(9069,M2)
      K1 = 9069/M2
      J0 = MOD(JSEED,M2)
      J1 = JSEED/M2
      DO 2 I=1,17
        JSEED = J0*K0
        J1 = MOD(JSEED/M2+J0*K1+J1*K0,M2/2)
        J0 = MOD(JSEED,M2)
    2   M(I) = J0+M2*J1
      J1=17
      I1=5
      RMAX = 1./FLOAT(M1)
C        Seed uniform (0,1) generator.  (Just a dummy call)
      RNOR=UNI(JD)
      DO 28 I=1,65
   28  W(I)=RMAX*V(I)
      GO TO 10
      END
