C
C   DRIVER FOR TESTING CMLIB ROUTINES
C     VCOST   VCOSTI
C     VSINT   VSINTI
C     VCOSQF  VCOSQB  VCOSQI
C     VSINQF  VSINQB  VSINQI
C     VRFFTF  VRFFTB  VRFFTI
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
      CALL XERMAX(1000)
C   TEST ALL ROUTINES
      CALL FFTQX1(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C
      IF(KPRINT.GE.1.AND.ITEST.NE.1) WRITE(LUN,2)
2     FORMAT(/' ***** WARNING -- AT LEAST ONE TEST FOR SUBLIBRARY  VFFT
     1HAS FAILED ***** ')
      IF(KPRINT.GE.1.AND.ITEST.EQ.1) WRITE(LUN,3)
3     FORMAT(/' ----- SUBLIBRARY VFFT PASSED ALL TESTS ----- ')
      END
      REAL FUNCTION R2MACH(I)
      COMMON/XXMULT/TIMES
      R2MACH=R1MACH(I)
      IF(I.EQ.1.OR. I.EQ.2) RETURN
      R2MACH = R2MACH * TIMES
      RETURN
      END
      SUBROUTINE FFTQX1(KPRINT,IPASS)
C
C  ------------------------------------------------------
C  QUICKTEST FOR VECTORIZED MULTIPLE FOURIER TRANSFORMS
C  ------------------------------------------------------
C
      PARAMETER (MAXN=20, MAXSEQ=20, NWORK=3*MAXN+15)
      REAL X(MAXSEQ,MAXN+1), XT(MAXSEQ,MAXN+1), WORK(NWORK)
      REAL SAVE(MAXSEQ,MAXN+1)
      COMMON/UNIT/LUN
C
      EXTERNAL VSINT,VCOST,VSINQF,VSINQB,VCOSQF,VCOSQB,
     +         VSINTI,VCOSTI,VSINQI,VCOSQI
C
      IPASS = 1
C
      T = UNI(2963)
      DO 10 J=1,MAXN
         DO 10 I=1,MAXSEQ
            X(I,J) = UNI(0)
   10 CONTINUE
C
      IF (KPRINT .EQ. 3) WRITE(LUN,2000)
C
      CALL TESTS('VCOST/VCOSTI',VCOSTI,VCOST,VCOST,
     +           KPRINT,IPASS,MAXN,MAXSEQ,X,XT,WORK,SAVE)
      CALL TESTS('VSINT/VSINTI',VSINTI,VSINT,VSINT,
     +           KPRINT,IPASS,MAXN,MAXSEQ,X,XT,WORK,SAVE)
      CALL TESTS('VCOSQF/VCOSQB/VCOSQI',VCOSQI,VCOSQF,VCOSQB,
     +           KPRINT,IPASS,MAXN,MAXSEQ,X,XT,WORK,SAVE)
      CALL TESTS('VSINQF/VSINQB/VSINQI',VSINQI,VSINQF,VSINQB,
     +           KPRINT,IPASS,MAXN,MAXSEQ,X,XT,WORK,SAVE)
C
      RETURN
C
 2000 FORMAT('1*****************************************************'
     +      /' *                                                   *'
     +      /' *        QUICKCHECK FOR VFFT SUBLIBRARY             *'
     +      /' *                                                   *'
     +      /' *****************************************************'
     +      / )
      END
      SUBROUTINE TESTS (NAME,FFTI,FFTF,FFTB,KPRINT,IPASS,
     +                  MAXN,MAXSEQ,X,XT,WORK,SAVE)
C
C-------------------------------------------------------------------
C  RUN A SET OF TESTS ON A FFTI/FFTF/FFTB TRIPLE
C-------------------------------------------------------------------
C
      LOGICAL FIRST
      CHARACTER*(*) NAME
      REAL X(MAXSEQ,*), XT(MAXSEQ,*), SAVE(MAXSEQ,*), WORK(*)
      EXTERNAL FFTI, FFTF, FFTB
      COMMON /UNIT/ LUN
C
      EPS = 100.0E0*R2MACH(4)
      FIRST = .TRUE.
      LNAME = LEN(NAME)
C
      IF (KPRINT .EQ. 3) WRITE(LUN,2000) NAME(1:LNAME)
 
      DO 500 N=1,MAXN
         DO 500 M=1,MAXSEQ
            CALL TEST(FFTI,FFTF,FFTB,X,MAXSEQ,M,N,XT,WORK,SAVE,ERROR)
            IF (ERROR .LT. EPS) THEN
C
C              ... TEST PASSED
C
               IF (KPRINT .EQ. 3) WRITE(LUN,2001) M,N,ERROR,'PASSED'
            ELSE
C
C              ... TEST FAILED
C
               IPASS = 0
               IF ((KPRINT .EQ. 2) .AND. FIRST) THEN
                  FIRST = .FALSE.
                  WRITE(LUN,2002) NAME(1:LNAME)
               ENDIF
               IF (KPRINT .GE. 2) WRITE(LUN,2001) M,N,ERROR,'FAILED'
            ENDIF
  500 CONTINUE
C
      RETURN
C
 2000 FORMAT(/' ****************************************************'
     +       /'  BEGINNING TEST OF ',A
     +       /' ****************************************************' /
     +       /'   NUMBER OF   LENGTH OF'
     +       /'   SEQUENCES   SEQUENCES       ERROR       STATUS'
     +       /' ----------------------------------------------------' /)
 2001 FORMAT(4X,I4,8X,I4,9X,1PE10.4,4X,A6)
 2002 FORMAT(/' ****************************************************'
     +       /'  FAILURES IN TESTS OF ',A
     +       /' ****************************************************' /
     +       /'   NUMBER OF   LENGTH OF'
     +       /'   SEQUENCES   SEQUENCES       ERROR       STATUS'
     +       /' ----------------------------------------------------' /)
      END
      SUBROUTINE TEST (FFTI,FFTF,FFTB,X,MDIM,M,N,XT,WORK,SAVE,ERROR)
C---------------------------------------------------------------------
C  TEST ONE SET OF FFT ROUTINES ON ONE DATASET
C---------------------------------------------------------------------
      REAL X(MDIM,*), XT(MDIM,*), SAVE(MDIM,*), WORK(*)
      DO 10 J=1,N
         DO 10 I=1,M
            SAVE(I,J) = X(I,J)
   10 CONTINUE
      CALL FFTI(N,WORK)
      CALL FFTF(M,N,SAVE,XT,MDIM,WORK)
      CALL FFTB(M,N,SAVE,XT,MDIM,WORK)
      DO 20 J=1,N
         DO 20 I=1,M
            SAVE(I,J) = ABS(SAVE(I,J) - X(I,J))
   20 CONTINUE
      ERROR = 0.0
      DO 30 J=1,N
         DO 30 I=1,M
            ERROR = MAX(SAVE(I,J),ERROR)
   30 CONTINUE
      RETURN
      END
