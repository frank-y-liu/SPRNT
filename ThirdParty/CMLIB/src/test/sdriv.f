C
C   DRIVER FOR TESTING CMLIB ROUTINES
C     SDRIV1
C     SDRIV2
C     SDRIV3
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
C   TEST SDRIV1,SDRIV2,SDRIV3
      CALL SDRIQX(KPRINT,IPASS)
      ITEST=ITEST*IPASS
C
      IF(KPRINT.GE.1.AND.ITEST.NE.1) WRITE(LUN,2)
2     FORMAT(/' ***** WARNING -- AT LEAST ONE TEST FOR SUBLIBRARY SDRIV
     1  HAS FAILED ***** ')
      IF(KPRINT.GE.1.AND.ITEST.EQ.1) WRITE(LUN,3)
3     FORMAT(/' ----- SUBLIBRARY SDRIV PASSED ALL TESTS ----- ')
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
      SUBROUTINE SDRIQX(KPRINT,IPASS)
c
c     Quick Check Program
c
c   Part of SDRIV1,2,3 package, single precision version,
c     to be used with the FORTRAN 77 compiler.
c   D. K. Kahaner, National Bureau of Standards, and
c   c. d. sutherland, los alamos national laboratory, july 5, 1983
c
      common/unit/lun
      external f, jacobn, fa, g
      real eps,ewt(1),hmax,t,tout,work(264),y(3)
      integer iwork(24),kprint,ipass
      data n /3/, ntask /1/, eps /.001e0/, ewt(1) /.1e0/, ierror /3/,
     8     mint /1/, miter /5/, impl /0/, ml /2/, mu /2/, mxord /5/,
     8     hmax /1.e0/, lenw /264/, leniw /24/, nde /0/,
     8     mxstep /10/, iroot /0/
      if(kprint.ge.2) then
      write(lun, '(// '' This program attempts to verify that all '',
     8  ''routines necessary for running'' / '' the sdriv package are'',
     8  '' available.  In order to do this, the error'' / '' handler '',
     8  ''routine must be accessed, leading to what seems to be an '',
     8  ''error'' / '' message. However, this message should be '',
     8  ''ignored if there follows a'' / '' statement to the effect '',
     8  ''that the test has actually passed.''  //)')
      endif
      call xsetf(-1)
      ipass=1
      t = 0.e0
      y(1) = 10.e0
      y(2) = 0.e0
      y(3) = 10.e0
      tout = .001e0
      mstate = 1
      call sdriv1 (n, t, y, tout, mstate, eps, work, lenw)
      t = 0.e0
      y(1) = 10.e0
      y(2) = 0.e0
      y(3) = 10.e0
      tout = .001e0
      mstate = 1
      call sdriv2 (n, t, y, f, tout, mstate, iroot, eps, ewt,
     8             mint, work, lenw, iwork, leniw, g)
      t = 0.e0
      y(1) = 10.e0
      y(2) = 0.e0
      y(3) = 10.e0
      nstate = 1
      tout = 1.e0
      call sdriv3 (n, t, y, f, nstate, tout, ntask, iroot, eps, ewt,
     8             ierror, mint, miter, impl, ml, mu, mxord,
     8             hmax, work, lenw, iwork, leniw, jacobn, fa, nde,
     8             mxstep, g)
      end
      subroutine f (n, t, y, yp)
      real t, y(*), yp(*)
      yp(1) = 1.e0 + (y(2) - y(1)) - y(1)*y(3)
      yp(2) = (y(1) - y(2)) - y(2)*y(3)
      yp(3) = 1.e0 - y(3)*(y(1) + y(2))
      end
      subroutine jacobn (n, t, y, dfdy, matdim, ml, mu)
      i = 1
      end
      subroutine fa (n, t, y, a, matdim, ml, mu, nde)
      i = 1
      end
      subroutine users (y, yh, ywt, save1, save2, t, h, el, impl, n,
     8  nde, iflag)
      i = 1
      end
      real function g (n, t, y, nroot)
      real y(*)
      g = 1.
      end
