      SUBROUTINE SDAINI(X,Y,YPRIME,NEQ,
     *   RES,JAC,H,WT,IDID,RPAR,IPAR,
     *   PHI,DELTA,E,WM,IWM,
     *   HMIN,UROUND,NONNEG)
C
C***BEGIN PROLOGUE  SDAINI
C***REFER TO  SDASSL
C***ROUTINES CALLED  SDANRM,SDAJAC,SDASLV
C***COMMON BLOCKS    SDA001
C***DATE WRITTEN   830315   (YYMMDD)
C***REVISION DATE  830315   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***END PROLOGUE SDAINI
C
C-------------------------------------------------------
C     sdaini takes one step of size h or smaller
C     with the backward euler method, to
C     find yprime at the initial time x. a modified
C     damped newton iteration is used to
C     solve the corrector iteration.
C
C     the initial guess yprime is used in the
C     prediction, and in forming the iteration
C     matrix, but is not involved in the
C     error test. this may have trouble
C     converging if the initial guess is no
C     good, or if g(xy,yprime) depends
C     nonlinearly on yprime.
C
C     the parameters represent:
C     x --         independent variable
C     y --         solution vector at x
C     yprime --    derivative of solution vector
C     neq --       number of equations
C     h --         stepsize. imder may use a stepsize
C                  smaller than h.
C     wt --        vector of weights for error
C                  criterion
C     idid --      completion code with the following meanings
C                  idid= 1 -- yprime was found successfully
C                  idid=-12 -- sdaini failed to find yprime
C     rpar,ipar -- real and integer parameter arrays
C                  that are not altered by sdaini
C     phi --       work space for sdaini
C     delta,e --   work space for sdaini
C     wm,iwm --    real and integer arrays storing
C                  matrix information
C
C-----------------------------------------------------------------
C
C
 
      LOGICAL CONVGD
      DIMENSION Y(*),YPRIME(*),WT(*)
      DIMENSION PHI(NEQ,*),DELTA(*),E(*)
      DIMENSION WM(*),IWM(*)
      DIMENSION RPAR(*),IPAR(*)
      EXTERNAL RES,JAC
      COMMON/SDA001/NPD,NTEMP,
     *  LML,LMU,LMXORD,LMTYPE,
     *  LNST,LNRE,LNJE,LETF,LCTF,LIPVT
 
      DATA MAXIT/10/,MJAC/5/
      DATA DAMP/0.75E0/
 
C
C
C---------------------------------------------------
C     block 1.
C     initializations.
C---------------------------------------------------
C
      IDID=1
      NEF=0
      NCF=0
      NSF=0
      YNORM=SDANRM(NEQ,Y,WT,RPAR,IPAR)
C
C     save y and yprime in phi
      DO 100 I=1,NEQ
         PHI(I,1)=Y(I)
100      PHI(I,2)=YPRIME(I)
 
C
C
C----------------------------------------------------
C     block 2.
C     do one backward euler step.
C----------------------------------------------------
C
C     set up for start of corrector iteration
200   CJ=1.0E0/H
      XNEW=X+H
C
C     predict solution and derivative
 
      DO 250 I=1,NEQ
250     Y(I)=Y(I)+H*YPRIME(I)
C
      JCALC=-1
      M=0
      CONVGD=.TRUE.
C
C
C     corrector loop.
300   IWM(LNRE)=IWM(LNRE)+1
      IRES=0
 
      CALL RES(XNEW,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF (IRES.LT.0) GO TO 430
C
C
C     evaluate the iteration matrix
      IF (JCALC.NE.-1) GO TO 310
      IWM(LNJE)=IWM(LNJE)+1
      JCALC=0
      CALL SDAJAC(NEQ,XNEW,Y,YPRIME,DELTA,CJ,H,
     *   IER,WT,E,WM,IWM,RES,IRES,
     *   UROUND,JAC,RPAR,IPAR)
 
      S=1000000.E0
      IF (IRES.LT.0) GO TO 430
      IF (IER.NE.0) GO TO 430
      NSF=0
 
C
C
C
C     multiply residual by damping factor
310   CONTINUE
      DO 320 I=1,NEQ
320      DELTA(I)=DELTA(I)*DAMP
 
C
C     compute a new iterate (back substitution)
C     store the correction in delta
 
      CALL SDASLV(NEQ,DELTA,WM,IWM)
 
C
C     update y and yprime
 
      DO 330 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
330      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
 
C
C     test for convergence of the iteration.
 
      DELNRM=SDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM.LE.100.E0*UROUND*YNORM)
     *   GO TO 400
 
      IF (M.GT.0) GO TO 340
         OLDNRM=DELNRM
         GO TO 350
 
340   RATE=(DELNRM/OLDNRM)**(1.0E0/FLOAT(M))
      IF (RATE.GT.0.90E0) GO TO 430
      S=RATE/(1.0E0-RATE)
 
350   IF (S*DELNRM .LE. 0.33E0) GO TO 400
C
C
C     the corrector has not yet converged. update
C     m and and test whether the maximum
C     number of iterations have been tried.
C     every mjac iterations, get a new
C     iteration matrix.
 
      M=M+1
      IF (M.GE.MAXIT) GO TO 430
 
      IF ((M/MJAC)*MJAC.EQ.M) JCALC=-1
 
      GO TO 300
 
C
C
C     the iteration has converged.
C     check nonnegativity constraints
400   IF (NONNEG.EQ.0) GO TO 450
      DO 410 I=1,NEQ
410      DELTA(I)=AMIN1(Y(I),0.0E0)
 
      DELNRM=SDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM.GT.0.33E0) GO TO 430
 
      DO 420 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
420      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
      GO TO 450
C
C
C     exits from corrector loop.
430   CONVGD=.FALSE.
450   IF (.NOT.CONVGD) GO TO 600
C
C
C
C-----------------------------------------------------
C     block 3.
C     the corrector iteration converged.
C     do error test.
C-----------------------------------------------------
C
 
      DO 510 I=1,NEQ
510      E(I)=Y(I)-PHI(I,1)
 
      ERR=SDANRM(NEQ,E,WT,RPAR,IPAR)
 
      IF (ERR.LE.1.0E0) RETURN
 
C
C
C
C--------------------------------------------------------
C     block 4.
C     the backward euler step failed. restore y
C     and yprime to their original values.
C     reduce stepsize and try again, if
C     possible.
C---------------------------------------------------------
C
 
600   CONTINUE
      DO 610 I=1,NEQ
         Y(I)=PHI(I,1)
610      YPRIME(I)=PHI(I,2)
 
      IF (CONVGD) GO TO 640
      IF (IER.EQ.0) GO TO 620
         NSF=NSF+1
         H=H*0.25E0
         IF (NSF.LT.3.AND.ABS(H).GE.HMIN) GO TO 690
         IDID=-12
         RETURN
620   IF (IRES.GT.-2) GO TO 630
         IDID=-12
         RETURN
630   NCF=NCF+1
      H=H*0.25E0
      IF (NCF.LT.10.AND.ABS(H).GE.HMIN) GO TO 690
         IDID=-12
         RETURN
 
640   NEF=NEF+1
      R=0.90E0/(2.0E0*ERR+0.0001E0)
      R=AMAX1(0.1E0,AMIN1(0.5E0,R))
      H=H*R
      IF (ABS(H).GE.HMIN.AND.NEF.LT.10) GO TO 690
         IDID=-12
         RETURN
690      GO TO 200
 
C-------------end of subroutine sdaini----------------------
      END
