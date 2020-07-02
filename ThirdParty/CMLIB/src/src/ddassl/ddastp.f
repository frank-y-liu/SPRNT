      SUBROUTINE DDASTP(X,Y,YPRIME,NEQ,
     *  RES,JAC,H,WT,JSTART,IDID,RPAR,IPAR,
     *  PHI,DELTA,E,WM,IWM,
     *  ALPHA,BETA,GAMMA,PSI,SIGMA,
     *  CJ,CJOLD,HOLD,S,HMIN,UROUND,
     *  IPHASE,JCALC,K,KOLD,NS,NONNEG)
C
C***BEGIN PROLOGUE  DDASTP
C***REFER TO  DDASSL
C***ROUTINES CALLED  DDANRM,DDAJAC,DDASLV,DDATRP
C***COMMON BLOCKS    DDA001
C***DATE WRITTEN   830315   (YYMMDD)
C***REVISION DATE  830315   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***END PROLOGUE  DDASTP
C
C
C-----------------------------------------------------------------------
C     dastep solves a system of differential/
C     algebraic equations of the form
C     g(x,y,yprime) = 0,  for one step (normally
C     from x to x+h).
C
C     the methods used are modified divided
C     difference,fixed leading coefficient
C     forms of backward differentiation
C     formulas. the code adjusts the stepsize
C     and order to control the local error per
C     step.
C
C
C     the parameters represent
C     x  --        independent variable
C     y  --        solution vector at x
C     yprime --    derivative of solution vector
C                  after successful step
C     neq --       number of equations to be integrated
C     res --       external user-supplied subroutine
C                  to evaluate the residual.  the call is
C                  call res(x,y,yprime,delta,ires,rpar,ipar)
C                  x,y,yprime are input.  delta is output.
C                  on input, ires=0.  res should alter ires only
C                  if it encounters an illegal value of y or a
C                  stop condition.  set ires=-1 if an input value
C                  of y is illegal, and dastep will try to solve
C                  the problem without getting ires = -1.  if
C                  ires=-2, dastep returns control to the calling
C                  program with idid = -11.
C     jac --       external user-supplied routine to evaluate
C                  the iteration matrix (this is optional)
C                  the call is of the form
C                  call jac(x,y,yprime,pd,cj,rpar,ipar)
C                  pd is the matrix of partial derivatives,
C                  pd=dg/dy+cj*dg/dyprime
C     h --         appropriate step size for next step.
C                  normally determined by the code
C     wt --        vector of weights for error criterion.
C     jstart --    integer variable set 0 for
C                  first step, 1 otherwise.
C     idid --      completion code with the following meanings%
C                  idid= 1 -- the step was completed successfully
C                  idid=-6 -- the error test failed repeatedly
C                  idid=-7 -- the corrector could not converge
C                  idid=-8 -- the iteration matrix is singular
C                  idid=-9 -- the corrector could not converge.
C                             there were repeated error test
C                             failures on this step.
C                  idid=-10-- the corrector could not converge
C                             because ires was equal to minus one
C                  idid=-11-- ires equal to -2 was encountered,
C                             and control is being returned to
C                             the calling program
C     rpar,ipar -- real and integer parameter arrays that
C                  are used for communication between the
C                  calling program and external user routines
C                  they are not altered by dastep
C     phi --       array of divided differences used by
C                  dastep. the length is neq*(k+1),where
C                  k is the maximum order
C     delta,e --   work vectors for dastep of length neq
C     wm,iwm --    real and integer arrays storing
C                  matrix information such as the matrix
C                  of partial derivatives,permutation
C                  vector,and various other information.
C
C     the other parameters are information
C     which is needed internally by dastep to
C     continue from step to step.
C
C-----------------------------------------------------------------------
C
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL CONVGD
      DIMENSION Y(*),YPRIME(*),WT(*)
      DIMENSION PHI(NEQ,*),DELTA(*),E(*)
      DIMENSION WM(*),IWM(*)
      DIMENSION PSI(*),ALPHA(*),BETA(*),GAMMA(*),SIGMA(*)
      DIMENSION RPAR(*),IPAR(*)
      EXTERNAL RES,JAC
      COMMON/DDA001/NPD,NTEMP,
     *   LML,LMU,LMXORD,LMTYPE,
     *   LNST,LNRE,LNJE,LETF,LCTF,LIPVT
      DATA MAXIT/4/
      DATA XRATE/0.25D0/
C
C
C
C
C
C-----------------------------------------------------------------------
C     block 1.
C     initialize. on the first call,set
C     the order to 1 and initialize
C     other variables.
C-----------------------------------------------------------------------
C
C     initializations for all calls
      IDID=1
      XOLD=X
      NCF=0
      NSF=0
      NEF=0
      IF(JSTART .NE. 0) GO TO 120
C
C     if this is the first step,perform
C     other initializations
      IWM(LETF) = 0
      IWM(LCTF) = 0
      K=1
      KOLD=0
      HOLD=0.0D0
      JSTART=1
      PSI(1)=H
      CJOLD = 1.0D0/H
      CJ = CJOLD
      S = 100.D0
      JCALC = -1
      DELNRM=1.0D0
      IPHASE = 0
      NS=0
120   CONTINUE
C
C
C
C
C
C-----------------------------------------------------------------------
C     block 2
C     compute coefficients of formulas for
C     this step.
C-----------------------------------------------------------------------
200   CONTINUE
      KP1=K+1
      KP2=K+2
      KM1=K-1
      XOLD=X
      IF(H.NE.HOLD.OR.K .NE. KOLD) NS = 0
      NS=MIN0(NS+1,KOLD+2)
      NSP1=NS+1
      IF(KP1 .LT. NS)GO TO 230
C
      BETA(1)=1.0D0
      ALPHA(1)=1.0D0
      TEMP1=H
      GAMMA(1)=0.0D0
      SIGMA(1)=1.0D0
      DO 210 I=2,KP1
         TEMP2=PSI(I-1)
         PSI(I-1)=TEMP1
         BETA(I)=BETA(I-1)*PSI(I-1)/TEMP2
         TEMP1=TEMP2+H
         ALPHA(I)=H/TEMP1
         SIGMA(I)=DBLE(I-1)*SIGMA(I-1)*ALPHA(I)
         GAMMA(I)=GAMMA(I-1)+ALPHA(I-1)/H
210      CONTINUE
      PSI(KP1)=TEMP1
230   CONTINUE
C
C     compute alphas, alpha0
      ALPHAS = 0.0D0
      ALPHA0 = 0.0D0
      DO 240 I = 1,K
        ALPHAS = ALPHAS - 1.0D0/DBLE(I)
        ALPHA0 = ALPHA0 - ALPHA(I)
240     CONTINUE
C
C     compute leading coefficient cj
      CJLAST = CJ
      CJ = -ALPHAS/H
C
C     compute variable stepsize error coefficient ck
      CK = DABS(ALPHA(KP1) + ALPHAS - ALPHA0)
      CK = DMAX1(CK,ALPHA(KP1))
C
C     decide whether new jacobian is needed
      TEMP1 = (1.0D0 - XRATE)/(1.0D0 + XRATE)
      TEMP2 = 1.0D0/TEMP1
      IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1
      IF (CJ .NE. CJLAST) S = 100.D0
C
C     change phi to phi star
      IF(KP1 .LT. NSP1) GO TO 280
      DO 270 J=NSP1,KP1
         DO 260 I=1,NEQ
260         PHI(I,J)=BETA(J)*PHI(I,J)
270      CONTINUE
280   CONTINUE
C
C     update time
      X=X+H
C
C
C
C
C
C-----------------------------------------------------------------------
C     block 3
C     predict the solution and derivative,
C     and solve the corrector equation
C-----------------------------------------------------------------------
C
C     first,predict the solution and derivative
300   CONTINUE
      DO 310 I=1,NEQ
         Y(I)=PHI(I,1)
310      YPRIME(I)=0.0D0
      DO 330 J=2,KP1
         DO 320 I=1,NEQ
            Y(I)=Y(I)+PHI(I,J)
320         YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)
330   CONTINUE
      PNORM = DDANRM (NEQ,Y,WT,RPAR,IPAR)
C
C
C
C     solve the corrector equation using a
C     modified newton scheme.
      CONVGD= .TRUE.
      M=0
      IWM(LNRE)=IWM(LNRE)+1
      IRES = 0
      CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
C
C
C     if indicated,reevaluate the
C     iteration matrix pd = dg/dy + cj*dg/dyprime
C     (where g(x,y,yprime)=0). set
C     jcalc to 0 as an indicator that
C     this has been done.
      IF(JCALC .NE. -1)GO TO 340
      IWM(LNJE)=IWM(LNJE)+1
      JCALC=0
      CALL DDAJAC(NEQ,X,Y,YPRIME,DELTA,CJ,H,
     * IER,WT,E,WM,IWM,RES,IRES,UROUND,JAC,RPAR,IPAR)
      CJOLD=CJ
      S = 100.D0
      IF (IRES .LT. 0) GO TO 380
      IF(IER .NE. 0)GO TO 380
      NSF=0
C
C
C     initialize the error accumulation vector e.
340   CONTINUE
      DO 345 I=1,NEQ
345      E(I)=0.0D0
C
      S = 100.E0
C
C
C     corrector loop.
350   CONTINUE
C
C     multiply residual by temp1 to accelerate convergence
      TEMP1 = 2.0D0/(1.0D0 + CJ/CJOLD)
      DO 355 I = 1,NEQ
355     DELTA(I) = DELTA(I) * TEMP1
C
C     compute a new iterate (back-substitution).
C     store the correction in delta.
      CALL DDASLV(NEQ,DELTA,WM,IWM)
C
C     update y,e,and yprime
      DO 360 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
         E(I)=E(I)-DELTA(I)
360      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
C
C     test for convergence of the iteration
      DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM .LE. 100.D0*UROUND*PNORM) GO TO 375
      IF (M .GT. 0) GO TO 365
         OLDNRM = DELNRM
         GO TO 367
365   RATE = (DELNRM/OLDNRM)**(1.0D0/DBLE(M))
      IF (RATE .GT. 0.90D0) GO TO 370
      S = RATE/(1.0D0 - RATE)
367   IF (S*DELNRM .LE. 0.33D0) GO TO 375
C
C     the corrector has not yet converged.
C     update m and test whether the
C     maximum number of iterations have
C     been tried.
      M=M+1
      IF(M.GE.MAXIT)GO TO 370
C
C     evaluate the residual
C     and go back to do another iteration
      IWM(LNRE)=IWM(LNRE)+1
      IRES = 0
      CALL RES(X,Y,YPRIME,DELTA,IRES,
     *  RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
      GO TO 350
C
C
C     the corrector failed to converge in maxit
C     iterations. if the iteration matrix
C     is not current,re-do the step with
C     a new iteration matrix.
370   CONTINUE
      IF(JCALC.EQ.0)GO TO 380
      JCALC=-1
      GO TO 300
C
C
C     the iteration has converged.  if nonnegativity of solution is
C     required, set the solution nonnegative, if the perturbation
C     to do it is small enough.  if the change is too large, then
C     consider the corrector iteration to have failed.
375   IF(NONNEG .EQ. 0) GO TO 390
      DO 377 I = 1,NEQ
377      DELTA(I) = DMIN1(Y(I),0.0D0)
      DELNRM = DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF(DELNRM .GT. 0.33D0) GO TO 380
      DO 378 I = 1,NEQ
378      E(I) = E(I) - DELTA(I)
      GO TO 390
C
C
C     exits from block 3
C     no convergence with current iteration
C     matrix,or singular iteration matrix
380   CONVGD= .FALSE.
390   JCALC = 1
      IF(.NOT.CONVGD)GO TO 600
C
C
C
C
C
C-----------------------------------------------------------------------
C     block 4
C     estimate the errors at orders k,k-1,k-2
C     as if constant stepsize was used. estimate
C     the local error at order k and test
C     whether the current step is successful.
C-----------------------------------------------------------------------
C
C     estimate errors at orders k,k-1,k-2
      ENORM = DDANRM(NEQ,E,WT,RPAR,IPAR)
      ERK = SIGMA(K+1)*ENORM
      TERK = FLOAT(K+1)*ERK
      EST = ERK
      KNEW=K
      IF(K .EQ. 1)GO TO 430
      DO 405 I = 1,NEQ
405     DELTA(I) = PHI(I,KP1) + E(I)
      ERKM1=SIGMA(K)*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKM1 = FLOAT(K)*ERKM1
      IF(K .GT. 2)GO TO 410
      IF(TERKM1 .LE. 0.5*TERK)GO TO 420
      GO TO 430
410   CONTINUE
      DO 415 I = 1,NEQ
415     DELTA(I) = PHI(I,K) + DELTA(I)
      ERKM2=SIGMA(K-1)*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKM2 = FLOAT(K-1)*ERKM2
      IF(DMAX1(TERKM1,TERKM2).GT.TERK)GO TO 430
C     lower the order
420   CONTINUE
      KNEW=K-1
      EST = ERKM1
C
C
C     calculate the local error for the current step
C     to see if the step was successful
430   CONTINUE
      ERR = CK * ENORM
      IF(ERR .GT. 1.0D0)GO TO 600
C
C
C
C
C
C-----------------------------------------------------------------------
C     block 5
C     the step is successful. determine
C     the best order and stepsize for
C     the next step. update the differences
C     for the next step.
C-----------------------------------------------------------------------
      IDID=1
      IWM(LNST)=IWM(LNST)+1
      KDIFF=K-KOLD
      KOLD=K
      HOLD=H
C
C
C     estimate the error at order k+1 unless%
C        already decided to lower order, or
C        already using maximum order, or
C        stepsize not constant, or
C        order raised in previous step
      IF(KNEW.EQ.KM1.OR.K.EQ.IWM(LMXORD))IPHASE=1
      IF(IPHASE .EQ. 0)GO TO 545
      IF(KNEW.EQ.KM1)GO TO 540
      IF(K.EQ.IWM(LMXORD)) GO TO 550
      IF(KP1.GE.NS.OR.KDIFF.EQ.1)GO TO 550
      DO 510 I=1,NEQ
510      DELTA(I)=E(I)-PHI(I,KP2)
      ERKP1 = (1.0D0/DBLE(K+2))*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKP1 = FLOAT(K+2)*ERKP1
      IF(K.GT.1)GO TO 520
      IF(TERKP1.GE.0.5D0*TERK)GO TO 550
      GO TO 530
520   IF(TERKM1.LE.DMIN1(TERK,TERKP1))GO TO 540
      IF(TERKP1.GE.TERK.OR.K.EQ.IWM(LMXORD))GO TO 550
C
C     raise order
530   K=KP1
      EST = ERKP1
      GO TO 550
C
C     lower order
540   K=KM1
      EST = ERKM1
      GO TO 550
C
C     if iphase = 0, increase order by one and multiply stepsize by
C     factor two
545   K = KP1
      HNEW = H*2.0D0
      H = HNEW
      GO TO 575
C
C
C     determine the appropriate stepsize for
C     the next step.
550   HNEW=H
      TEMP2=K+1
      R=(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      IF(R .LT. 2.0D0) GO TO 555
      HNEW = 2.0D0*H
      GO TO 560
555   IF(R .GT. 1.0D0) GO TO 560
      R = DMAX1(0.5D0,DMIN1(0.9D0,R))
      HNEW = H*R
560   H=HNEW
C
C
C     update differences for next step
575   CONTINUE
      IF(KOLD.EQ.IWM(LMXORD))GO TO 585
      DO 580 I=1,NEQ
580      PHI(I,KP2)=E(I)
585   CONTINUE
      DO 590 I=1,NEQ
590      PHI(I,KP1)=PHI(I,KP1)+E(I)
      DO 595 J1=2,KP1
         J=KP1-J1+1
         DO 595 I=1,NEQ
595      PHI(I,J)=PHI(I,J)+PHI(I,J+1)
      RETURN
C
C
C
C
C
C-----------------------------------------------------------------------
C     block 6
C     the step is unsuccessful. restore x,psi,phi
C     determine appropriate stepsize for
C     continuing the integration, or exit with
C     an error flag if there have been many
C     failures.
C-----------------------------------------------------------------------
600   IPHASE = 1
C
C     restore x,phi,psi
      X=XOLD
      IF(KP1.LT.NSP1)GO TO 630
      DO 620 J=NSP1,KP1
         TEMP1=1.0D0/BETA(J)
         DO 610 I=1,NEQ
610         PHI(I,J)=TEMP1*PHI(I,J)
620      CONTINUE
630   CONTINUE
      DO 640 I=2,KP1
640      PSI(I-1)=PSI(I)-H
C
C
C     test whether failure is due to corrector iteration
C     or error test
      IF(CONVGD)GO TO 660
      IWM(LCTF)=IWM(LCTF)+1
C
C
C     the newton iteration failed to converge with
C     a current iteration matrix.  determine the cause
C     of the failure and take appropriate action.
      IF(IER.EQ.0)GO TO 650
C
C     the iteration matrix is singular. reduce
C     the stepsize by a factor of 4. if
C     this happens three times in a row on
C     the same step, return with an error flag
      NSF=NSF+1
      R = 0.25D0
      H=H*R
      IF (NSF .LT. 3 .AND. DABS(H) .GE. HMIN) GO TO 690
      IDID=-8
      GO TO 675
C
C
C     the newton iteration failed to converge for a reason
C     other than a singular iteration matrix.  if ires = -2, then
C     return.  otherwise, reduce the stepsize and try again, unless
C     too many failures have occured.
650   CONTINUE
      IF (IRES .GT. -2) GO TO 655
      IDID = -11
      GO TO 675
655   NCF = NCF + 1
      R = 0.25D0
      H = H*R
      IF (NCF .LT. 10 .AND. DABS(H) .GE. HMIN) GO TO 690
      IDID = -7
      IF (IRES .LT. 0) IDID = -10
      IF (NEF .GE. 3) IDID = -9
      GO TO 675
C
C
C     the newton scheme converged,and the cause
C     of the failure was the error estimate
C     exceeding the tolerance.
660   NEF=NEF+1
      IWM(LETF)=IWM(LETF)+1
      IF (NEF .GT. 1) GO TO 665
C
C     on first error test failure, keep current order or lower
C     order by one.  compute new stepsize based on differences
C     of the solution.
      K = KNEW
      TEMP2 = K + 1
      R = 0.90D0*(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      R = DMAX1(0.25D0,DMIN1(0.9D0,R))
      H = H*R
      IF (DABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     on second error test failure, use the current order or
C     decrease order by one.  reduce the stepsize by a factor of
C     one quarter.
665   IF (NEF .GT. 2) GO TO 670
      K = KNEW
      H = 0.25D0*H
      IF (DABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     on third and subsequent error test failures, set the order to
C     one and reduce the stepsize by a factor of one quarter
670   K = 1
      H = 0.25D0*H
      IF (DABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C
C
C
C     for all crashes, restore y to its last value,
C     interpolate to find yprime at last x, and return
675   CONTINUE
      CALL DDATRP(X,X,Y,YPRIME,NEQ,K,PHI,PSI)
      RETURN
C
C
C     go back and try this step again
690   GO TO 200
C
C------end of subroutine dastep------
      END
