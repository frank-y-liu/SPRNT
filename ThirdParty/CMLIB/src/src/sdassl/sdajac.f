      SUBROUTINE SDAJAC(NEQ,X,Y,YPRIME,DELTA,CJ,H,
     *  IER,WT,E,WM,IWM,RES,IRES,UROUND,JAC,RPAR,IPAR)
C
C***BEGIN PROLOGUE  SDAJAC
C***REFER TO  SDASSL
C***ROUTINES CALLED  SGEFA,SGBFA
C***COMMON BLOCKS    SDA001
C***DATE WRITTEN   830315   (YYMMDD)
C***REVISION DATE  830315   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***END PROLOGUE  SDAJAC
C-----------------------------------------------------------------------
C     this routine computes the iteration matrix
C     pd=dg/dy+cj*dg/dyprime (where g(x,y,yprime)=0).
C     here pd is computed by the user-supplied
C     routine jac if iwm(mtype) is 1 or 4, and
C     it is computed by numerical finite differencing
C     if iwm(mtype)is 2 or 5
C     the parameters have the following meanings.
C     y        = array containing predicted values
C     yprime   = array containing predicted derivatives
C     delta    = residual evaluated at (x,y,yprime)
C                (used only if iwm(mtype)=2 or 5)
C     cj       = scalar parameter defining iteration matrix
C     h        = current stepsize in integration
C     ier      = variable which is .ne. 0
C                if iteration matrix is singular,
C                and 0 otherwise.
C     wt       = vector of weights for computing norms
C     e        = work space (temporary) of length neq
C     wm       = real work space for matrices. on
C                output it contains the lu decomposition
C                of the iteration matrix.
C     iwm      = integer work space containing
C                matrix information
C     res      = name of the external user-supplied routine
C                to evaluate the residual function g(x,y,yprime)
C     ires     = flag which is equal to zero if no illegal values
C                in res, and less than zero otherwise.  (if ires
C                is less than zero, the matrix was not completed)
C                in this case (if ires .lt. 0), then ier = 0.
C     uround   = the unit roundoff error of the machine being used.
C     jac      = name of the external user-supplied routine
C                to evaluate the iteration matrix (this routine
C                is only used if iwm(mtype) is 1 or 4)
C-----------------------------------------------------------------------
C
      EXTERNAL RES,JAC
      DIMENSION Y(*),YPRIME(*),DELTA(*),WT(*),E(*)
      DIMENSION WM(*),IWM(*),RPAR(*),IPAR(*)
      COMMON/SDA001/NPD,NTEMP,
     *  LML,LMU,LMXORD,LMTYPE,
     *  LNST,LNRE,LNJE,LETF,LCTF,LIPVT
C
      IER = 0
      NPDM1=NPD-1
      MTYPE=IWM(LMTYPE)
      GO TO (100,200,300,400,500),MTYPE
C
C
C     dense user-supplied matrix
100   LENPD=NEQ*NEQ
      DO 110 I=1,LENPD
110      WM(NPDM1+I)=0.0E0
      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      GO TO 230
C
C
C     dense finite-difference-generated matrix
200   IRES=0
      NROW=NPDM1
      SQUR = SQRT(UROUND)
      DO 210 I=1,NEQ
         DEL=SQUR*AMAX1(ABS(Y(I)),ABS(H*YPRIME(I)),
     *     ABS(WT(I)))
         DEL=SIGN(DEL,H*YPRIME(I))
         DEL=(Y(I)+DEL)-Y(I)
         YSAVE=Y(I)
         YPSAVE=YPRIME(I)
         Y(I)=Y(I)+DEL
         YPRIME(I)=YPRIME(I)+CJ*DEL
         CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
         IF (IRES .LT. 0) RETURN
         DELINV=1.0E0/DEL
         DO 220 L=1,NEQ
220      WM(NROW+L)=(E(L)-DELTA(L))*DELINV
      NROW=NROW+NEQ
      Y(I)=YSAVE
      YPRIME(I)=YPSAVE
210   CONTINUE
C
C
C     do dense-matrix lu decomposition on pd
230      CALL SGEFA(WM(NPD),NEQ,NEQ,IWM(LIPVT),IER)
      RETURN
C
C
C     dummy section for iwm(mtype)=3
300   RETURN
C
C
C     banded user-supplied matrix
400   LENPD=(2*IWM(LML)+IWM(LMU)+1)*NEQ
      DO 410 I=1,LENPD
410      WM(NPDM1+I)=0.0E0
      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      MEBAND=2*IWM(LML)+IWM(LMU)+1
      GO TO 550
C
C
C     banded finite-difference-generated matrix
500   MBAND=IWM(LML)+IWM(LMU)+1
      MBA=MIN0(MBAND,NEQ)
      MEBAND=MBAND+IWM(LML)
      MEB1=MEBAND-1
      MSAVE=(NEQ/MBAND)+1
      ISAVE=NTEMP-1
      IPSAVE=ISAVE+MSAVE
      IRES=0
      SQUR=SQRT(UROUND)
      DO 540 J=1,MBA
         DO 510 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          WM(ISAVE+K)=Y(N)
          WM(IPSAVE+K)=YPRIME(N)
          DEL=SQUR*AMAX1(ABS(Y(N)),ABS(H*YPRIME(N)),
     *      ABS(WT(N)))
          DEL=SIGN(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          Y(N)=Y(N)+DEL
510       YPRIME(N)=YPRIME(N)+CJ*DEL
      CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) RETURN
      DO 530 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          Y(N)=WM(ISAVE+K)
          YPRIME(N)=WM(IPSAVE+K)
          DEL=SQUR*AMAX1(ABS(Y(N)),ABS(H*YPRIME(N)),
     *      ABS(WT(N)))
          DEL=SIGN(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          DELINV=1.0E0/DEL
          I1=MAX0(1,(N-IWM(LMU)))
          I2=MIN0(NEQ,(N+IWM(LML)))
          II=N*MEB1-IWM(LML)+NPDM1
          DO 520 I=I1,I2
520         WM(II+I)=(E(I)-DELTA(I))*DELINV
530      CONTINUE
540   CONTINUE
C
C
C     do lu decomposition of banded pd
550   CALL SGBFA(WM(NPD),MEBAND,NEQ,
     *    IWM(LML),IWM(LMU),IWM(LIPVT),IER)
      RETURN
C------end of subroutine sdajac------
      END
