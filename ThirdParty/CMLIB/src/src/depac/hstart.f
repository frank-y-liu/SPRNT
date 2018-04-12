      SUBROUTINE HSTART(F,NEQ,A,B,Y,YPRIME,ETOL,MORDER,SMALL,BIG,SPY,
     1   PV,YP,SF,RPAR,IPAR,H)
C***BEGIN PROLOGUE  HSTART
C***REFER TO  DEABM,DEBDF,DERKF
C
C   HSTART computes a starting step size to be used in solving initial
C   value problems in ordinary differential equations.
C **********************************************************************
C  Abstract
C
C     Subroutine HSTART computes a starting step size to be used by an
C     initial value method in solving ordinary differential equations.
C     It is based on an estimate of the local Lipschitz constant for the
C     differential equation   (lower bound on a norm of the Jacobian) ,
C     a bound on the differential equation  (first derivative) , and
C     a bound on the partial derivative of the equation with respect to
C     the independent variable.
C     (all approximated near the initial point A)
C
C     Subroutine HSTART uses a function subprogram VNORM for computing
C     a vector norm. The maximum norm is presently utilized though it
C     can easily be replaced by any other vector norm. It is presumed
C     that any replacement norm routine would be carefully coded to
C     prevent unnecessary underflows or overflows from occurring, and
C     also, would not alter the vector or number of components.
C
C **********************************************************************
C  On input you must provide the following
C
C      F -- This is a subroutine of the form
C                               F(X,U,UPRIME,RPAR,IPAR)
C             which defines the system of first order differential
C             equations to be solved. For the given values of X and the
C             vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must
C             evaluate the NEQ components of the system of differential
C             equations  dU/DX=F(X,U)  and store the derivatives in the
C             array UPRIME(*), that is,  UPRIME(I) = * dU(I)/DX *  for
C             equations I=1,...,NEQ.
C
C             Subroutine F must NOT alter X or U(*). You must declare
C             the name F in an external statement in your program that
C             calls HSTART. You must dimension U and UPRIME in F.
C
C             RPAR and IPAR are real and integer parameter arrays which
C             you can use for communication between your program and
C             subroutine F. They are not used or altered by HSTART. If
C             you do not need RPAR or IPAR, ignore these parameters by
C             treating them as dummy arguments. If you do choose to use
C             them, dimension them in your program and in F as arrays
C             of appropriate length.
C
C      NEQ -- This is the number of (first order) differential equations
C             to be integrated.
C
C      A -- This is the initial point of integration.
C
C      B -- This is a value of the independent variable used to define
C             the direction of integration. A reasonable choice is to
C             set  B  to the first point at which a solution is desired.
C             you can also use  B, if necessary, to restrict the length
C             of the first integration step because the algorithm will
C             not compute a starting step length which is bigger than
C             ABS(B-A), unless  B  has been chosen too close to  A.
C             (it is presumed that HSTART has been called with  B
C             different from  A  on the machine being used. Also see the
C             discussion about the parameter  SMALL.)
C
C      Y(*) -- This is the vector of initial values of the NEQ solution
C             components at the initial point  A.
C
C      YPRIME(*) -- This is the vector of derivatives of the NEQ
C             solution components at the initial point  A.
C             (defined by the differential equations in subroutine F)
C
C      ETOL -- This is the vector of error tolerances corresponding to
C             the NEQ solution components. It is assumed that all
C             elements are positive. Following the first integration
C             step, the tolerances are expected to be used by the
C             integrator in an error test which roughly requires that
C                        ABS(local error) .LE. ETOL
C             for each vector component.
C
C      MORDER -- This is the order of the formula which will be used by
C             the initial value method for taking the first integration
C             step.
C
C      SMALL -- This is a small positive machine dependent constant
C             which is used for protecting against computations with
C             numbers which are too small relative to the precision of
C             floating point arithmetic.  SMALL  should be set to
C             (approximately) the smallest positive real number such
C             that  (1.+SMALL) .GT. 1.  on the machine being used. The
C             quantity  SMALL**(3/8)  is used in computing increments of
C             variables for approximating derivatives by differences.
C             Also the algorithm will not compute a starting step length
C             which is smaller than  100*SMALL*ABS(A).
C
C      BIG -- This is a large positive machine dependent constant which
C             is used for preventing machine overflows. A reasonable
C             choice is to set big to (approximately) the square root of
C             the largest real number which can be held in the machine.
C
C      SPY(*),PV(*),YP(*),SF(*) -- These are real work arrays of length
C             NEQ which provide the routine with needed storage space.
C
C      RPAR,IPAR -- These are parameter arrays, of real and integer
C             type, respectively, which can be used for communication
C             between your program and the F subroutine. They are not
C             used or altered by HSTART.
C
C **********************************************************************
C  On Output  (after the return from HSTART),
C
C      H -- Is an appropriate starting step size to be attempted by the
C             differential equation method.
C
C           All parameters in the call list remain unchanged except for
C           the working arrays SPY(*),PV(*),YP(*), and SF(*).
C
C **********************************************************************
C***ROUTINES CALLED  VNORM
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***END PROLOGUE  HSTART
C
C
      DIMENSION Y(NEQ),YPRIME(NEQ),ETOL(NEQ),
     1          SPY(NEQ),PV(NEQ),YP(NEQ),SF(NEQ)       ,RPAR(*),IPAR(*)
      EXTERNAL F
C
C.......................................................................
C
C***FIRST EXECUTABLE STATEMENT  HSTART
      DX=B-A
      ABSDX=ABS(DX)
      RELPER=SMALL**0.375
C
C.......................................................................
C
C     COMPUTE AN APPROXIMATE BOUND (DFDXB) ON THE PARTIAL
C     DERIVATIVE OF THE EQUATION WITH RESPECT TO THE
C     INDEPENDENT VARIABLE. PROTECT AGAINST AN OVERFLOW.
C     ALSO COMPUTE A BOUND (FBND) ON THE FIRST DERIVATIVE LOCALLY.
C
      DA=SIGN(AMAX1(AMIN1(RELPER*ABS(A),ABSDX),100.*SMALL*ABS(A)),DX)
      IF (DA .EQ. 0.) DA=RELPER*DX
      CALL F(A+DA,Y,SF,RPAR,IPAR)
      DO 10 J=1,NEQ
   10   YP(J)=SF(J)-YPRIME(J)
      DELF=VNORM(YP,NEQ)
      DFDXB=BIG
      IF (DELF .LT. BIG*ABS(DA)) DFDXB=DELF/ABS(DA)
      FBND=VNORM(SF,NEQ)
C
C.......................................................................
C
C     COMPUTE AN ESTIMATE (DFDUB) OF THE LOCAL LIPSCHITZ CONSTANT FOR
C     THE SYSTEM OF DIFFERENTIAL EQUATIONS. THIS ALSO REPRESENTS AN
C     ESTIMATE OF THE NORM OF THE JACOBIAN LOCALLY.
C     THREE ITERATIONS (TWO WHEN NEQ=1) ARE USED TO ESTIMATE THE
C     LIPSCHITZ CONSTANT BY NUMERICAL DIFFERENCES. THE FIRST
C     PERTURBATION VECTOR IS BASED ON THE INITIAL DERIVATIVES AND
C     DIRECTION OF INTEGRATION. THE SECOND PERTURBATION VECTOR IS
C     FORMED USING ANOTHER EVALUATION OF THE DIFFERENTIAL EQUATION.
C     THE THIRD PERTURBATION VECTOR IS FORMED USING PERTURBATIONS BASED
C     ONLY ON THE INITIAL VALUES. COMPONENTS THAT ARE ZERO ARE ALWAYS
C     CHANGED TO NON-ZERO VALUES (EXCEPT ON THE FIRST ITERATION). WHEN
C     INFORMATION IS AVAILABLE, CARE IS TAKEN TO ENSURE THAT COMPONENTS
C     OF THE PERTURBATION VECTOR HAVE SIGNS WHICH ARE CONSISTENT WITH
C     THE SLOPES OF LOCAL SOLUTION CURVES.
C     ALSO CHOOSE THE LARGEST BOUND (FBND) FOR THE FIRST DERIVATIVE.
C
C                       PERTURBATION VECTOR SIZE IS HELD CONSTANT FOR
C                       ALL ITERATIONS. COMPUTE THIS CHANGE FROM THE
C                               SIZE OF THE VECTOR OF INITIAL VALUES.
      DELY=RELPER*VNORM(Y,NEQ)
      IF (DELY .EQ. 0.) DELY=RELPER
      DELY=SIGN(DELY,DX)
      DELF=VNORM(YPRIME,NEQ)
      FBND=AMAX1(FBND,DELF)
      IF (DELF .EQ. 0.) GO TO 30
C                       USE INITIAL DERIVATIVES FOR FIRST PERTURBATION
      DO 20 J=1,NEQ
        SPY(J)=YPRIME(J)
   20   YP(J)=YPRIME(J)
      GO TO 50
C                       CANNOT HAVE A NULL PERTURBATION VECTOR
   30 DO 40 J=1,NEQ
        SPY(J)=0.
   40   YP(J)=1.
      DELF=VNORM(YP,NEQ)
C
   50 DFDUB=0.
      LK=MIN0(NEQ+1,3)
      DO 140 K=1,LK
C                       DEFINE PERTURBED VECTOR OF INITIAL VALUES
        DO 60 J=1,NEQ
   60     PV(J)=Y(J)+DELY*(YP(J)/DELF)
        IF (K .EQ. 2) GO TO 80
C                       EVALUATE DERIVATIVES ASSOCIATED WITH PERTURBED
C                       VECTOR  AND  COMPUTE CORRESPONDING DIFFERENCES
        CALL F(A,PV,YP,RPAR,IPAR)
        DO 70 J=1,NEQ
   70     PV(J)=YP(J)-YPRIME(J)
        GO TO 100
C                       USE A SHIFTED VALUE OF THE INDEPENDENT VARIABLE
C                                             IN COMPUTING ONE ESTIMATE
   80   CALL F(A+DA,PV,YP,RPAR,IPAR)
        DO 90 J=1,NEQ
   90     PV(J)=YP(J)-SF(J)
C                       CHOOSE LARGEST BOUNDS ON THE FIRST DERIVATIVE
C                                      AND A LOCAL LIPSCHITZ CONSTANT
  100   FBND=AMAX1(FBND,VNORM(YP,NEQ))
        DELF=VNORM(PV,NEQ)
        IF (DELF .GE. BIG*ABS(DELY)) GO TO 150
        DFDUB=AMAX1(DFDUB,DELF/ABS(DELY))
        IF (K .EQ. LK) GO TO 160
C                       CHOOSE NEXT PERTURBATION VECTOR
        IF (DELF .EQ. 0.) DELF=1.
        DO 130 J=1,NEQ
          IF (K .EQ. 2) GO TO 110
          DY=ABS(PV(J))
          IF (DY .EQ. 0.) DY=DELF
          GO TO 120
  110     DY=Y(J)
          IF (DY .EQ. 0.) DY=DELY/RELPER
  120     IF (SPY(J) .EQ. 0.) SPY(J)=YP(J)
          IF (SPY(J) .NE. 0.) DY=SIGN(DY,SPY(J))
  130     YP(J)=DY
  140   DELF=VNORM(YP,NEQ)
C
C                       PROTECT AGAINST AN OVERFLOW
  150 DFDUB=BIG
C
C.......................................................................
C
C     COMPUTE A BOUND (YDPB) ON THE NORM OF THE SECOND DERIVATIVE
C
  160 YDPB=DFDXB+DFDUB*FBND
C
C.......................................................................
C
C     DEFINE THE TOLERANCE PARAMETER UPON WHICH THE STARTING STEP SIZE
C     IS TO BE BASED.  A VALUE IN THE MIDDLE OF THE ERROR TOLERANCE
C     RANGE IS SELECTED.
C
      TOLMIN=BIG
      TOLSUM=0.
      DO 170 K=1,NEQ
        TOLEXP=ALOG10(ETOL(K))
        TOLMIN=AMIN1(TOLMIN,TOLEXP)
 170    TOLSUM=TOLSUM+TOLEXP
      TOLP=10.**(0.5*(TOLSUM/FLOAT(NEQ)+TOLMIN)/FLOAT(MORDER+1))
C
C.......................................................................
C
C     COMPUTE A STARTING STEP SIZE BASED ON THE ABOVE FIRST AND SECOND
C     DERIVATIVE INFORMATION
C
C                       RESTRICT THE STEP LENGTH TO BE NOT BIGGER THAN
C                       ABS(B-A).   (UNLESS  B  IS TOO CLOSE TO  A)
      H=ABSDX
C
      IF (YDPB .NE. 0.  .OR.  FBND .NE. 0.) GO TO 180
C
C                       BOTH FIRST DERIVATIVE TERM (FBND) AND SECOND
C                                    DERIVATIVE TERM (YDPB) ARE ZERO
      IF (TOLP .LT. 1.) H=ABSDX*TOLP
      GO TO 200
C
  180 IF (YDPB .NE. 0.) GO TO 190
C
C                       ONLY SECOND DERIVATIVE TERM (YDPB) IS ZERO
      IF (TOLP .LT. FBND*ABSDX) H=TOLP/FBND
      GO TO 200
C
C                       SECOND DERIVATIVE TERM (YDPB) IS NON-ZERO
  190 SRYDPB=SQRT(0.5*YDPB)
      IF (TOLP .LT. SRYDPB*ABSDX) H=TOLP/SRYDPB
C
C                       FURTHER RESTRICT THE STEP LENGTH TO BE NOT
C                                                 BIGGER THAN  1/DFDUB
  200 IF (H*DFDUB .GT. 1.) H=1./DFDUB
C
C                       FINALLY, RESTRICT THE STEP LENGTH TO BE NOT
C                       SMALLER THAN  100*SMALL*ABS(A).  HOWEVER, IF
C                       A=0. AND THE COMPUTED H UNDERFLOWED TO ZERO,
C                       THE ALGORITHM RETURNS  SMALL*ABS(B)  FOR THE
C                                                       STEP LENGTH.
      H=AMAX1(H,100.*SMALL*ABS(A))
      IF (H .EQ. 0.) H=SMALL*ABS(B)
C
C                       NOW SET DIRECTION OF INTEGRATION
      H=SIGN(H,DX)
C
      RETURN
      END
