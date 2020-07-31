      SUBROUTINE SOS(FNC,NEQ,X,RTOLX,ATOLX,TOLF,IFLAG,RW,LRW,IW,LIW)
C***BEGIN PROLOGUE  SOS
C***DATE WRITTEN   801001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  F2A
C***KEYWORDS  BROWNS METHOD,NEWTONS METHOD,NONLINEAR EQUATIONS,ROOTS,
C             SOLUTIONS
C***AUTHOR  WATTS, H. A., (SNLA)
C***PURPOSE  SOS Solves a square system of nonlinear equations.
C***DESCRIPTION
C
C     SOS solves a system of NEQ simultaneous nonlinear equations in
C     NEQ unknowns.  That is, it solves the problem   F(X)=0
C     where X is a vector with components  X(1),...,X(NEQ)  and  F
C     is a vector of nonlinear functions.  Each equation is of the form
C
C               F (X(1),...,X(NEQ))=0     for K=1,...,NEQ.
C                K
C
C     The algorithm is based on an iterative method which is a
C     variation of Newton's method using Gaussian elimination
C     in a manner similar to the Gauss-Seidel process.  Convergence
C     is roughly quadratic.  All partial derivatives required by
C     the algorithm are approximated by first difference quotients.
C     The convergence behavior of this code is affected by the
C     ordering of the equations, and it is advantageous to place linear
C     and mildly nonlinear equations first in the ordering.
C
C     Actually, SOS is merely an interfacing routine for
C     calling subroutine SOSEQS which embodies the solution
C     algorithm.  The purpose of this is to add greater
C     flexibility and ease of use for the prospective user.
C
C     SOSEQS calls the accompanying routine SOSSOL, which solves special
C     triangular linear systems by back-substitution.
C
C     The user must supply a function subprogram which evaluates the
C     K-th equation only (K specified by SOSEQS) for each call
C     to the subprogram.
C
C     SOS represents an implementation of the mathematical algorithm
C     described in the references below.  It is a modification of the
C     code SOSNLE written by H. A. Watts in 1973.
C
C     References
C       1. K. M. Brown, Algorithm 316, Comm. A.C.M., Vol. 10, 1967,
C          pp. 728-729.
C       2. K. M. Brown, A Quadratically Convergent Newton-Like Method
C          Based Upon Gaussian Elimination, SIAM J. Num. Anal., Vol. 6,
C          1969, pp. 560-569.
C
C
C   -Input-
C
C     FNC -name of the function program which evaluates the equations.
C          This name must be in an EXTERNAL statement in the calling
C          program.  The user must supply FNC in the form FNC(X,K),
C          where X is the solution vector (which must be dimensioned
C          in FNC) and FNC returns the value of the K-th function.
C
C     NEQ -number of equations to be solved.
C
C     X   -solution vector.  Initial guesses must be supplied.
C
C     RTOLX -relative error tolerance used in the convergence criteria.
C          Each solution component X(I) is checked by an accuracy test
C          of the form   ABS(X(I)-XOLD(I)) .LE. RTOLX*ABS(X(I))+ATOLX,
C          where XOLD(I) represents the previous iteration value.
C          RTOLX must be non-negative.
C
C     ATOLX -absolute error tolerance used in the convergence criteria.
C          ATOLX must be non-negative.  If the user suspects some
C          solution component may be zero, he should set ATOLX to an
C          appropriate (depends on the scale of the remaining variables)
C          positive value for better efficiency.
C
C     TOLF -residual error tolerance used in the convergence criteria.
C          Convergence will be indicated if all residuals (values of the
C          functions or equations) are not bigger than TOLF in
C          magnitude.  Note that extreme care must be given in assigning
C          an appropriate value for TOLF because this convergence test
C          is dependent on the scaling of the equations.  An
C          inappropriate value can cause premature termination of the
C          iteration process.
C
C     IFLAG -optional input indicator.  You must set  IFLAG=-1  if you
C          want to use any of the optional input items listed below.
C          Otherwise set it to zero.
C
C     RW  -a real work array which is split apart by SOS and used
C          internally by SOSEQS.
C
C     LRW -dimension of the RW array.  LRW must be at least
C                    1 + 6*NEQ + NEQ*(NEQ+1)/2
C
C     IW  -an integer work array which is split apart by SOS and used
C          internally by SOSEQS.
C
C     LIW -dimension of the IW array. LIW must be at least  3 + NEQ.
C
C   -Optional Input-
C
C     IW(1) -internal printing parameter.  You must set  IW(1)=-1  if
C          you want the intermediate solution iterates to be printed.
C
C     IW(2) -iteration limit.  The maximum number of allowable
C          iterations can be specified, if desired.  To override the
C          default value of 50, set IW(2) to the number wanted.
C
C     Remember, if you tell the code that you are using one of the
C               options (by setting IFLAG=-1), you must supply values
C               for both IW(1) and IW(2).
C
C **********************************************************************
C
C   -Output-
C
C     X   -solution vector.
C
C     IFLAG -status indicator
C
C                         *** Convergence to a Solution ***
C
C          1 means satisfactory convergence to a solution was achieved.
C            Each solution component X(I) satisfies the error tolerance
C            test   ABS(X(I)-XOLD(I)) .LE. RTOLX*ABS(X(I))+ATOLX.
C
C          2 means procedure converged to a solution such that all
C            residuals are at most TOLF in magnitude,
C            ABS(FNC(X,I)) .LE. TOLF.
C
C          3 means that conditions for both IFLAG=1 and IFLAG=2 hold.
C
C          4 means possible numerical convergence.  Behavior indicates
C            limiting precision calculations as a result of user asking
C            for too much accuracy or else convergence is very slow.
C            Residual norms and solution increment norms have
C            remained roughly constant over several consecutive
C            iterations.
C
C                         *** Task Interrupted ***
C
C          5 means the allowable number of iterations has been met
C            without obtaining a solution to the specified accuracy.
C            Very slow convergence may be indicated.  Examine the
C            approximate solution returned and see if the error
C            tolerances seem appropriate.
C
C          6 means the allowable number of iterations has been met and
C            the iterative process does not appear to be converging.
C            A local minimum may have been encountered or there may be
C            limiting precision difficulties.
C
C          7 means that the iterative scheme appears to be diverging.
C            Residual norms and solution increment norms have
C            increased over several consecutive iterations.
C
C                         *** Task Cannot Be Continued ***
C
C          8 means that a Jacobian-related matrix was singular.
C
C          9 means improper input parameters.
C
C          *** IFLAG should be examined after each call to      ***
C          *** SOS with the appropriate action being taken.     ***
C
C
C     RW(1) -contains a norm of the residual.
C
C     IW(3) -contains the number of iterations used by the process.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  SOSEQS,XERRWV
C***END PROLOGUE  SOS
C
C
      EXTERNAL FNC
C
      DIMENSION X(NEQ), RW(LRW), IW(LIW)
C
C***FIRST EXECUTABLE STATEMENT  SOS
      INPFLG = IFLAG
C
C     CHECK FOR VALID INPUT
C
      IF (NEQ .GT. 0) GO TO 10
      CALL XERRWV( 'SOS -- THE NUMBER OF EQUATIONS MUST BE A POSITIVE IN
     1TEGER.                   YOU HAVE CALLED THE CODE WITH  NEQ = (I1)
     2.', 119, 1, 1, 1, NEQ, 0, 0, 0., 0.)
      IFLAG = 9
C
   10 IF (RTOLX .GE. 0. .AND. ATOLX .GE. 0.) GO TO 20
      CALL XERRWV(   'SOS -- THE ERROR TOLERANCES FOR THE SOLUTION ITERA
     1TES CANNOT BE              NEGATIVE. YOU HAVE CALLED THE CODE WITH
     2  RTOLX = (R1)   AND           ATOLX = (R2).', 160, 2, 1, 0, 0, 0,
     3  2, RTOLX, ATOLX)
      IFLAG = 9
C
   20 IF (TOLF .GE. 0.) GO TO 30
      CALL XERRWV( 'SOS -- THE RESIDUAL ERROR TOLERANCE MUST BE NON-NEGA
     1TIVE.                    YOU HAVE CALLED THE CODE WITH  TOLF = (R1
     2).', 120, 3, 1, 0, 0, 0, 1, TOLF, 0.)
      IFLAG = 9
C
   30 IPRINT = 0
      MXIT = 50
      IF (INPFLG.NE.(-1)) GO TO 40
      IF (IW(1) .EQ. (-1)) IPRINT = -1
      MXIT = IW(2)
      IF (MXIT .GT. 0) GO TO 40
      CALL XERRWV(   'SOS -- YOU HAVE TOLD THE CODE TO USE OPTIONAL INPU
     1T ITEMS BY SETTING         IFLAG=-1. HOWEVER YOU HAVE CALLED THE C
     2ODE WITH THE MAXIMUM            ALLOWABLE NUMBER OF ITERATIONS SET
     3 TO IW(2) = (I1).', 199, 4, 1, 1, MXIT, 0, 0, 0., 0.)
      IFLAG = 9
C
   40 NC = (NEQ*(NEQ+1))/2
      IF (LRW .GE. 1+6*NEQ+NC) GO TO 50
      CALL XERRWV(    'SOS -- DIMENSION OF THE RW ARRAY MUST BE AT LEAST
     1                                   1 + 6*NEQ + NEQ*(NEQ+1)/2     .
     2                                YOU HAVE CALLED THE CODE WITH  LRW
     3= (I1). '      , 189, 5, 1, 1, LRW, 0, 0, 0., 0.)
      IFLAG = 9
C
   50 IF (LIW .GE. 3+NEQ) GO TO 60
      CALL XERRWV(    'SOS -- DIMENSION OF THE IW ARRAY MUST BE AT LEAST
     1   3 + NEQ.                 YOU HAVE CALLED THE CODE WITH  LIW = (
     2I1).', 119, 6, 1, 1, LIW, 0, 0, 0., 0.)
      IFLAG = 9
C
   60 IF (IFLAG .EQ. 9) RETURN
C
      NCJS = 6
      NSRRC = 4
      NSRI = 5
C
      K1 = NC + 2
      K2 = K1 + NEQ
      K3 = K2 + NEQ
      K4 = K3 + NEQ
      K5 = K4 + NEQ
      K6 = K5 + NEQ
C
      CALL SOSEQS(FNC, NEQ, X, RTOLX, ATOLX, TOLF, IFLAG, MXIT, NCJS,
     1            NSRRC, NSRI, IPRINT, RW(1), RW(2), NC, RW(K1),
     2            RW(K2), RW(K3), RW(K4), RW(K5), RW(K6), IW(4))
C
      IW(3) = MXIT
      RETURN
      END
