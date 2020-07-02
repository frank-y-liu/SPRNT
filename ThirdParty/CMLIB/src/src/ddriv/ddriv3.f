      SUBROUTINE DDRIV3 (N,T,Y,F,NSTATE,TOUT,NTASK,NROOT,EPS,EWT,IERROR,
     8   MINT,MITER,IMPL,ML,MU,MXORD,HMAX,WORK,LENW,IWORK,LENIW,JACOBN,
     8   FA,NDE,MXSTEP,G)
C***BEGIN PROLOGUE  DDRIV3
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850924   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD,
C             DOUBLE PRECISION
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***PURPOSE  The function of DDRIV3 is to solve N ordinary differential
C            equations of the form dY(I)/dT = F(Y(I),T), given the
C            initial conditions Y(I) = YI.  The program has options to
C            allow the solution of both stiff and non-stiff differential
C            equations.  Other important options are available.  DDRIV3
C            uses double precision arithmetic.
C***DESCRIPTION
C
C  I.  ABSTRACT  .......................................................
C
C    The primary function of DDRIV3 is to solve N ordinary differential
C    equations of the form dY(I)/dT = F(Y(I),T), given the initial
C    conditions Y(I) = YI.  The program has options to allow the
C    solution of both stiff and non-stiff differential equations.  In
C    addition, DDRIV3 may be used to solve:
C      1. The initial value problem, A*dY(I)/dT = F(Y(I),T), where A is
C         a non-singular matrix depending on Y and T.
C      2. The hybrid differential/algebraic initial value problem,
C         A*dY(I)/dT = F(Y(I),T), where A is a vector (whose values may
C         depend upon Y and T) some of whose components will be zero
C         corresponding to those equations which are algebraic rather
C         than differential.
C    DDRIV3 is to be called once for each output point of T.
C
C  II.  PARAMETERS  ....................................................
C
C       (REMEMBER--To run DDRIV3 correctly in double precision, ALL
C       non-integer arguments in the call sequence, including
C       arrays, MUST be declared double precision.)
C    The user should use parameter names in the call sequence of DDRIV3
C    for those quantities whose value may be altered by DDRIV3.  The
C    parameters in the call sequence are:
C
C    N      = (Input) The number of dependent functions whose solution
C             is desired.  N must not be altered during a problem.
C
C    T      = The independent variable.  On input for the first call, T
C             is the initial point.  On output, T is the point at which
C             the solution is given.
C
C    Y      = The vector of dependent variables.  Y is used as input on
C             the first call, to set the initial values.  On output, Y
C             is the computed solution vector.  This array Y is passed
C             in the call sequence of the user-provided routines F,
C             JACOBN, FA, USERS, and G.
C
C    F      = A subroutine supplied by the user.  The name must be
C             declared EXTERNAL in the user's calling program.  This
C             subroutine is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   DOUBLE PRECISION Y(*), YDOT(*)
C                     .
C                     .
C                   YDOT(1) = ...
C                     .
C                     .
C                   YDOT(N) = ...
C                   END (Sample)
C             This computes YDOT = F(Y,T), the right hand side of the
C             differential equations.  Here Y is a vector of length at
C             least N.  The actual length of Y is determined by the
C             user's declaration in the program which calls DDRIV3.
C             Thus the dimensioning of Y in F, while required by FORTRAN
C             convention, does not actually allocate any storage.  When
C             this subroutine is called, the first N components of Y are
C             intermediate approximations to the solution components.
C             The user should not alter these values.  Here YDOT is a
C             vector of length N.  The user should only compute YDOT(I)
C             for I from 1 to N.
C
C    NSTATE = An integer describing the status of integration.  The
C             meaning of NSTATE is as follows:
C               1  (Input) Means the first call to the routine.  This
C                  value must be set by the user.  On all subsequent
C                  calls the value of NSTATE should be tested by the
C                  user, but must not be altered.  (As a convenience to
C                  the user who may wish to put out the initial
C                  conditions, DDRIV3 can be called with NSTATE=1, and
C                  TOUT=T.  In this case the program will return with
C                  NSTATE unchanged, i.e., NSTATE=1.)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  MXSTEP steps without reaching TOUT.  The user can
C                  continue the integration by simply calling DDRIV3
C                  again.
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling DDRIV3
C                  again.
C               5  (Output) A root was found at a point less than TOUT.
C                  The user can continue the integration toward TOUT by
C                  simply calling DDRIV3 again.
C
C    TOUT   = (Input) The point at which the solution is desired.  The
C             position of TOUT relative to T on the first call
C             determines the direction of integration.
C
C    NTASK  = (Input) An index specifying the manner of returning the
C             solution, according to the following:
C               NTASK = 1  Means DDRIV3 will integrate past TOUT and
C                          interpolate the solution.  This is the most
C                          efficient mode.
C               NTASK = 2  Means DDRIV3 will return the solution after
C                          each internal integration step, or at TOUT,
C                          whichever comes first.  In the latter case,
C                          the program integrates exactly to TOUT.
C               NTASK = 3  Means DDRIV3 will adjust its internal step to
C                          reach TOUT exactly (useful if a singularity
C                          exists beyond TOUT.)
C
C    NROOT  = (Input) The number of equations whose roots are desired.
C             If NROOT is zero, the root search is not active.  This
C             option is useful for obtaining output at points which are
C             not known in advance, but depend upon the solution, e.g.,
C             when some solution component takes on a specified value.
C             The root search is carried out using the user-written
C             function G (see description of G below.)  DDRIV3 attempts
C             to find the value of T at which one of the equations
C             changes sign.  DDRIV3 can find at most one root per
C             equation per internal integration step, and will then
C             return the solution either at TOUT or at a root, whichever
C             occurs first in the direction of integration.  The index
C             of the equation whose root is being reported is stored in
C             the sixth element of IWORK.
C             NOTE: NROOT is never altered by this program.
C
C    EPS    = On input, the requested relative accuracy in all solution
C             components.  EPS = 0 is allowed.  On output, the adjusted
C             relative accuracy if the input value was too small.  The
C             value of EPS should be set as large as is reasonable,
C             because the amount of work done by DDRIV3 increases as EPS
C             decreases.
C
C    EWT    = (Input) Problem zero, i.e., the smallest, nonzero,
C             physically meaningful value for the solution.  (Array,
C             possibly of length one.  See following description of
C             IERROR.)  Setting EWT smaller than necessary can adversely
C             affect the running time.
C
C    IERROR = (Input) Error control indicator.  A value of 3 is
C             suggested for most problems.  Other choices and detailed
C             explanations of EWT and IERROR are given below for those
C             who may need extra flexibility.
C
C             These last three input quantities EPS, EWT and IERROR
C             control the accuracy of the computed solution.  EWT and
C             IERROR are used internally to compute an array YWT.  One
C             step error estimates divided by YWT(I) are kept less than
C             EPS in root mean square norm.
C                 IERROR (Set by the user) =
C                 1  Means YWT(I) = 1. (Absolute error control)
C                                   EWT is ignored.
C                 2  Means YWT(I) = ABS(Y(I)),  (Relative error control)
C                                   EWT is ignored.
C                 3  Means YWT(I) = MAX(ABS(Y(I)), EWT(1)).
C                 4  Means YWT(I) = MAX(ABS(Y(I)), EWT(I)).
C                    This choice is useful when the solution components
C                    have differing scales.
C                 5  Means YWT(I) = EWT(I).
C             If IERROR is 3, EWT need only be dimensioned one.
C             If IERROR is 4 or 5, the user must dimension EWT at least
C             N, and set its values.
C
C    MINT   = (Input) The integration method indicator.
C               MINT = 1  Means the Adams methods, and is used for
C                         non-stiff problems.
C               MINT = 2  Means the stiff methods of Gear (i.e., the
C                         backward differentiation formulas), and is
C                         used for stiff problems.
C               MINT = 3  Means the program dynamically selects the
C                         Adams methods when the problem is non-stiff
C                         and the Gear methods when the problem is
C                         stiff.  When using the Adams methods, the
C                         program uses a value of MITER=0; when using
C                         the Gear methods, the program uses the value
C                         of MITER provided by the user.  Only a value
C                         of IMPL = 0 and a value of MITER = 1, 2, 4, or
C                         5 is allowed for this option.  The user may
C                         not alter the value of MINT or MITER without
C                         restarting, i.e., setting NSTATE to 1.
C
C    MITER  = (Input) The iteration method indicator.
C               MITER = 0  Means functional iteration.  This value is
C                          suggested for non-stiff problems.
C               MITER = 1  Means chord method with analytic Jacobian.
C                          In this case, the user supplies subroutine
C                          JACOBN (see description below).
C               MITER = 2  Means chord method with Jacobian calculated
C                          internally by finite differences.
C               MITER = 3  Means chord method with corrections computed
C                          by the user-written routine named USERS.
C                          This option allows all matrix algebra and
C                          storage decisions to be made by the user.
C                          The routine USERS is called by DDRIV3 when
C                          certain linear systems must be solved.  The
C                          user may choose any method to form, store and
C                          solve these systems in order to obtain the
C                          solution result that is returned to DDRIV3.
C                          In particular, this allows sparse matrix
C                          methods to be used.
C                          The call sequence for this routine is
C
C                           SUBROUTINE USERS (Y, YH, YWT, SAVE1, SAVE2,
C                          8              T, H, EL, IMPL, N, NDE, IFLAG)
C                           DOUBLE PRECISION Y(*), YH(*), YWT(*),
C                          8        SAVE1(*), SAVE2(*), T, H, EL
C
C                          The input variable IFLAG indicates what
C                          action is to be taken.  Subroutine USERS
C                          should perform the following operations,
C                          depending on the value of IFLAG and IMPL.
C
C                          IFLAG = 0
C                            IMPL = 0.  USERS is not called.
C                            IMPL = 1 or 2.  Solve the system
C                                A*X = SAVE2,
C                              returning the result in SAVE2.  The array
C                              SAVE1 can be used as a work array.
C
C                          IFLAG = 1
C                            IMPL = 0.  Compute, decompose and store the
C                              matrix (I - H*EL*J), where I is the
C                              identity matrix and J is the Jacobian
C                              matrix of the right hand side.  The array
C                              SAVE1 can be used as a work array.
C                            IMPL = 1 or 2. Compute, decompose and store
C                              the matrix (A - H*EL*J).  The array SAVE1
C                              can be used as a work array.
C
C                          IFLAG = 2
C                            IMPL = 0.   Solve the system
C                                (I - H*EL*J)*X = H*SAVE2 - YH - SAVE1,
C                              returning the result in SAVE2.
C                            IMPL = 1, or 2.  Solve the system
C                              (A - H*EL*J)*X = H*SAVE2 - A*(YH + SAVE1)
C                              returning the result in SAVE2.
C                            The array SAVE1 should not be altered.
C
C                          When using a value of MITER = 3, the
C                          subroutine FA is not required, even if IMPL
C                          is not 0.  For further information on using
C                          this option, see section IV-F below.
C
C               MITER = 4  Means the same as MITER = 1 but the A and
C                          Jacobian matrices are assumed to be banded.
C               MITER = 5  Means the same as MITER = 2 but the A and
C                          Jacobian matrices are assumed to be banded.
C
C    IMPL   = (Input) The implicit method indicator.
C               IMPL = 0 Means solving dY(I)/dT = F(Y(I),T).
C               IMPL = 1 Means solving A*dY(I)/dT = F(Y(I),T),
C                        non-singular A (see description of FA below.)
C                        Only MINT = 1 or 2, and MITER = 1, 2, 3, 4, or
C                        5 are allowed for this option.
C               IMPL = 2 Means solving certain systems of hybrid
C                        differential/algebraic equations (see
C                        description of FA below.)  Only MINT = 2 and
C                        MITER = 1, 2, 3, 4, or 5, are allowed for this
C                        option.
C               The value of IMPL must not be changed during a problem.
C
C    ML     = (Input) The lower half-bandwidth in the case of a banded
C             A or Jacobian matrix.  (I.e., maximum(R-C) for nonzero
C             A(R,C).)
C
C    MU     = (Input) The upper half-bandwidth in the case of a banded
C             A or Jacobian matrix.  (I.e., maximum(C-R).)
C
C    MXORD  = (Input) The maximum order desired. This is .LE. 12 for
C             the Adams methods and .LE. 5 for the Gear methods.  Normal
C             value is 12 and 5, respectively.  If MINT is 3, the
C             maximum order used will be MIN(MXORD, 12) when using the
C             Adams methods, and MIN(MXORD, 5) when using the Gear
C             methods.  MXORD must not be altered during a problem.
C
C    HMAX   = (Input) The maximum magnitude of the step size that will
C             be used for the problem.  This is useful for ensuring that
C             important details are not missed.  If this is not the
C             case, a large value, such as the interval length, is
C             suggested.
C
C    WORK
C    LENW   = (Input)
C             WORK is an array of LENW double precision words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       DOUBLE PRECISION WORK(...)
C             The following table gives the required minimum value for
C             the length of WORK, depending on the value of IMPL and
C             MITER.  LENW should be set to the value used.  The
C             contents of WORK should not be disturbed between calls to
C             DDRIV3.
C
C      IMPL =   0                   1                   2
C              ---------------------------------------------------------
C MITER =  0   (MXORD+4)*N +       Not allowed         Not allowed
C              2*NROOT + 204
C
C         1,2  N*N+(MXORD+4)*N     2*N*N+(MXORD+4)*N   N*N+(MXORD+5)*N
C              + 2*NROOT + 204     + 2*NROOT + 204     + 2*NROOT + 204
C
C          3   (MXORD+4)*N +       (MXORD+4)*N +       (MXORD+4)*N +
C              2*NROOT + 204       2*NROOT + 204       2*NROOT + 204
C
C         4,5  (2*ML+MU)*N +       (4*ML+2*MU)*N +     (2*ML+MU)*N +
C              (MXORD+5)*N +       (MXORD+6)*N +       (MXORD+6)*N +
C              2*NROOT + 204       2*NROOT + 204       2*NROOT + 204
C              ---------------------------------------------------------
C
C    IWORK
C    LENIW  = (Input)
C             IWORK is an integer array of length LENIW used internally
C             for temporary storage.  The user must allocate space for
C             this array in the calling program by a statement such as
C                       INTEGER IWORK(...)
C             The length of IWORK should be at least
C               21      if MITER is 0 or 3, or
C               N+21    if MITER is 1, 2, 4, or 5, or MINT is 3,
C             and LENIW should be set to the value used.  The contents
C             of IWORK should not be disturbed between calls to DDRIV3.
C
C    JACOBN = A subroutine supplied by the user, if MITER is 1 or 4.
C             If this is the case, the name must be declared EXTERNAL in
C             the user's calling program.  Given a system of N
C             differential equations, it is meaningful to speak about
C             the partial derivative of the I-th right hand side with
C             respect to the J-th dependent variable.  In general there
C             are N*N such quantities.  Often however the equations can
C             be ordered so that the I-th differential equation only
C             involves dependent variables with index near I, e.g., I+1,
C             I-2.  Such a system is called banded.  If, for all I, the
C             I-th equation depends on at most the variables
C               Y(I-ML), Y(I-ML+1), ... , Y(I), Y(I+1), ... , Y(I+MU)
C             then we call ML+MU+1 the bandwith of the system.  In a
C             banded system many of the partial derivatives above are
C             automatically zero.  For the cases MITER = 1, 2, 4, and 5,
C             some of these partials are needed.  For the cases
C             MITER = 2 and 5 the necessary derivatives are
C             approximated numerically by DDRIV3, and we only ask the
C             user to tell DDRIV3 the value of ML and MU if the system
C             is banded.  For the cases MITER = 1 and 4 the user must
C             derive these partials algebraically and encode them in
C             subroutine JACOBN.  By computing these derivatives the
C             user can often save 20-30 per cent of the computing time.
C             Usually, however, the accuracy is not much affected and
C             most users will probably forego this option.  The optional
C             user-written subroutine JACOBN has the form:
C                   SUBROUTINE JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
C                   DOUBLE PRECISION Y(*), DFDY(MATDIM,*)
C                     .
C                     .
C                     Calculate values of DFDY
C                     .
C                     .
C                   END (Sample)
C             Here Y is a vector of length at least N.  The actual
C             length of Y is determined by the user's declaration in the
C             program which calls DDRIV3.  Thus the dimensioning of Y in
C             JACOBN, while required by FORTRAN convention, does not
C             actually allocate any storage.  When this subroutine is
C             called, the first N components of Y are intermediate
C             approximations to the solution components.  The user
C             should not alter these values.  If the system is not
C             banded (MITER=1), the partials of the I-th equation with
C             respect to the J-th dependent function are to be stored in
C             DFDY(I,J).  Thus partials of the I-th equation are stored
C             in the I-th row of DFDY.  If the system is banded
C             (MITER=4), then the partials of the I-th equation with
C             respect to Y(J) are to be stored in DFDY(K,J), where
C             K=I-J+MU+1.
C
C    FA     = A subroutine supplied by the user if IMPL is 1 or 2, and
C             MITER is not 3.  If so, the name must be declared EXTERNAL
C             in the user's calling program.  This subroutine computes
C             the array A, where A*dY(I)/dT = F(Y(I),T).
C             There are two cases:
C
C               IMPL=1.
C               Subroutine FA is of the form:
C                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
C                   DOUBLE PRECISION Y(*), A(MATDIM,*)
C                     .
C                     .
C                     Calculate ALL values of A
C                     .
C                     .
C                   END (Sample)
C               In this case A is assumed to be a nonsingular matrix,
C               with the same structure as DFDY (see JACOBN description
C               above).  Programming considerations prevent complete
C               generality.  If MITER is 1 or 2, A is assumed to be full
C               and the user must compute and store all values of
C               A(I,J), I,J=1, ... ,N.  If MITER is 4 or 5, A is assumed
C               to be banded with lower and upper half bandwidth ML and
C               MU.  The left hand side of the I-th equation is a linear
C               combination of dY(I-ML)/dT, dY(I-ML+1)/dT, ... ,
C               dY(I)/dT, ... , dY(I+MU-1)/dT, dY(I+MU)/dT.  Thus in the
C               I-th equation, the coefficient of dY(J)/dT is to be
C               stored in A(K,J), where K=I-J+MU+1.
C               NOTE: The array A will be altered between calls to FA.
C
C               IMPL=2.
C               Subroutine FA is of the form:
C                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
C                   DOUBLE PRECISION Y(*), A(*)
C                     .
C                     .
C                     Calculate non-zero values of A(1),...,A(NDE)
C                     .
C                     .
C                   END (Sample)
C               In this case it is assumed that the system is ordered by
C               the user so that the differential equations appear
C               first, and the algebraic equations appear last.  The
C               algebraic equations must be written in the form:
C               0 = F(Y(I),T).  When using this option it is up to the
C               user to provide initial values for the Y(I) that satisfy
C               the algebraic equations as well as possible.  It is
C               further assumed that A is a vector of length NDE.  All
C               of the components of A, which may depend on T, Y(I),
C               etc., must be set by the user to non-zero values.
C             Here Y is a vector of length at least N.  The actual
C             length of Y is determined by the user's declaration in the
C             program which calls DDRIV3.  Thus the dimensioning of Y in
C             FA, while required by FORTRAN convention, does not
C             actually allocate any storage.  When this subroutine is
C             called, the first N components of Y are intermediate
C             approximations to the solution components.  The user
C             should not alter these values.  FA is always called
C             immediately after calling F, with the same values of T
C             and Y.
C
C    NDE    = (Input) The number of differential equations.  This is
C             required only for IMPL = 2, with NDE .LT. N.
C
C    MXSTEP = (Input) The maximum number of internal steps allowed on
C             one call to DDRIV3.
C
C    G      = A double precision FORTRAN function supplied by the user
C             if NROOT is not 0.  In this case, the name must be
C             declared EXTERNAL in the user's calling program.  G is
C             repeatedly called with different values of IROOT to obtain
C             the value of each of the NROOT equations for which a root
C             is desired.  G is of the form:
C                   DOUBLE PRECISION FUNCTION G (N, T, Y, IROOT)
C                   DOUBLE PRECISION Y(*)
C                   GO TO (10, ...), IROOT
C              10   G = ...
C                     .
C                     .
C                   END (Sample)
C             Here, Y is a vector of length at least N, whose first N
C             components are the solution components at the point T.
C             The user should not alter these values.  The actual length
C             of Y is determined by the user's declaration in the
C             program which calls DDRIV3.  Thus the dimensioning of Y in
C             G, while required by FORTRAN convention, does not actually
C             allocate any storage.
C
C***LONG DESCRIPTION
C
C  III.  OTHER COMMUNICATION TO THE USER  ..............................
C
C    A. The solver communicates to the user through the parameters
C       above.  In addition it writes diagnostic messages through the
C       standard error handling program XERROR.  That program will
C       terminate the user's run if it detects a probable problem setup
C       error, e.g., insufficient storage allocated by the user for the
C       WORK array.  Messages are written on the standard error message
C       file.  At installations which have this error handling package
C       the user should determine the standard error handling file from
C       the local documentation.  Otherwise the short but serviceable
C       routine, XERROR, available with this package, can be used.  That
C       program writes on logical unit 6 to transmit messages.  A
C       complete description of XERROR is given in the Sandia
C       Laboratories report SAND78-1189 by R. E. Jones.  Following is a
C       list of possible errors.  Unless otherwise noted, all messages
C       come from DDRIV3:
C
C        No.  Type         Message
C        ---  ----         -------
C         1   Fatal        From DDRIV2: The integration method flag has
C                          an illegal value.
C         2   Warning      The output point is inconsistent with the
C                          value of NTASK and T.
C         3   Warning      Number of steps to reach TOUT exceeds MXSTEP.
C         4   Recoverable  Requested accuracy is too stringent.
C         5   Warning      Step size is below the roundoff level.
C         6   Fatal        EPS is less than zero.
C         7   Fatal        N is not positive.
C         8   Fatal        Insufficient work space provided.
C         9   Fatal        Improper value for MINT, MITER and/or IMPL.
C        10   Fatal        The IWORK array is too small.
C        11   Fatal        The step size has gone to zero.
C        12   Fatal        Excessive amount of work.
C        13   Fatal        For IMPL=1 or 2, the matrix A is singular.
C        14   Fatal        MXORD is not positive.
C        15   Fatal        From DDRIV1: N is greater than 200.
C        16   Fatal        From DDRIV1: The WORK array is too small.
C
C    B. The first three elements of WORK and the first five elements of
C       IWORK will contain the following statistical data:
C         AVGH     The average step size used.
C         HUSED    The step size last used (successfully).
C         AVGORD   The average order used.
C         IMXERR   The index of the element of the solution vector that
C                  contributed most to the last error test.
C         NQUSED   The order last used (successfully).
C         NSTEP    The number of steps taken.
C         NFE      The number of evaluations of the right hand side.
C         NJE      The number of evaluations of the Jacobian matrix.
C
C  IV.  REMARKS  .......................................................
C
C    A. Other routines used:
C         DDNTP, DDZRO, DDSTP, DDNTL, DDPST, DDCOR, DDCST,
C         DDPSC, and DDSCL;
C         DGEFA, DGESL, DGBFA, DGBSL, and DNRM2 (from LINPACK)
C         D1MACH (from the Bell Laboratories Machine Constants Package)
C         XERROR (from the SLATEC Common Math Library)
C       The last seven routines above, not having been written by the
C       present authors, are not explicitly part of this package.
C
C    B. On any return from DDRIV3 all information necessary to continue
C       the calculation is contained in the call sequence parameters,
C       including the work arrays.  Thus it is possible to suspend one
C       problem, integrate another, and then return to the first.
C
C    C. There are user-written routines which are only required by
C       DDRIV3 when certain parameters are set.  Thus a message warning
C       of unsatisfied externals may be issued during the load or link
C       phase.  This message should never refer to F.  This message can
C       be ignored if: it refers to JACOBN and MITER is not 1 or 4, or
C       it refers to FA and IMPL is 0 or MITER is 3, or it refers to
C       USERS and MITER is not 3, or it refers to G and NROOT is 0.
C
C    D. If this package is to be used in an overlay situation, the user
C       must declare in the primary overlay the variables in the call
C       sequence to DDRIV3.
C
C    E. Changing parameters during an integration.
C       The value of NROOT, EPS, EWT, IERROR, MINT, MITER, or HMAX may
C       be altered by the user between calls to DDRIV3.  For example, if
C       too much accuracy has been requested (the program returns with
C       NSTATE = 4 and an increased value of EPS) the user may wish to
C       increase EPS further.  In general, prudence is necessary when
C       making changes in parameters since such changes are not
C       implemented until the next integration step, which is not
C       necessarily the next call to DDRIV3.  This can happen if the
C       program has already integrated to a point which is beyond the
C       new point TOUT.
C
C    F. As the price for complete control of matrix algebra, the DDRIV3
C       USERS option puts all responsibility for Jacobian matrix
C       evaluation on the user.  It is often useful to approximate
C       numerically all or part of the Jacobian matrix.  However this
C       must be done carefully.  The FORTRAN sequence below illustrates
C       the method we recommend.  It can be inserted directly into
C       subroutine USERS to approximate Jacobian elements in rows I1
C       to I2 and columns J1 to J2.
C              DOUBLE PRECISION DFDY(N,N), EPSJ, H, R, D1MACH,
C             8     SAVE1(N), SAVE2(N), T, UROUND, Y(N), YJ, YWT(N)
C              UROUND = D1MACH(4)
C              EPSJ = UROUND**(1.D0/3.D0)
C              DO 30 J = J1,J2
C                R = EPSJ*MAX(ABS(YWT(J)), ABS(Y(J)))
C                IF (R .EQ. 0.D0) R = EPSJ
C                YJ = Y(J)
C                Y(J) = Y(J) + R
C                CALL F (N, T, Y, SAVE1)
C                Y(J) = YJ
C                DO 20 I = I1,I2
C         20       DFDY(I,J) = (SAVE1(I) - SAVE2(I))/R
C         30     CONTINUE
C       Many problems give rise to structured sparse Jacobians, e.g.,
C       block banded.  It is possible to approximate them with fewer
C       function evaluations than the above procedure uses; see Curtis,
C       Powell and Reid, J. Inst. Maths Applics, (1974), Vol. 13,
C       pp. 117-119.
C
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  DDSTP,SDNTP,SDZRO,DGEFA,DGESL,DGBFA,DGBSL,DNRM2,
C                    D1MACH,XERROR
C***END PROLOGUE  DDRIV3
      EXTERNAL F, JACOBN, FA, G
      DOUBLE PRECISION AE, BIG, EPS, EWT(*), G, GLAST, H, HMAX, HSIGN,
     8     NROUND, RE, D1MACH, SIZE, DNRM2, SUM, T, TLAST, TOUT, TROOT,
     8     UROUND, WORK(*), Y(*)
      INTEGER IWORK(*)
      LOGICAL CONVRG
      CHARACTER MSG*205
      PARAMETER(NROUND = 20.D0)
      PARAMETER(IAVGH = 1, IHUSED = 2, IAVGRD = 3,
     8          IEL = 4, IH = 160, IHMAX = 161, IHOLD = 162,
     8          IHSIGN = 163, IRC = 164, IRMAX = 165, IT = 166,
     8          ITOUT = 167, ITQ = 168, ITREND = 204, IYH = 205,
     8          INDMXR = 1, INQUSD = 2, INSTEP = 3, INFE = 4, INJE = 5,
     8          INROOT = 6, ICNVRG = 7, IJROOT = 8, IJTASK = 9,
     8          IMNTLD = 10, IMTRLD = 11, INQ = 12, INRTLD = 13,
     8          INDTRT = 14, INWAIT = 15, IMNT = 16, IMTRSV = 17,
     8          IMTR = 18, IMXRDS = 19, IMXORD = 20)
      PARAMETER(INDPRT = 21, INDPVT = 22)
C***FIRST EXECUTABLE STATEMENT  DDRIV3
      UROUND = D1MACH (4)
      IF (NROOT .NE. 0) THEN
        AE = D1MACH(1)
        RE = UROUND
      END IF
      IF (EPS .LT. 0.D0) THEN
        WRITE(MSG, '(''DDRIV36FE Illegal input.  EPS,'', D16.8,
     8  '', is negative.'')') EPS
        CALL XERROR(MSG, 60, 6, 2)
        RETURN
      END IF
      IF (N .LE. 0) THEN
        WRITE(MSG, '(''DDRIV37FE Illegal input.  Number of equations,'',
     8  I8, '', is not positive.'')') N
        CALL XERROR(MSG, 72, 7, 2)
        RETURN
      END IF
      IF (MXORD .LE. 0) THEN
        WRITE(MSG, '(''DDRIV314FE Illegal input.  Maximum order,'', I8,
     8  '', is not positive.'')') MXORD
        CALL XERROR(MSG, 67, 14, 2)
        RETURN
      END IF
      IF ((MINT .LT. 1 .OR. MINT .GT. 3) .OR. (MINT .EQ. 3 .AND.
     8  (MITER .EQ. 0 .OR. MITER .EQ. 3 .OR. IMPL .NE. 0))
     8  .OR. (MITER .LT. 0 .OR. MITER .GT. 5) .OR.
     8  (IMPL .NE. 0 .AND. IMPL .NE. 1 .AND. IMPL .NE. 2) .OR.
     8  ((IMPL .EQ. 1 .OR. IMPL .EQ. 2) .AND. MITER .EQ. 0) .OR.
     8  (IMPL .EQ. 2 .AND. MINT .EQ. 1)) THEN
        WRITE(MSG, '(''DDRIV39FE Illegal input.  Improper value for '',
     8  ''MINT, MITER and/or IMPL.'')')
        CALL XERROR(MSG, 69, 9, 2)
        RETURN
      END IF
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        LIWCHK = INDPVT - 1
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2 .OR. MITER .EQ. 4 .OR.
     8  MITER .EQ. 5) THEN
        LIWCHK = INDPVT + N - 1
      END IF
      IF (LENIW .LT. LIWCHK) THEN
        WRITE(MSG, '(''DDRIV310FE Illegal input.  Insufficient '',
     8  ''storage allocated for the IWORK array.  Based on the '')')
        WRITE(MSG(94:), '(''value of the input parameters involved, '',
     8  ''the required storage is'', I8)') LIWCHK
        CALL XERROR(MSG, 164, 10, 2)
        RETURN
      END IF
C                                                Allocate the WORK array
C                                         IYH is the index of YH in WORK
      IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN
        MAXORD = MIN(MXORD, 12)
      ELSE IF (MINT .EQ. 2) THEN
        MAXORD = MIN(MXORD, 5)
      END IF
      IDFDY = IYH + (MAXORD + 1)*N
C                                             IDFDY is the index of DFDY
C
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3)  THEN
        IYWT = IDFDY
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2)  THEN
        IYWT = IDFDY + N*N
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5)  THEN
        IYWT = IDFDY + (2*ML + MU + 1)*N
      END IF
C                                               IYWT is the index of YWT
      ISAVE1 = IYWT + N
C                                           ISAVE1 is the index of SAVE1
      ISAVE2 = ISAVE1 + N
C                                           ISAVE2 is the index of SAVE2
      IGNOW = ISAVE2 + N
C                                             IGNOW is the index of GNOW
      ITROOT = IGNOW + NROOT
C                                           ITROOT is the index of TROOT
      IA = ITROOT + NROOT
C                                                   IA is the index of A
      IF (IMPL .EQ. 0 .OR. MITER .EQ. 3) THEN
        LENCHK = IA - 1
      ELSE IF (IMPL .EQ. 1 .AND. (MITER .EQ. 1 .OR. MITER .EQ. 2)) THEN
        LENCHK = IA - 1 + N*N
      ELSE IF (IMPL .EQ. 1 .AND. (MITER .EQ. 4 .OR. MITER .EQ. 5)) THEN
        LENCHK = IA - 1 + (2*ML + MU + 1)*N
      ELSE IF (IMPL .EQ. 2 .AND. MITER .NE. 3) THEN
        LENCHK = IA - 1 + N
      END IF
      IF (LENW .LT. LENCHK) THEN
        WRITE(MSG, '(''DDRIV38FE Illegal input.  Insufficient '',
     8  ''storage allocated for the WORK array.  Based on the '')')
        WRITE(MSG(92:), '(''value of the input parameters involved, '',
     8  ''the required storage is'', I8)') LENCHK
        CALL XERROR(MSG, 162, 8, 2)
        RETURN
      END IF
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        MATDIM = 1
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        MATDIM = N
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        MATDIM = 2*ML + MU + 1
      END IF
      IF (IMPL .EQ. 0 .OR. IMPL .EQ. 1) THEN
        NDECOM = N
      ELSE IF (IMPL .EQ. 2) THEN
        NDECOM = NDE
      END IF
      IF (NSTATE .EQ. 1) THEN
C                                                  Initialize parameters
        IF (T .EQ. TOUT) RETURN
        IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN
          IWORK(IMXORD) = MIN(MXORD, 12)
        ELSE IF (MINT .EQ. 2) THEN
          IWORK(IMXORD) = MIN(MXORD, 5)
        END IF
        IWORK(IMXRDS) = MXORD
        IF (MINT .EQ. 1 .OR. MINT .EQ. 2) THEN
          IWORK(IMNT) = MINT
          IWORK(IMTR) = MITER
          IWORK(IMNTLD) = MINT
          IWORK(IMTRLD) = MITER
        ELSE IF (MINT .EQ. 3) THEN
          IWORK(IMNT) = 1
          IWORK(IMTR) = 0
          IWORK(IMNTLD) = IWORK(IMNT)
          IWORK(IMTRLD) = IWORK(IMTR)
          IWORK(IMTRSV) = MITER
        END IF
        WORK(IHMAX) = HMAX
        H = (TOUT - T)*(1.D0 - 4.D0*UROUND)
        H = SIGN(MIN(ABS(H), HMAX), H)
        WORK(IH) = H
        HSIGN = SIGN(1.D0, H)
        WORK(IHSIGN) = HSIGN
        IWORK(IJTASK) = 0
        WORK(IAVGH) = 0.D0
        WORK(IAVGRD) = 0.D0
        IWORK(INQUSD) = 0
        IWORK(INSTEP) = 0
        IWORK(INFE) = 0
        IWORK(INJE) = 0
        WORK(IT) = T
        IWORK(ICNVRG) = 0
        IWORK(INDPRT) = 0
C                                                 Set initial conditions
        DO 30 I = 1,N
          JYH = I + IYH - 1
 30       WORK(JYH) = Y(I)
        GO TO 180
      END IF
C                                             On a continuation, check
C                                             that output points have
C                                             been or will be overtaken.
      IF (IWORK(ICNVRG) .EQ. 1) THEN
        CONVRG = .TRUE.
      ELSE
        CONVRG = .FALSE.
      END IF
      T = WORK(IT)
      H = WORK(IH)
      HSIGN = WORK(IHSIGN)
      IF (IWORK(IJTASK) .EQ. 0) GO TO 180
C
C                                   IWORK(IJROOT) flags unreported
C                                   roots, and is set to the value of
C                                   NTASK when a root was last selected.
C                                   It is set to zero when all roots
C                                   have been reported.  IWORK(INROOT)
C                                   contains the index and WORK(ITOUT)
C                                   contains the value of the root last
C                                   selected to be reported.
C                                   IWORK(INRTLD) contains the value of
C                                   NROOT and IWORK(INDTRT) contains
C                                   the value of ITROOT when the array
C                                   of roots was last calculated.
      IF(NROOT .NE. 0) THEN
        JROOT = IWORK(IJROOT)
        IF (JROOT .GT. 0) THEN
C                                      TOUT has just been reported.
C                                      If TROOT .LE. TOUT, report TROOT.
          IF (NSTATE .NE. 5) THEN
            IF (TOUT*HSIGN .GE. WORK(ITOUT)*HSIGN) THEN
              TROOT = WORK(ITOUT)
              CALL DDNTP(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
              T = TROOT
              NSTATE = 5
              GO TO 580
            END IF
C                                         A root has just been reported.
C                                         Select the next root.
          ELSE
            TROOT = T
            IROOT = 0
            DO 50 I = 1,IWORK(INRTLD)
              JTROOT = IWORK(INDTRT) + I - 1
              IF (WORK(JTROOT)*HSIGN .LE. TROOT*HSIGN) THEN
C
C                                              Check for multiple roots.
C
                IF (WORK(JTROOT) .EQ. WORK(ITOUT) .AND.
     8          I .GT. IWORK(INROOT)) THEN
                  IROOT = I
                  TROOT = WORK(JTROOT)
                  GO TO 60
                END IF
                IF (WORK(JTROOT)*HSIGN .GT. WORK(ITOUT)*HSIGN) THEN
                  IROOT = I
                  TROOT = WORK(JTROOT)
                END IF
              END IF
 50           CONTINUE
 60         IWORK(INROOT) = IROOT
            WORK(ITOUT) = TROOT
            IWORK(IJROOT) = NTASK
            IF (NTASK .EQ. 1) THEN
              IF (IROOT .EQ. 0) THEN
                IWORK(IJROOT) = 0
              ELSE
                IF (TOUT*HSIGN .GE. TROOT*HSIGN) THEN
                  CALL DDNTP(H, 0, N, IWORK(INQ), T, TROOT,WORK(IYH),Y)
                  NSTATE = 5
                  T = TROOT
                  GO TO 580
                END IF
              END IF
            ELSE IF (NTASK .EQ. 2 .OR. NTASK .EQ. 3) THEN
C
C                                     If there are no more roots, or the
C                                     user has altered TOUT to be less
C                                     than a root, set IJROOT to zero.
C
              IF (IROOT .EQ. 0 .OR. (TOUT*HSIGN .LT. TROOT*HSIGN)) THEN
                IWORK(IJROOT) = 0
              ELSE
                CALL DDNTP(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH), Y)
                NSTATE = 5
                T = TROOT
                GO TO 580
              END IF
            END IF
          END IF
        END IF
      END IF
C
      IF (NTASK .EQ. 1) THEN
        NSTATE = 2
        IF (T*HSIGN .GE. TOUT*HSIGN) THEN
          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
      ELSE IF (NTASK .EQ. 2) THEN
C                                                      Check if TOUT has
C                                                      been reset .LT. T
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
          WRITE(MSG, '(''DDRIV32WRN With NTASK='', I1, '' on input, '',
     8    ''T,'', D16.8, '', was beyond TOUT,'', D16.8, ''.  Solution'',
     8    '' obtained by interpolation.'')') NTASK, T, TOUT
          CALL XERROR(MSG, 124, 2, 0)
          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          NSTATE = 2
          GO TO 580
        END IF
C                                   Determine if TOUT has been overtaken
C
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
          NSTATE = 2
          GO TO 560
        END IF
C                                             If there are no more roots
C                                             to report, report T.
        IF (NSTATE .EQ. 5) THEN
          NSTATE = 2
          GO TO 560
        END IF
        NSTATE = 2
C                                                       See if TOUT will
C                                                       be overtaken.
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
          H = TOUT - T
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
          WORK(IH) = H
          IF (H .EQ. 0.D0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        NSTATE = 2
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
          WRITE(MSG, '(''DDRIV32WRN With NTASK='', I1, '' on input, '',
     8    ''T,'', D16.8, '', was beyond TOUT,'', D16.8, ''.  Solution'',
     8    '' obtained by interpolation.'')') NTASK, T, TOUT
          CALL XERROR(MSG, 124, 2, 0)
          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
          GO TO 560
        END IF
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
          H = TOUT - T
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
          WORK(IH) = H
          IF (H .EQ. 0.D0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      END IF
C                         Implement changes in MINT, MITER, and/or HMAX.
C
      IF ((MINT .NE. IWORK(IMNTLD) .OR. MITER .NE. IWORK(IMTRLD)) .AND.
     8  MINT .NE. 3 .AND. IWORK(IMNTLD) .NE. 3) IWORK(IJTASK) = -1
      IF (HMAX .NE. WORK(IHMAX)) THEN
        H = SIGN(MIN(ABS(H), HMAX), H)
        IF (H .NE. WORK(IH)) THEN
          IWORK(IJTASK) = -1
          WORK(IH) = H
        END IF
        WORK(IHMAX) = HMAX
      END IF
C
 180  NSTEPL = IWORK(INSTEP)
      DO 190 I = 1,N
        JYH = IYH + I - 1
 190    Y(I) = WORK(JYH)
      IF (NROOT .NE. 0) THEN
        DO 200 I = 1,NROOT
          JGNOW = IGNOW + I - 1
 200      WORK(JGNOW) = G (N, T, Y, I)
      END IF
      IF (IERROR .EQ. 1) THEN
        DO 230 I = 1,N
          JYWT = I + IYWT - 1
 230      WORK(JYWT) = 1.D0
        GO TO 410
      ELSE IF (IERROR .EQ. 5) THEN
        DO 250 I = 1,N
          JYWT = I + IYWT - 1
 250      WORK(JYWT) = EWT(I)
        GO TO 410
      END IF
C                                       Reset YWT array.  Looping point.
 260  IF (IERROR .EQ. 2) THEN
        DO 280 I = 1,N
          IF (Y(I) .EQ. 0.D0) GO TO 290
          JYWT = I + IYWT - 1
 280      WORK(JYWT) = ABS(Y(I))
        GO TO 410
 290    IF (IWORK(IJTASK) .EQ. 0) THEN
          CALL F (N, T, Y, WORK(ISAVE2))
          IWORK(INFE) = IWORK(INFE) + 1
          IF (MITER .EQ. 3 .AND. IMPL .NE. 0) THEN
            IFLAG = 0
            CALL USERS(Y,WORK(IYH),WORK(IYWT),WORK(ISAVE1),WORK(ISAVE2),
     8                 T, H, WORK(IEL), IMPL, N, NDECOM, IFLAG)
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
              CALL DGEFA (WORK(IA), MATDIM, N, IWORK(INDPVT), INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGESL(WORK(IA),MATDIM,N,IWORK(INDPVT),WORK(ISAVE2),0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              JAML = IA + ML
              CALL FA (N, T, Y, WORK(JAML), MATDIM, ML, MU, NDECOM)
              CALL DGBFA (WORK(IA),MATDIM,N,ML,MU,IWORK(INDPVT),INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGBSL (WORK(IA), MATDIM, N, ML, MU, IWORK(INDPVT),
     8                    WORK(ISAVE2), 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (N, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
            DO 340 I = 1,NDECOM
              JA = I + IA - 1
              JSAVE2 = I + ISAVE2 - 1
              IF(WORK(JA) .EQ. 0.D0) GO TO 690
 340          WORK(JSAVE2) = WORK(JSAVE2)/WORK(JA)
          END IF
        END IF
        DO 360 J = I,N
          JYWT = J + IYWT - 1
          IF (Y(J) .NE. 0.D0) THEN
            WORK(JYWT) = ABS(Y(J))
          ELSE
            IF (IWORK(IJTASK) .EQ. 0) THEN
              JSAVE2 = J + ISAVE2 - 1
              WORK(JYWT) = ABS(H*WORK(JSAVE2))
            ELSE
              JHYP = J + IYH + N - 1
              WORK(JYWT) = ABS(WORK(JHYP))
            END IF
          END IF
          IF (WORK(JYWT) .EQ. 0.D0) WORK(JYWT) = UROUND
 360      CONTINUE
      ELSE IF (IERROR .EQ. 3) THEN
        DO 380 I = 1,N
          JYWT = I + IYWT - 1
 380      WORK(JYWT) = MAX(EWT(1), ABS(Y(I)))
      ELSE IF (IERROR .EQ. 4) THEN
        DO 400 I = 1,N
          JYWT = I + IYWT - 1
 400      WORK(JYWT) = MAX(EWT(I), ABS(Y(I)))
      END IF
C
 410  DO 420 I = 1,N
        JYWT = I + IYWT - 1
        JSAVE2 = I + ISAVE2 - 1
 420    WORK(JSAVE2) = Y(I)/WORK(JYWT)
      SUM = DNRM2(N, WORK(ISAVE2), 1)/SQRT(DBLE(N))
      IF (EPS .LT. SUM*UROUND) THEN
        EPS = SUM*UROUND*(1.D0 + 10.D0*UROUND)
        WRITE(MSG, '(''DDRIV34REC At T,'', D16.8, '', the requested '',
     8  ''accuracy, EPS, was not obtainable with the machine '',
     8  ''precision.  EPS has been increased to'')') T
        WRITE(MSG(137:), '(D16.8)') EPS
        CALL XERROR(MSG, 152, 4, 1)
        NSTATE = 4
        GO TO 560
      END IF
      IF (ABS(H) .GE. UROUND*ABS(T)) THEN
        IWORK(INDPRT) = 0
      ELSE IF (IWORK(INDPRT) .EQ. 0) THEN
        WRITE(MSG, '(''DDRIV35WRN At T,'', D16.8, '', the step size,'',
     8  D16.8, '', is smaller than the roundoff level of T.  '')') T, H
        WRITE(MSG(109:), '(''This may occur if there is an abrupt '',
     8  ''change in the right hand side of the differential '',
     8  ''equations.'')')
        CALL XERROR(MSG, 205, 5, 0)
        IWORK(INDPRT) = 1
      END IF
      IF (NTASK.NE.2) THEN
        IF ((IWORK(INSTEP)-NSTEPL) .GT. MXSTEP) THEN
          WRITE(MSG, '(''DDRIV33WRN At T,'', D16.8, '', '', I8,
     8    '' steps have been taken without reaching TOUT,'', D16.8)')
     8    T, MXSTEP, TOUT
          CALL XERROR(MSG, 103, 3, 0)
          NSTATE = 3
          GO TO 560
        END IF
      END IF
C
C     CALL DDSTP (EPS, F, FA, HMAX, IMPL, JACOBN, MATDIM, MAXORD,
C    8            MINT, MITER, ML, MU, N, NDE, YWT, UROUND,
C    8            AVGH, AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD,
C    8            NFE, NJE, NQUSED, NSTEP, T, Y, YH,  A, CONVRG,
C    8            DFDY, EL, HOLD, IPVT, JSTATE, NQ, NWAIT, RC,
C    8            RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG, MTRSV, MXRDSV)
C
      CALL DDSTP (EPS, F, FA, WORK(IHMAX), IMPL, JACOBN, MATDIM,
     8            IWORK(IMXORD), IWORK(IMNT), IWORK(IMTR), ML, MU, N,
     8            NDECOM, WORK(IYWT), UROUND, WORK(IAVGH), WORK(IAVGRD),
     8            WORK(IH), WORK(IHUSED), IWORK(IJTASK), IWORK(IMNTLD),
     8            IWORK(IMTRLD), IWORK(INFE), IWORK(INJE),
     8            IWORK(INQUSD), IWORK(INSTEP), WORK(IT), Y, WORK(IYH),
     8            WORK(IA), CONVRG, WORK(IDFDY), WORK(IEL), WORK(IHOLD),
     8            IWORK(INDPVT), JSTATE, IWORK(INQ), IWORK(INWAIT),
     8            WORK(IRC), WORK(IRMAX), WORK(ISAVE1), WORK(ISAVE2),
     8            WORK(ITQ), WORK(ITREND), MINT, IWORK(IMTRSV),
     8            IWORK(IMXRDS))
      T = WORK(IT)
      H = WORK(IH)
      GO TO (470, 670, 680, 690), JSTATE
 470  IWORK(IJTASK) = 1
C                                 Determine if a root has been overtaken
      IF (ABS(H) .GE. UROUND*ABS(T) .AND. NROOT .NE. 0) THEN
        IROOT = 0
        DO 500 I = 1,NROOT
          JTROOT = ITROOT + I - 1
          JGNOW = IGNOW + I - 1
          GLAST = WORK(JGNOW)
          WORK(JGNOW) = G (N, T, Y, I)
          IF (GLAST*WORK(JGNOW) .GT. 0.D0) THEN
            WORK(JTROOT) = T + H
          ELSE
            IF (WORK(JGNOW) .EQ. 0.D0) THEN
              WORK(JTROOT) = T
              IROOT = I
            ELSE
              IF (GLAST .EQ. 0.D0) THEN
                WORK(JTROOT) = T + H
              ELSE
                TLAST = T - WORK(IHUSED)
                IROOT = I
                TROOT = T
                CALL DDZRO (AE, G, H, N, IWORK(INQ), IROOT, RE, T,
     8                      WORK(IYH), UROUND,  TROOT, TLAST,
     8                      WORK(JGNOW), GLAST,  WORK(ISAVE1))
                WORK(JTROOT) = TROOT
              END IF
            END IF
          END IF
 500      CONTINUE
        IF (IROOT .EQ. 0) THEN
          IWORK(IJROOT) = 0
C                                                  Select the first root
        ELSE
          IWORK(IJROOT) = NTASK
          IWORK(INRTLD) = NROOT
          IWORK(INDTRT) = ITROOT
          TROOT = T + H
          DO 510 I = 1,NROOT
            JTROOT = ITROOT + I - 1
            IF (WORK(JTROOT)*HSIGN .LT. TROOT*HSIGN) THEN
              TROOT = WORK(JTROOT)
              IROOT = I
            END IF
 510        CONTINUE
          IWORK(INROOT) = IROOT
          WORK(ITOUT) = TROOT
          IF (TROOT*HSIGN .LE. TOUT*HSIGN) THEN
            CALL DDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
            NSTATE = 5
            T = TROOT
            GO TO 580
          END IF
        END IF
      END IF
C                               Test for NTASK condition to be satisfied
      NSTATE = 2
      IF (NTASK .EQ. 1) THEN
        IF (T*HSIGN .LT. TOUT*HSIGN) GO TO 260
        CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
        T = TOUT
        GO TO 580
C                               TOUT is assumed to have been attained
C                               exactly if T is within twenty roundoff
C                               units of TOUT, relative to max(TOUT, T).
      ELSE IF (NTASK .EQ. 2) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
            WORK(IH) = H
            IF (H .EQ. 0.D0) GO TO 670
            IWORK(IJTASK) = -1
          END IF
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
            WORK(IH) = H
            IF (H .EQ. 0.D0) GO TO 670
            IWORK(IJTASK) = -1
          END IF
          GO TO 260
        END IF
      END IF
C                                      All returns are made through this
C                                      section.  IMXERR is determined.
 560  DO 570 I = 1,N
        JYH = I + IYH - 1
 570    Y(I) = WORK(JYH)
 580  IF (CONVRG) THEN
        IWORK(ICNVRG) = 1
      ELSE
        IWORK(ICNVRG) = 0
      END IF
      IF (IWORK(IJTASK) .EQ. 0) RETURN
      BIG = 0.D0
      IMXERR = 1
      IWORK(INDMXR) = IMXERR
      DO  590 I = 1,N
C                                            SIZE = ABS(ERROR(I)/YWT(I))
        JYWT = I + IYWT - 1
        JERROR = I + ISAVE1 - 1
        SIZE = ABS(WORK(JERROR)/WORK(JYWT))
        IF (BIG .LT. SIZE) THEN
          BIG = SIZE
          IMXERR = I
          IWORK(INDMXR) = IMXERR
        END IF
 590    CONTINUE
      RETURN
C                                        Fatal errors are processed here
C
 670  WRITE(MSG, '(''DDRIV311FE At T,'', D16.8, '', the attempted '',
     8  ''step size has gone to zero.  Often this occurs if the '',
     8  ''problem setup is incorrect.'')') T
      CALL XERROR(MSG, 129, 11, 2)
      RETURN
C
 680  WRITE(MSG, '(''DDRIV312FE At T,'', D16.8, '', the step size has'',
     8  '' been reduced about 50 times without advancing the '')') T
      WRITE(MSG(103:), '(''solution.  Often this occurs if the '',
     8  ''problem setup is incorrect.'')')
      CALL XERROR(MSG, 165, 12, 2)
      RETURN
C
 690  WRITE(MSG, '(''DDRIV313FE At T,'', D16.8, '', while solving'',
     8  '' A*YDOT = F, A is singular.'')') T
      CALL XERROR(MSG, 74, 13, 2)
      RETURN
      END
