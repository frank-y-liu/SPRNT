      SUBROUTINE DDRIV1 (N,T,Y,TOUT,MSTATE,EPS,WORK,LENW)
C***BEGIN PROLOGUE  DDRIV1
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850924   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD,
C             DOUBLE PRECISION
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***PURPOSE  The function of DDRIV1 is to solve N (200 or fewer)
C            ordinary differential equations of the form
C            dY(I)/dT = F(Y(I),T), given the initial conditions
C            Y(I) = YI.  DDRIV1 uses double precision arithmetic.
C***DESCRIPTION
C
C  I.  CHOOSING THE CORRECT ROUTINE  ...................................
C
C     SDRIV
C     DDRIV
C     CDRIV
C           These are the generic names for three packages for solving
C           initial value problems for ordinary differential equations.
C           SDRIV uses single precision arithmetic.  DDRIV uses double
C           precision arithmetic.  CDRIV allows complex-valued
C           differential equations, integrated with respect to a single,
C           real, independent variable.
C
C    As an aid in selecting the proper program, the following is a
C    discussion of the important options or restrictions associated with
C    each program:
C
C      A. DDRIV1 should be tried first for those routine problems with
C         no more than 200 differential equations.  Internally this
C         routine has two important technical defaults:
C           1. Numerical approximation of the Jacobian matrix of the
C              right hand side is used.
C           2. The stiff solver option is used.
C         Most users of DDRIV1 should not have to concern themselves
C         with these details.
C
C      B. DDRIV2 should be considered for those problems for which
C         DDRIV1 is inadequate (SDRIV2 has no explicit restriction on
C         the number of differential equations.)  For example, DDRIV1
C         may have difficulty with problems having zero initial
C         conditions and zero derivatives.  In this case DDRIV2, with an
C         appropriate value of the parameter EWT, should perform more
C         efficiently.  DDRIV2 provides three important additional
C         options:
C           1. The nonstiff equation solver (as well as the stiff
C              solver) is available.
C           2. The root-finding option is available.
C           3. The program can dynamically select either the non-stiff
C              or the stiff methods.
C         Internally this routine also defaults to the numerical
C         approximation of the Jacobian matrix of the right hand side.
C
C      C. DDRIV3 is the most flexible, and hence the most complex, of
C         the programs.  Its important additional features include:
C           1. The ability to exploit band structure in the Jacobian
C              matrix.
C           2. The ability to solve some implicit differential
C              equations, i.e., those having the form:
C                   A(Y,T)*dY/dT = F(Y,T).
C           3. The option of integrating in the one step mode.
C           4. The option of allowing the user to provide a routine
C              which computes the analytic Jacobian matrix of the right
C              hand side.
C           5. The option of allowing the user to provide a routine
C              which does all the matrix algebra associated with
C              corrections to the solution components.
C
C  II.  ABSTRACT  ......................................................
C
C    The function of DDRIV1 is to solve N (200 or fewer) ordinary
C    differential equations of the form dY(I)/dT = F(Y(I),T), given the
C    initial conditions Y(I) = YI.  DDRIV1 is to be called once for each
C    output point.
C
C  III.  PARAMETERS  ...................................................
C
C       (REMEMBER--To run DDRIV1 correctly in double precision, ALL
C       non-integer arguments in the call sequence, including
C       arrays, MUST be declared double precision.)
C    The user should use parameter names in the call sequence of DDRIV1
C    for those quantities whose value may be altered by DDRIV1.  The
C    parameters in the call sequence are:
C
C    N      = (Input) The number of differential equations, N .LE. 200
C
C    T      = The independent variable.  On input for the first call, T
C             is the initial point.  On output, T is the point at which
C             the solution is given.
C
C    Y      = The vector of dependent variables.  Y is used as input on
C             the first call, to set the initial values.  On output, Y
C             is the computed solution vector.  This array Y is passed
C             in the call sequence of the user-provided routine F.
C
C    TOUT   = (Input) The point at which the solution is desired.
C
C    MSTATE = An integer describing the status of integration.  The user
C             must initialize MSTATE to +1 or -1.  If MSTATE is
C             positive, the routine will integrate past TOUT and
C             interpolate the solution.  This is the most efficient
C             mode.  If MSTATE is negative, the routine will adjust its
C             internal step to reach TOUT exactly (useful if a
C             singularity exists beyond TOUT.)  The meaning of the
C             magnitude of MSTATE:
C               1  (Input) Means the first call to the routine.  This
C                  value must be set by the user.  On all subsequent
C                  calls the value of MSTATE should be tested by the
C                  user.  Unless DDRIV1 is to be reinitialized, only the
C                  sign of MSTATE may be changed by the user.  (As a
C                  convenience to the user who may wish to put out the
C                  initial conditions, DDRIV1 can be called with
C                  MSTATE=+1(-1), and TOUT=T.  In this case the program
C                  will return with MSTATE unchanged, i.e.,
C                  MSTATE=+1(-1).)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  1000 steps without reaching TOUT.  The user can
C                  continue the integration by simply calling DDRIV1
C                  again.
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling DDRIV1
C                  again.
C
C    EPS    = On input, the requested relative accuracy in all solution
C             components.  On output, the adjusted relative accuracy if
C             the input value was too small.  The value of EPS should be
C             set as large as is reasonable, because the amount of work
C             done by DDRIV1 increases as EPS decreases.
C
C    WORK
C    LENW   = (Input)
C             WORK is an array of LENW double precision words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       DOUBLE PRECISION WORK(...)
C             The length of WORK should be at least N*N + 10*N + 225
C             and LENW should be set to the value used.  The contents of
C             WORK should not be disturbed between calls to DDRIV1.
C
C***LONG DESCRIPTION
C
C  IV.  USAGE  .........................................................
C
C                   PROGRAM SAMPLE
C                   DOUBLE PRECISION Y(...), WORK(...)
C                   OPEN(FILE='TAPE6', UNIT=6, STATUS='NEW')
C                   N = ...                       Number of equations
C                   T = ...                       Initial point
C                   DO 10 I = 1,N
C              10     Y(I) = ...                  Set initial conditions
C                   TOUT = T
C                   MSTATE = 1
C                   EPS = ...
C                   LENW = ...
C              20   CALL DDRIV1 (N, T, Y, TOUT, MSTATE, EPS, WORK, LENW)
C                   IF (MSTATE .GT. 2) STOP
C                   WRITE(6, 100) TOUT, (Y(I), I=1,N)
C                   TOUT = TOUT + 1.
C                   IF (TOUT .LE. 10.) GO TO 20
C              100  FORMAT(...)
C                   END (Sample)
C
C    The user must write a subroutine called F to evaluate the right
C    hand side of the differential equations.  It is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   DOUBLE PRECISION Y(*), YDOT(*)
C                     .
C                     .
C                   YDOT(1) = ...
C                     .
C                     .
C                   YDOT(N) = ...
C                   END (Sample)
C    This computes YDOT = F(Y,T), the right hand side of the
C    differential equations.  Here Y is a vector of length at least N.
C    The actual length of Y is determined by the user's declaration in
C    the program which calls DDRIV1.  Thus the dimensioning of Y in F,
C    while required by FORTRAN convention, does not actually allocate
C    any storage.  When this subroutine is called, the first N
C    components of Y are intermediate approximations to the solution
C    components.  The user should not alter these values.  Here YDOT is
C    a vector of length N.  The user should only compute YDOT(I) for I
C    from 1 to N.
C
C  V.  OTHER COMMUNICATION TO THE USER  ................................
C
C    The solver communicates to the user through the parameters above.
C    In addition it writes diagnostic messages through the standard
C    error handling program XERROR.  That program will terminate the
C    user's run if it detects a probable problem setup error, e.g.,
C    insufficient storage allocated by the user for the WORK array.  For
C    further information see section III-A of the writeup for DDRIV3.
C
C  VI.  REMARKS  .......................................................
C
C    A. There are user-written routines which are only required by
C       DDRIV2 or DDRIV3 when certain parameters are set.  Thus a
C       message warning of unsatisfied externals may be issued during
C       the load or link phase.  This message can be ignored unless it
C       refers to F.
C
C    For other information, see section IV of the writeup for DDRIV3.
C
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  DDRIV3,D1MACH,XERROR
C***END PROLOGUE  DDRIV1
      EXTERNAL F, JACOBN, FA, G
      DOUBLE PRECISION EPS, EWT, G, HMAX, D1MACH, T, TOUT,
     8     WORK(*), Y(*)
      PARAMETER(MXN = 200, IDLIW = 21, MXLIW = IDLIW + MXN)
      INTEGER IWORK(MXLIW)
      CHARACTER MSG*103
      PARAMETER(NROOT = 0, EWT = 1.D0, IERROR = 2, MINT = 2, MITER = 2,
     8          IMPL = 0, ML = 0, MU = 0, MXORD = 5, NDE = 0,
     8          MXSTEP = 1000)
C***FIRST EXECUTABLE STATEMENT  DDRIV1
      IF (N .GT. MXN) THEN
        WRITE(MSG, '(''DDRIV115FE Illegal input.  The number of '',
     8  ''equations,'', I8, '', is greater than the maximum allowed.'')
     8  ') N
        CALL XERROR(MSG, 97, 15, 2)
        RETURN
      END IF
      IF (MSTATE .GT. 0) THEN
        NSTATE = MSTATE
        NTASK = 1
      ELSE
        NSTATE = - MSTATE
        NTASK = 3
      END IF
      HMAX = SQRT(D1MACH(2))
      LENIW = N + IDLIW
      LENWCM = LENW - LENIW
      IF (LENWCM .LT. (N*N + 9*N + 204)) THEN
        LNWCHK = N*N + 9*N + 204 + LENIW
        WRITE(MSG, '(''DDRIV116FE Insufficient storage allocated for '',
     8  ''the work array.  The required storage is at least'', I8)')
     8  LNWCHK
        CALL XERROR(MSG, 103, 16, 2)
        RETURN
      END IF
      IF (NSTATE .NE. 1) THEN
        DO 20 I = 1,LENIW
          II = I + LENWCM
 20       IWORK(I) = INT(WORK(II))
      END IF
      CALL DDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS, EWT,
     8             IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,
     8             LENWCM, IWORK, LENIW, JACOBN, FA, NDE, MXSTEP, G)
      DO 40 I = 1,LENIW
        II = LENWCM + I
 40     WORK(II) = DBLE(IWORK(I))
      IF (MSTATE .GE. 0) THEN
        MSTATE = NSTATE
      ELSE
        MSTATE = - NSTATE
      END IF
      END
