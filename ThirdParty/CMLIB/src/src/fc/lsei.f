      SUBROUTINE LSEI(W,MDW,ME,MA,MG,N,PRGOPT,X,RNORME,RNORML,MODE,WS,
     1   IP)
C***BEGIN PROLOGUE  LSEI
C***DATE WRITTEN   790701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C   000601  Modified AMAX1 to generic MAX (JEC)
C
C***CATEGORY NO.  K1A2A,D9
C***KEYWORDS  CONSTRAINED LEAST SQUARES,CURVE FITTING,DATA FITTING,
C             EQUALITY CONSTRAINTS,INEQUALITY CONSTRAINTS,
C             QUADRATIC PROGRAMMING
C***AUTHOR  HANSON, R. J., (SNLA)
C           HASKELL, K. H., (SNLA)
C***PURPOSE  Solve a linearly constrained least squares problem with
C            equality and inequality constraints, and optionally compute
C            a covariance matrix.
C***DESCRIPTION
C
C     DIMENSION W(MDW,N+1),PRGOPT(*),X(N),
C     WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2)
C     above, K=MAX(MA+MG,N).
C
C     Written by R. J. Hanson and K. H. Haskell.  For further math.
C     and algorithmic details see Sandia Laboratories Tech. Repts.
C     SAND-77-0552, (1978), and SAND-78-1290, (1979), or Math. Prog.,
C     Vol. 21 (1981) pp. 98-118 and ACM Trans. on Math. Software,
C     Sept. 1982.
C
C     Abstract
C
C     This subprogram solves a linearly constrained least squares
C     problem with both equality and inequality constraints, and, if the
C     user requests, obtains a covariance matrix of the solution
C     parameters.
C
C     Suppose there are given matrices E, A and G of respective
C     dimensions ME by N, MA by N and MG by N, and vectors F, B and H of
C     respective lengths ME, MA and MG.  This subroutine solves the
C     linearly constrained least squares problem
C
C                   EX = F, (E ME by N) (equations to be exactly
C                                       satisfied)
C                   AX = B, (A MA by N) (equations to be
C                                       approximately satisfied,
C                                       least squares sense)
C                   GX .GE. H,(G MG by N) (inequality constraints)
C
C     The inequalities GX .GE. H mean that every component of the
C     product GX must be .GE. the corresponding component of H.
C
C     In case the equality constraints cannot be satisfied, a
C     generalized inverse solution residual vector length is obtained
C     for F-EX.  This is the minimal length possible for F-EX.
C
C
C     Any values ME .GE. 0, MA .GE. 0, or MG .GE. 0 are permitted.  The
C     rank of the matrix E is estimated during the computation.  We call
C     this value KRANKE.  It is an output parameter in IP(1) defined
C     below.  Using a generalized inverse solution of EX=F, a reduced
C     least squares problem with inequality constraints is obtained.
C     The tolerances used in these tests for determining the rank
C     of E and the rank of the reduced least squares problem are
C     given in Sandia Tech. Rept. SAND-78-1290.  They can be
C     modified by the user if new values are provided in
C     the option list of the array PRGOPT(*).
C
C     The user must dimension all arrays appearing in the call list..
C     W(MDW,N+1),PRGOPT(*),X(N),WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2)
C     where K=MAX(MA+MG,N).  This allows for a solution of a range of
C     problems in the given working space.  The dimension of WS(*)
C     given is a necessary overestimate.  Once a particular problem
C     has been run, the output parameter IP(3) gives the actual
C     dimension required for that problem.
C
C     The parameters for LSEI( ) are
C
C     Input..
C
C     W(*,*),MDW,   The array W(*,*) is doubly subscripted with
C     ME,MA,MG,N    first dimensioning parameter equal to MDW.
C                   For this discussion let us call M = ME+MA+MG.  Then
C                   MDW must satisfy MDW .GE. M.  The condition
C                   MDW .LT. M is an error.
C
C                   The array W(*,*) contains the matrices and vectors
C
C                                  (E  F)
C                                  (A  B)
C                                  (G  H)
C
C                   in rows and columns 1,...,M and 1,...,N+1
C                   respectively.
C
C                   The integers ME, MA, and MG are the
C                   respective matrix row dimensions
C                   of E, A and G.  Each matrix has N columns.
C
C     PRGOPT(*)    This real-valued array is the option vector.
C                  If the user is satisfied with the nominal
C                  subprogram features set
C
C                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
C
C                  Otherwise PRGOPT(*) is a linked list consisting of
C                  groups of data of the following form
C
C                  LINK
C                  KEY
C                  DATA SET
C
C                  The parameters LINK and KEY are each one word.
C                  The DATA SET can be comprised of several words.
C                  The number of items depends on the value of KEY.
C                  The value of LINK points to the first
C                  entry of the next group of data within
C                  PRGOPT(*).  The exception is when there are
C                  no more options to change.  In that
C                  case, LINK=1 and the values KEY and DATA SET
C                  are not referenced.  The general layout of
C                  PRGOPT(*) is as follows.
C
C               ...PRGOPT(1)=link1 (link to first entry of next group)
C               .  PRGOPT(2)=key1 (key to the option change)
C               .  PRGOPT(3)=data value (data value for this change)
C               .       .
C               .       .
C               .       .
C               ...PRGOPT(LINK1)=link2 (link to the first entry of
C               .                       next group)
C               .  PRGOPT(LINK1+1)=key2 (key to the option change)
C               .  PRGOPT(LINK1+2)=data value
C               ...     .
C               .       .
C               .       .
C               ...PRGOPT(LINK)=1 (no more options to change)
C
C                  Values of LINK that are nonpositive are errors.
C                  A value of LINK .GT. NLINK=100000 is also an error.
C                  This helps prevent using invalid but positive
C                  values of LINK that will probably extend
C                  beyond the program limits of PRGOPT(*).
C                  Unrecognized values of KEY are ignored.  The
C                  order of the options is arbitrary and any number
C                  of options can be changed with the following
C                  restriction.  To prevent cycling in the
C                  processing of the option array, a count of the
C                  number of options changed is maintained.
C                  Whenever this count exceeds NOPT=1000, an error
C                  message is printed and the subprogram returns.
C
C                  Options..
C
C                  KEY=1
C                         Compute in W(*,*) the N by N
C                  covariance matrix of the solution variables
C                  as an output parameter.  Nominally the
C                  covariance matrix will not be computed.
C                  (This requires no user input.)
C                  The data set for this option is a single value.
C                  It must be nonzero when the covariance matrix
C                  is desired.  If it is zero, the covariance
C                  matrix is not computed.  When the covariance matrix
C                  is computed, the first dimensioning parameter
C                  of the array W(*,*) must satisfy MDW .GE. MAX0(M,N).
C
C                  KEY=10
C                         Suppress scaling of the inverse of the
C                  normal matrix by the scale factor RNORM**2/
C                  MAX(1, no. of degrees of freedom).  This option
C                  only applies when the option for computing the
C                  covariance matrix (KEY=1) is used.  With KEY=1 and
C                  KEY=10 used as options the unscaled inverse of the
C                  normal matrix is returned in W(*,*).
C                  The data set for this option is a single value.
C                  When it is nonzero no scaling is done.  When it is
C                  zero scaling is done.  The nominal case is to do
C                  scaling so if option (KEY=1) is used alone, the
C                  matrix will be scaled on output.
C
C                  KEY=2
C                         Scale the nonzero columns of the
C                         entire data matrix.
C                  (E)
C                  (A)
C                  (G)
C
C                  to have length one.  The data set for this
C                  option is a single value.  It must be
C                  nonzero if unit length column scaling
C                  is desired.
C
C                  KEY=3
C                         Scale columns of the entire data matrix
C                  (E)
C                  (A)
C                  (G)
C
C                  with a user-provided diagonal matrix.
C                  The data set for this option consists
C                  of the N diagonal scaling factors, one for
C                  each matrix column.
C
C                  KEY=4
C                         Change the rank determination tolerance for
C                  the equality constraint equations from
C                  the nominal value of SQRT(SRELPR).  This quantity can
C                  be no smaller than SRELPR, the arithmetic-
C                  storage precision.  The quantity SRELPR is the
C                  largest positive number such that T=1.+SRELPR
C                  satisfies T .EQ. 1.  The quantity used
C                  here is internally restricted to be at
C                  least SRELPR.  The data set for this option
C                  is the new tolerance.
C
C                  KEY=5
C                         Change the rank determination tolerance for
C                  the reduced least squares equations from
C                  the nominal value of SQRT(SRELPR).  This quantity can
C                  be no smaller than SRELPR, the arithmetic-
C                  storage precision.  The quantity used
C                  here is internally restricted to be at
C                  least SRELPR.  The data set for this option
C                  is the new tolerance.
C
C                  For example, suppose we want to change
C                  the tolerance for the reduced least squares
C                  problem, compute the covariance matrix of
C                  the solution parameters, and provide
C                  column scaling for the data matrix.  For
C                  these options the dimension of PRGOPT(*)
C                  must be at least N+9.  The Fortran statements
C                  defining these options would be as follows:
C
C                  PRGOPT(1)=4 (link to entry 4 in PRGOPT(*))
C                  PRGOPT(2)=1 (covariance matrix key)
C                  PRGOPT(3)=1 (covariance matrix wanted)
C
C                  PRGOPT(4)=7 (link to entry 7 in PRGOPT(*))
C                  PRGOPT(5)=5 (least squares equas.  tolerance key)
C                  PRGOPT(6)=... (new value of the tolerance)
C
C                  PRGOPT(7)=N+9 (link to entry N+9 in PRGOPT(*))
C                  PRGOPT(8)=3 (user-provided column scaling key)
C
C                  CALL SCOPY(N,D,1,PRGOPT(9),1)  (Copy the N
C                    scaling factors from the user array D(*)
C                    to PRGOPT(9)-PRGOPT(N+8))
C
C                  PRGOPT(N+9)=1 (no more options to change)
C
C                  The contents of PRGOPT(*) are not modified
C                  by the subprogram.
C                  The options for WNNLS( ) can also be included
C                  in this array.  The values of KEY recognized
C                  by WNNLS( ) are 6, 7 and 8.  Their functions
C                  are documented in the usage instructions for
C                  subroutine WNNLS( ).  Normally these options
C                  do not need to be modified when using LSEI( ).
C
C     IP(1),       The amounts of working storage actually
C     IP(2)        allocated for the working arrays WS(*) and
C                  IP(*), respectively.  These quantities are
C                  compared with the actual amounts of storage
C                  needed by LSEI( ).  Insufficient storage
C                  allocated for either WS(*) or IP(*) is an
C                  error.  This feature was included in LSEI( )
C                  because miscalculating the storage formulas
C                  for WS(*) and IP(*) might very well lead to
C                  subtle and hard-to-find execution errors.
C
C                  The length of WS(*) must be at least
C
C                  LW = 2*(ME+N)+K+(MG+2)*(N+7)
C
C                  where K = max(MA+MG,N)
C                  This test will not be made if IP(1).LE.0.
C
C                  The length of IP(*) must be at least
C
C                  LIP = MG+2*N+2
C                  This test will not be made if IP(2).LE.0.
C
C     Output..
C
C     X(*),RNORME,  The array X(*) contains the solution parameters
C     RNORML        if the integer output flag MODE = 0 or 1.
C                   The definition of MODE is given directly below.
C                   When MODE = 0 or 1, RNORME and RNORML
C                   respectively contain the residual vector
C                   Euclidean lengths of F - EX and B - AX.  When
C                   MODE=1 the equality constraint equations EX=F
C                   are contradictory, so RNORME .NE. 0.  The residual
C                   vector F-EX has minimal Euclidean length.  For
C                   MODE .GE. 2, none of these parameters are
C                   defined.
C
C     MODE          Integer flag that indicates the subprogram
C                   status after completion.  If MODE .GE. 2, no
C                   solution has been computed.
C
C                   MODE =
C
C                   0  Both equality and inequality constraints
C                      are compatible and have been satisfied.
C
C                   1  Equality constraints are contradictory.
C                      A generalized inverse solution of EX=F was used
C                      to minimize the residual vector length F-EX.
C                      In this sense, the solution is still meaningful.
C
C                   2  Inequality constraints are contradictory.
C
C                   3  Both equality and inequality constraints
C                      are contradictory.
C
C                   The following interpretation of
C                   MODE=1,2 or 3 must be made.  The
C                   sets consisting of all solutions
C                   of the equality constraints EX=F
C                   and all vectors satisfying GX .GE. H
C                   have no points in common.  (In
C                   particular this does not say that
C                   each individual set has no points
C                   at all, although this could be the
C                   case.)
C
C                   4  Usage error occurred.  The value
C                      of MDW is .LT. ME+MA+MG, MDW is
C                      .LT. N and a covariance matrix is
C                      requested, or the option vector
C                      PRGOPT(*) is not properly defined,
C                      or the lengths of the working arrays
C                      WS(*) and IP(*), when specified in
C                      IP(1) and IP(2) respectively, are not
C                      long enough.
C
C     W(*,*)        The array W(*,*) contains the N by N symmetric
C                   covariance matrix of the solution parameters,
C                   provided this was requested on input with
C                   the option vector PRGOPT(*) and the output
C                   flag is returned with MODE = 0 or 1.
C
C     IP(*)         The integer working array has three entries
C                   that provide rank and working array length
C                   information after completion.
C
C                      IP(1) = rank of equality constraint
C                              matrix.  Define this quantity
C                              as KRANKE.
C
C                      IP(2) = rank of reduced least squares
C                              problem.
C
C                      IP(3) = the amount of storage in the
C                              working array WS(*) that was
C                              actually used by the subprogram.
C                              The formula given above for the length
C                              of WS(*) is a necessary overestimate.
C                              If exactly the same problem matrices
C                              are used in subsequent executions,
C                              the declared dimension of WS(*) can
C                              be reduced to this output value.
C     User Designated
C     Working Arrays..
C
C     WS(*),IP(*)              These are respectively type real
C                              and type integer working arrays.
C                              Their required minimal lengths are
C                              given above.
C***REFERENCES  K.H. HASKELL AND R.J. HANSON, *AN ALGORITHM FOR
C                 LINEAR LEAST SQUARES PROBLEMS WITH EQUALITY AND
C                 NONNEGATIVITY CONSTRAINTS*, SAND77-0552, JUNE 1978.
C               K.H. HASKELL AND R.J. HANSON, *SELECTED ALGORITHMS FOR
C                 THE LINEARLY CONSTRAINED LEAST SQUARES PROBLEM--
C                 A USERS GUIDE*, SAND78-1290, AUGUST 1979.
C               K.H. HASKELL AND R.J. HANSON, *AN ALGORITHM FOR
C                 LINEAR LEAST SQUARES PROBLEMS WITH EQUALITY AND
C                 NONNEGATIVITY CONSTRAINTS*, MATH. PROG. 21 (1981),
C                 PP. 98-118.
C               R.J. HANSON AND K.H. HASKELL, *TWO ALGORITHMS FOR THE
C                 LINEARLY CONSTRAINED LEAST SQUARES PROBLEM*, ACM
C                 TRANS. ON MATH. SOFTWARE, SEPT. 1982.
C***ROUTINES CALLED  H12,LSI,SASUM,SAXPY,SCOPY,SDOT,SNRM2,SSCAL,SSWAP,
C                    XERROR,XERRWV
C***END PROLOGUE  LSEI
C
C     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
C     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
C     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
C     (START EDITING AT LINE WITH C++ IN COLS. 1-3.)
C     /REAL (12 BLANKS)/DOUBLE PRECISION/,/SASUM/DASUM/,/SDOT/DDOT/,
C     /SNRM2/DNRM2/,/ SQRT/ DSQRT/,/ ABS/ DABS/,
C     /SCOPY/DCOPY/,/SSCAL/DSCAL/,/SAXPY/DAXPY/,/SSWAP/DSWAP/,/E0/D0/,
C     /, DUMMY/,SNGL(DUMMY)/,/SRELPR/DRELPR/
C
C
C     SUBROUTINES CALLED
C
C     LSI           PART OF THIS PACKAGE.  SOLVES A
C                   CONSTRAINED LEAST SQUARES PROBLEM WITH
C                   INEQUALITY CONSTRAINTS.
C
C++
C     SDOT,SSCAL,   SUBROUTINES FROM THE BLAS PACKAGE.
C     SAXPY,SASUM,  SEE TRANS. MATH. SOFT., VOL. 5, NO. 3, P. 308.
C     SCOPY,SNRM2,
C     SSWAP
C
C     H12           SUBROUTINE TO CONSTRUCT AND APPLY A
C                   HOUSEHOLDER TRANSFORMATION.
C
C     XERROR        FROM SLATEC ERROR PROCESSING PACKAGE.
C                   THIS IS DOCUMENTED IN SANDIA TECH. REPT.,
C                   SAND78-1189.
C
C     SUBROUTINE LSEI(W,MDW,ME,MA,MG,N,PRGOPT,X,RNORME,RNORML,MODE,WS,
C    1 IP)
C
C     REVISED MAY, 1982.
C
      REAL             W(MDW,*), PRGOPT(*), X(*), WS(*), RNORME, RNORML
      INTEGER IP(3)
      REAL             DUMMY, ENORM, SRELPR, FNORM, GAM, HALF, ONE, RB,
     1RN, RNMAX, SIZE, SN, SNMAX, T, TAU, UJ, UP, VJ, XNORM, XNRME, ZERO
      REAL             SASUM, SDOT, SNRM2, SQRT, ABS
      LOGICAL COV
      DATA ZERO /0.E0/, SRELPR /0.E0/, HALF /0.5E0/, ONE /1.E0/
C
C     CHECK THAT ENOUGH STORAGE WAS ALLOCATED IN WS(*) AND IP(*).
C***FIRST EXECUTABLE STATEMENT  LSEI
      IF (.NOT.(IP(1).GT.0)) GO TO 20
      LCHK = 2*(ME+N) + MAX0(MA+MG,N) + (MG+2)*(N+7)
      IF (.NOT.(IP(1).LT.LCHK)) GO TO 10
      MODE = 4
      NERR = 2
      IOPT = 1
      CALL XERRWV( 'LSEI( ), INSUFFICIENT STORAGE ALLOCATED FOR WS(*), N
     1EED LW=I1 BELOW', 67, NERR, IOPT, 1, LCHK, 0, 0, DUMMY, DUMMY)
      RETURN
   10 CONTINUE
   20 IF (.NOT.(IP(2).GT.0)) GO TO 40
      LCHK = MG + 2*N + 2
      IF (.NOT.(IP(2).LT.LCHK)) GO TO 30
      MODE = 4
      NERR = 2
      IOPT = 1
      CALL XERRWV( 'LSEI( ), INSUFFICIENT STORAGE ALLOCATED FOR IP(*), N
     1EED LIP=I1 BELOW', 68, NERR, IOPT, 1, LCHK, 0, 0, DUMMY, DUMMY)
      RETURN
   30 CONTINUE
C
C     COMPUTE MACHINE PRECISION=SRELPR ONLY WHEN NECESSARY.
   40 IF (.NOT.(SRELPR.EQ.ZERO)) GO TO 70
c*** changed by RF Boisvert, 19-Feb-92  (fails on HP 9000 Series 300)
      srelpr = r1mach(4)
c      SRELPR = ONE
c   50 IF (ONE+SRELPR.EQ.ONE) GO TO 60
c      SRELPR = SRELPR*HALF
c      GO TO 50
c   60 SRELPR = SRELPR + SRELPR
C
C     COMPUTE NUMBER OF POSSIBLE RIGHT MULTIPLYING HOUSEHOLDER
C     TRANSFORMATIONS.
   70 M = ME + MA + MG
      MODE = 0
      IF (N.LE.0 .OR. M.LE.0) RETURN
      IF (.NOT.(MDW.LT.M)) GO TO 80
      NERR = 1
      IOPT = 1
      CALL XERROR( 'LSEI( ), MDW.LT.ME+MA+MG IS AN ERROR', 36, NERR,
     1 IOPT)
      MODE = 4
      RETURN
   80 NP1 = N + 1
      KRANKE = MIN0(ME,N)
      N1 = 2*KRANKE + 1
      N2 = N1 + N
C
C     PROCESS-OPTION-VECTOR
      ASSIGN 90 TO IGO990
      GO TO 480
   90 IF (.NOT.(COV .AND. MDW.LT.N)) GO TO 100
      NERR = 2
      IOPT = 1
      CALL XERROR( 'LSEI( ), MDW.LT.N, WHEN COV MATRIX NEEDED, IS AN ERR
     1OR', 54,    NERR, IOPT)
      MODE = 4
      RETURN
  100 L = KRANKE
C
C     COMPUTE NORM OF EQUALITY CONSTRAINT MATRIX AND RT SIDE.
      ENORM = ZERO
      DO 110 J=1,N
        ENORM = MAX(ENORM,SASUM(ME,W(1,J),1))
  110 CONTINUE
      FNORM = SASUM(ME,W(1,NP1),1)
      IF (.NOT.(L.GT.0)) GO TO 190
      SNMAX = ZERO
      RNMAX = ZERO
      DO 180 I=1,L
C
C     COMPUTE MAXIMUM RATIO OF VECTOR LENGTHS. PARTITION
C     IS AT COL. I.
        DO 150 K=I,ME
          SN = SDOT(N-I+1,W(K,I),MDW,W(K,I),MDW)
          RN = SDOT(I-1,W(K,1),MDW,W(K,1),MDW)
          IF (.NOT.(RN.EQ.ZERO .AND. SN.GT.SNMAX)) GO TO 120
          SNMAX = SN
          IMAX = K
          GO TO 140
  120     IF (.NOT.(K.EQ.I .OR. (SN*RNMAX.GT.RN*SNMAX))) GO TO 130
          SNMAX = SN
          RNMAX = RN
          IMAX = K
  130     CONTINUE
  140     CONTINUE
  150   CONTINUE
C
C     INTERCHANGE ROWS IF NECESSARY.
        IF (I.NE.IMAX) CALL SSWAP(NP1, W(I,1), MDW, W(IMAX,1), MDW)
        IF (.NOT.(SNMAX.GT.TAU**2*RNMAX)) GO TO 160
C
C     ELIMINATE ELEMS I+1,...,N IN ROW I.
        CALL H12(1, I, I+1, N, W(I,1), MDW, WS(I), W(I+1,1), MDW, 1,
     1   M-I)
        GO TO 170
  160   KRANKE = I - 1
        GO TO 200
  170   CONTINUE
  180 CONTINUE
  190 CONTINUE
  200 CONTINUE
C
C     SAVE DIAG. TERMS OF LOWER TRAP. MATRIX.
      CALL SCOPY(KRANKE, W, MDW+1, WS(KRANKE+1), 1)
C
C     USE HOUSEHOLDER TRANS FROM LEFT TO ACHIEVE KRANKE BY KRANKE UPPER
C     TRIANGULAR FORM.
      IF (.NOT.(KRANKE.GT.0 .AND. KRANKE.LT.ME)) GO TO 220
      DO 210 KK=1,KRANKE
        K = KRANKE + 1 - KK
C
C     APPLY TRANFORMATION TO MATRIX COLS. 1,...,K-1.
        CALL H12(1, K, KRANKE+1, ME, W(1,K), 1, UP, W, 1, MDW, K-1)
C
C     APPLY TO RT SIDE VECTOR.
        CALL H12(2, K, KRANKE+1, ME, W(1,K), 1, UP, W(1,NP1), 1, 1, 1)
  210 CONTINUE
  220 IF (.NOT.(KRANKE.GT.0)) GO TO 240
C
C     SOLVE FOR VARIABLES 1,...,KRANKE IN NEW COORDINATES.
      CALL SCOPY(KRANKE, W(1,NP1), 1, X, 1)
      DO 230 I=1,KRANKE
        X(I) = (X(I)-SDOT(I-1,W(I,1),MDW,X,1))/W(I,I)
  230 CONTINUE
C
C     COMPUTE RESIDUALS FOR REDUCED PROBLEM.
  240 MEP1 = ME + 1
      RNORML = ZERO
      IF (.NOT.(ME.LT.M)) GO TO 270
      DO 260 I=MEP1,M
        W(I,NP1) = W(I,NP1) - SDOT(KRANKE,W(I,1),MDW,X,1)
        SN = SDOT(KRANKE,W(I,1),MDW,W(I,1),MDW)
        RN = SDOT(N-KRANKE,W(I,KRANKE+1),MDW,W(I,KRANKE+1),MDW)
        IF (.NOT.(RN.LE.TAU**2*SN .AND. KRANKE.LT.N)) GO TO 250
        W(I,KRANKE+1) = ZERO
        CALL SCOPY(N-KRANKE, W(I,KRANKE+1), 0, W(I,KRANKE+1), MDW)
  250   CONTINUE
  260 CONTINUE
C
C     COMPUTE EQUAL. CONSTRAINT EQUAS. RESIDUAL LENGTH.
  270 RNORME = SNRM2(ME-KRANKE,W(KRANKE+1,NP1),1)
C
C     MOVE REDUCED PROBLEM DATA UPWARD IF KRANKE.LT.ME.
      IF (.NOT.(KRANKE.LT.ME)) GO TO 290
      DO 280 J=1,NP1
        CALL SCOPY(M-ME, W(ME+1,J), 1, W(KRANKE+1,J), 1)
  280 CONTINUE
C
C     COMPUTE SOLN OF REDUCED PROBLEM.
  290 CALL LSI(W(KRANKE+1,KRANKE+1), MDW, MA, MG, N-KRANKE, PRGOPT,
     1 X(KRANKE+1), RNORML, MODE, WS(N2), IP(2))
      IF (.NOT.(ME.GT.0)) GO TO 330
C
C     TEST FOR CONSISTENCY OF EQUALITY CONSTRAINTS.
      MDEQC = 0
      XNRME = SASUM(KRANKE,W(1,NP1),1)
      IF (RNORME.GT.TAU*(ENORM*XNRME+FNORM)) MDEQC = 1
      MODE = MODE + MDEQC
C
C     CHECK IF SOLN TO EQUAL. CONSTRAINTS SATISFIES INEQUAL.
C     CONSTRAINTS WHEN THERE ARE NO DEGREES OF FREEDOM LEFT.
      IF (.NOT.(KRANKE.EQ.N .AND. MG.GT.0)) GO TO 320
      XNORM = SASUM(N,X,1)
      MAPKE1 = MA + KRANKE + 1
      MEND = MA + KRANKE + MG
      DO 310 I=MAPKE1,MEND
        SIZE = SASUM(N,W(I,1),MDW)*XNORM + ABS(W(I,NP1))
        IF (.NOT.(W(I,NP1).GT.TAU*SIZE)) GO TO 300
        MODE = MODE + 2
        GO TO 450
  300   CONTINUE
  310 CONTINUE
  320 CONTINUE
  330 IF (.NOT.(KRANKE.GT.0)) GO TO 420
C
C     REPLACE DIAG. TERMS OF LOWER TRAP. MATRIX.
      CALL SCOPY(KRANKE, WS(KRANKE+1), 1, W, MDW+1)
C
C     REAPPLY TRANS TO PUT SOLN IN ORIGINAL COORDINATES.
      DO 340 II=1,KRANKE
        I = KRANKE + 1 - II
        CALL H12(2, I, I+1, N, W(I,1), MDW, WS(I), X, 1, 1, 1)
  340 CONTINUE
C
C     COMPUTE COV MATRIX OF EQUAL. CONSTRAINED PROBLEM.
      IF (.NOT.(COV)) GO TO 410
      DO 400 JJ=1,KRANKE
        J = KRANKE + 1 - JJ
        IF (.NOT.(J.LT.N)) GO TO 390
        RB = WS(J)*W(J,J)
        IF (RB.NE.ZERO) RB = ONE/RB
        JP1 = J + 1
        DO 350 I=JP1,N
          W(I,J) = SDOT(N-J,W(I,JP1),MDW,W(J,JP1),MDW)*RB
  350   CONTINUE
        GAM = SDOT(N-J,W(JP1,J),1,W(J,JP1),MDW)*RB
        GAM = HALF*GAM
        CALL SAXPY(N-J, GAM, W(J,JP1), MDW, W(JP1,J), 1)
        DO 370 I=JP1,N
          DO 360 K=I,N
            W(I,K) = W(I,K) + W(J,I)*W(K,J) + W(I,J)*W(J,K)
            W(K,I) = W(I,K)
  360     CONTINUE
  370   CONTINUE
        UJ = WS(J)
        VJ = GAM*UJ
        W(J,J) = UJ*VJ + UJ*VJ
        DO 380 I=JP1,N
          W(J,I) = UJ*W(I,J) + VJ*W(J,I)
  380   CONTINUE
        CALL SCOPY(N-J, W(J,JP1), MDW, W(JP1,J), 1)
  390   CONTINUE
  400 CONTINUE
  410 CONTINUE
C
C     APPLY THE SCALING TO THE COVARIANCE MATRIX.
  420 IF (.NOT.(COV)) GO TO 440
      DO 430 I=1,N
        L = N1 + I
        CALL SSCAL(N, WS(L-1), W(I,1), MDW)
        CALL SSCAL(N, WS(L-1), W(1,I), 1)
  430 CONTINUE
  440 CONTINUE
  450 CONTINUE
C
C     RESCALE SOLN. VECTOR.
      IF (.NOT.(MODE.LE.1)) GO TO 470
      DO 460 J=1,N
        L = N1 + J
        X(J) = X(J)*WS(L-1)
  460 CONTINUE
  470 IP(1) = KRANKE
      IP(3) = IP(3) + 2*KRANKE + N
      RETURN
  480 CONTINUE
C     TO PROCESS-OPTION-VECTOR
C
C     THE NOMINAL TOLERANCE USED IN THE CODE
C     FOR THE EQUALITY CONSTRAINT EQUATIONS.
      TAU = SQRT(SRELPR)
C
C     THE NOMINAL COLUMN SCALING USED IN THE CODE IS
C     THE IDENTITY SCALING.
      WS(N1) = ONE
      CALL SCOPY(N, WS(N1), 0, WS(N1), 1)
C
C     NO COVARIANCE MATRIX IS NOMINALLY COMPUTED.
      COV = .FALSE.
C
C     DEFINE BOUND FOR NUMBER OF OPTIONS TO CHANGE.
      NOPT = 1000
      NTIMES = 0
C
C     DEFINE BOUND FOR POSITIVE VALUES OF LINK.
      NLINK = 100000
      LAST = 1
      LINK = PRGOPT(1)
      IF (.NOT.(LINK.LE.0 .OR. LINK.GT.NLINK)) GO TO 490
      NERR = 3
      IOPT = 1
      CALL XERROR( 'LSEI( ) THE OPTION VECTOR IS UNDEFINED', 38, NERR,
     1 IOPT)
      MODE = 4
      RETURN
  490 IF (.NOT.(LINK.GT.1)) GO TO 540
      NTIMES = NTIMES + 1
      IF (.NOT.(NTIMES.GT.NOPT)) GO TO 500
      NERR = 3
      IOPT = 1
      CALL XERROR( 'LSEI( ). THE LINKS IN THE OPTION VECTOR ARE CYCLING.
     1', 52,      NERR, IOPT)
      MODE = 4
      RETURN
  500 KEY = PRGOPT(LAST+1)
      IF (KEY.EQ.1) COV = PRGOPT(LAST+2).NE.ZERO
      IF (.NOT.(KEY.EQ.2 .AND. PRGOPT(LAST+2).NE.ZERO)) GO TO 520
      DO 510 J=1,N
        T = SNRM2(M,W(1,J),1)
        IF (T.NE.ZERO) T = ONE/T
        L = N1 + J
        WS(L-1) = T
  510 CONTINUE
  520 IF (KEY.EQ.3) CALL SCOPY(N, PRGOPT(LAST+2), 1, WS(N1), 1)
      IF (KEY.EQ.4) TAU = MAX(SRELPR,PRGOPT(LAST+2))
      NEXT = PRGOPT(LINK)
      IF (.NOT.(NEXT.LE.0 .OR. NEXT.GT.NLINK)) GO TO 530
      NERR = 3
      IOPT = 1
      CALL XERROR( 'LSEI( ) THE OPTION VECTOR IS UNDEFINED', 38, NERR,
     1 IOPT)
      MODE = 4
      RETURN
  530 LAST = LINK
      LINK = NEXT
      GO TO 490
  540 DO 550 J=1,N
        L = N1 + J
        CALL SSCAL(M, WS(L-1), W(1,J), 1)
  550 CONTINUE
      GO TO 560
  560 GO TO IGO990, (90)
      END
