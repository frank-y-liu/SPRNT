      SUBROUTINE DPOFS(A,LDA,N,V,ITASK,IND,WORK)
C***BEGIN PROLOGUE  DPOFS
C***DATE WRITTEN   800514   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY (YYMMDD)
C   000601  Changed references to IDINT to generic INT   (RFB)
C***CATEGORY NO.  D2B1B
C***KEYWORDS  DOUBLE PRECISION,LINEAR EQUATIONS,POSITIVE DEFINITE,
C             SYMMETRIC
C***AUTHOR  VOORHEES, E., (LANL)
C***PURPOSE  DPOFS solves a POSITIVE DEFINITE SYMMETRIC
C            double precision NXN system of linear equations.
C***DESCRIPTION
C
C    Subroutine DPOFS solves a  positive definite symmetric
C    NxN system of double precision linear equations using
C    LINPACK subroutines DPOCO and DPOSL.  That is, if A is an
C    NxN double precision positive definite symmetric matrix and if
C    X and B are double precision N-vectors, then DPOFS solves
C    the equation
C
C                          A*X=B.
C
C    The matrix A is first factored into upper and lower tri-
C    angular matrices R and R-TRANPOSE.  These factors are used to
C    find the solution vector X.  An approximate condition number is
C    calculated to provide a rough estimate of the number of
C    digits of accuracy in the computed solution.
C
C    If the equation A*X=B is to be solved for more than one vector
C    B, the factoring of A does not need to be performed again and
C    the option only to solve (ITASK .GT. 1) will be faster for
C    the succeeding solutions.  In this case, the contents of A,
C    LDA, and N must not have been altered by the user following
C    factorization (ITASK=1).  IND will not be changed by DPOFS
C    in this case.
C
C  Argument Description ***
C
C    A      DOUBLE PRECISION(LDA,N)
C             on entry, the doubly subscripted array with dimension
C               (LDA,N) which contains the coefficient matrix.  Only
C               the upper triangle, including the diagonal, of the
C               coefficient matrix need be entered and will subse-
C               quently be referenced and changed by the routine.
C             on return, A contains in its upper triangle an upper
C               triangular matrix R such that A = (R-TRANPOSE) * R .
C    LDA    INTEGER
C             the leading dimension of the array A.  LDA must be great-
C             er than or equal to N.  (terminal error message IND=-1)
C    N      INTEGER
C             the order of the matrix A.  N must be greater
C             than or equal to 1.  (terminal error message IND=-2)
C    V      DOUBLE PRECISION(N)
C             on entry, the singly subscripted array(vector) of di-
C               mension N which contains the right hand side B of a
C               system of simultaneous linear equations  A*X=B.
C             on return, V contains the solution vector, X .
C    ITASK  INTEGER
C             If ITASK = 1, the matrix A is factored and then the
C               linear equation is solved.
C             If ITASK .GT. 1, the equation is solved using the existing
C               factored matrix A.
C             If ITASK .LT. 1, then terminal error message IND=-3 is
C               printed.
C    IND    INTEGER
C             GT. 0  IND is a rough estimate of the number of digits
C                     of accuracy in the solution, X.
C             LT. 0  See error message corresponding to IND below.
C    WORK   DOUBLE PRECISION(N)
C             a singly subscripted array of dimension at least N.
C
C  Error Messages Printed ***
C
C    IND=-1  Terminal   N is greater than LDA.
C    IND=-2  Terminal   N is less than 1.
C    IND=-3  Terminal   ITASK is less than 1.
C    IND=-4  Terminal   The matrix A is computationally singular or
C                         is not positive definite.  A solution
C                         has not been computed.
C    IND=-10 Warning    The solution has no apparent significance.
C                         The solution may be inaccurate or the
C                         matrix A may be poorly scaled.
C
C               Note-  The above Terminal(*fatal*) Error Messages are
C                      designed to be handled by XERRWV in which
C                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
C                      for warning error messages from XERROR.  Unless
C                      the user provides otherwise, an error message
C                      will be printed followed by an abort.
C***REFERENCES  SUBROUTINE DPOFS WAS DEVELOPED BY GROUP C-3, LOS ALAMOS
C                 SCIENTIFIC LABORATORY, LOS ALAMOS, NM 87545.
C                 THE LINPACK SUBROUTINES USED BY DPOFS ARE DESCRIBED IN
C                 DETAIL IN THE *LINPACK USERS GUIDE* PUBLISHED BY
C                 THE SOCIETY FOR INDUSTRIAL AND APPLIED MATHEMATICS
C                 (SIAM) DATED 1979.
C***ROUTINES CALLED  D1MACH,DPOCO,DPOSL,XERROR,XERRWV
C***END PROLOGUE  DPOFS
C
      INTEGER LDA,N,ITASK,IND,INFO
      DOUBLE PRECISION A(LDA,N),V(N),WORK(N),D1MACH
      DOUBLE PRECISION RCOND
C***FIRST EXECUTABLE STATEMENT  DPOFS
      IF (LDA.LT.N)  GO TO 101
      IF (N.LE.0)  GO TO 102
      IF (ITASK.LT.1) GO TO 103
      IF (ITASK.GT.1) GO TO 20
C
C     FACTOR MATRIX A INTO R
      CALL DPOCO(A,LDA,N,RCOND,WORK,INFO)
C
C     CHECK FOR POSITIVE DEFINITE MATRIX
      IF(INFO.NE.0) GO TO 104
C
C     COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
      IND=-INT(DLOG10(D1MACH(4)/RCOND))
C
C     CHECK FOR IND GREATER THAN ZERO
      IF (IND.GT.0)  GO TO 20
      IND=-10
      CALL XERROR( 'DPOFS ERROR (IND=-10) -- SOLUTION MAY HAVE NO SIGNIF
     1ICANCE',58,-10,0)
C
C     SOLVE AFTER FACTORING
   20 CALL DPOSL(A,LDA,N,V)
      RETURN
C
C     IF LDA.LT.N, IND=-1, TERMINAL XERRWV MESSAGE
  101 IND=-1
      CALL XERRWV( 'DPOFS ERROR (IND=-1) -- LDA=I1 IS LESS THAN N=I2',
     148,-1,1,2,LDA,N,0,0,0)
      RETURN
C
C     IF N.LT.1, IND=-2, TERMINAL XERRWV MESSAGE
  102 IND=-2
      CALL XERRWV( 'DPOFS ERROR (IND=-2) -- N=I1 IS LESS THAN 1',
     143,-2,1,1,N,0,0,0,0)
      RETURN
C
C     IF ITASK.LT.1, IND=-3, TERMINAL XERRWV MESSAGE
  103 IND=-3
      CALL XERRWV( 'DPOFS ERROR (IND=-3) -- ITASK=I1 IS LESS THAN 1',
     147,-3,1,1,ITASK,0,0,0,0)
      RETURN
C
C     IF SINGULAR MATRIX, IND=-4, TERMINAL XERRWV MESSAGE
  104 IND=-4
      CALL XERRWV( 'DPOFS ERROR (IND=-4)  SINGULAR OR NOT POS DEF--NO SO
     1LUTION',58,-4,1,0,0,0,0,0,0)
      RETURN
C
      END
