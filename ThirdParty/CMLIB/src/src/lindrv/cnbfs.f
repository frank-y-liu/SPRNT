      SUBROUTINE CNBFS(ABE,LDA,N,ML,MU,V,ITASK,IND,WORK,IWORK)
C***BEGIN PROLOGUE  CNBFS
C***DATE WRITTEN   800813   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2C2
C***KEYWORDS  BAND,COMPLEX,LINEAR EQUATIONS,NONSYMMETRIC
C***AUTHOR  VOORHEES, E., (LANL)
C***PURPOSE  CNBFS solves a GENERAL NONSYMMETRIC BANDED
C            COMPLEX NXN system of linear equations.
C***DESCRIPTION
C
C    Subroutine CNBFS solves a general nonsymmetric banded NxN
C    system of single precision complex linear equations using
C    SLATEC subroutines CNBCO and CNBSL.  These are adaptations
C    of the LINPACK subroutines CGBCO and CGBSL which require
C    a different format for storing the matrix elements.  If
C    A  is an NxN complex matrix and if  X  and  B  are complex
C    N-vectors, then CNBFS solves the equation
C
C                          A*X=B.
C
C    A band matrix is a matrix whose nonzero elements are all
C    fairly near the main diagonal, specifically  A(I,J) = 0
C    if  I-J is greater than  ML  or  J-I  is greater than
C    MU .  The integers ML and MU are called the lower and upper
C    band widths and  M = ML+MU+1  is the total band width.
C    CNBFS uses less time and storage than the corresponding
C    program for general matrices (CGEFS) if 2*ML+MU .LT. N .
C
C    The matrix A is first factored into upper and lower tri-
C    angular matrices U and L using partial pivoting.  These
C    factors and the pivoting information are used to find the
C    solution vector X.  An approximate condition number is
C    calculated to provide a rough estimate of the number of
C    digits of accuracy in the computed solution.
C
C    If the equation A*X=B is to be solved for more than one vector
C    B, the factoring of A does not need to be performed again and
C    the option to only solve (ITASK .GT. 1) will be faster for
C    the succeeding solutions.  In this case, the contents of A,
C    LDA, N and IWORK must not have been altered by the user follow-
C    ing factorization (ITASK=1).  IND will not be changed by CNBFS
C    in this case.
C
C
C    Band Storage
C
C          If  A  is a band matrix, the following program segment
C          will set up the input.
C
C                  ML = (band width below the diagonal)
C                  MU = (band width above the diagonal)
C                  DO 20 I = 1, N
C                     J1 = MAX0(1, I-ML)
C                     J2 = MIN0(N, I+MU)
C                     DO 10 J = J1, J2
C                        K = J - I + ML + 1
C                        ABE(I,K) = A(I,J)
C               10    CONTINUE
C               20 CONTINUE
C
C          This uses columns  1  through  ML+MU+1  of ABE .
C          Furthermore,  ML  additional columns are needed in
C          ABE  starting with column  ML+MU+2  for elements
C          generated during the triangularization.  The total
C          number of columns needed in  ABE  is  2*ML+MU+1 .
C
C    Example:  If the original matrix is
C
C          11 12 13  0  0  0
C          21 22 23 24  0  0
C           0 32 33 34 35  0
C           0  0 43 44 45 46
C           0  0  0 54 55 56
C           0  0  0  0 65 66
C
C     then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABE should contain
C
C           * 11 12 13  +     , * = not used
C          21 22 23 24  +     , + = used for pivoting
C          32 33 34 35  +
C          43 44 45 46  +
C          54 55 56  *  +
C          65 66  *  *  +
C
C
C  Argument Description ***
C
C    ABE    COMPLEX(LDA,NC)
C             on entry, contains the matrix in band storage as
C               described above.  NC  must not be less than
C               2*ML+MU+1 .  The user is cautioned to specify  NC
C               with care since it is not an argument and cannot
C               be checked by CNBFS.  The rows of the original
C               matrix are stored in the rows of  ABE  and the
C               diagonals of the original matrix are stored in
C               columns  1  through  ML+MU+1  of  ABE .
C             on return, contains an upper triangular matrix U and
C               the multipliers necessary to construct a matrix L
C               so that A=L*U.
C    LDA    INTEGER
C             the leading dimension of array ABE.  LDA must be great-
C             er than or equal to N.  (terminal error message IND=-1)
C    N      INTEGER
C             the order of the matrix A.  N must be greater
C             than or equal to 1 .  (terminal error message IND=-2)
C    ML     INTEGER
C             the number of diagonals below the main diagonal.
C             ML  must not be less than zero nor greater than or
C             equal to  N .  (terminal error message IND=-5)
C    MU     INTEGER
C             the number of diagonals above the main diagonal.
C             MU  must not be less than zero nor greater than or
C             equal to  N .  (terminal error message IND=-6)
C    V      COMPLEX(N)
C             on entry, the singly subscripted array(vector) of di-
C               mension N which contains the right hand side B of a
C               system of simultaneous linear equations A*X=B.
C             on return, V contains the solution vector, X .
C    ITASK  INTEGER
C             if ITASK = 1, the matrix A is factored and then the
C               linear equation is solved.
C             if ITASK .GT. 1, the equation is solved using the existing
C               factored matrix A and IWORK.
C             if ITASK .LT. 1, then terminal error message IND=-3 is
C               printed.
C    IND    INTEGER
C             GT. 0  IND is a rough estimate of the number of digits
C                     of accuracy in the solution, X.
C             LT. 0  see error message corresponding to IND below.
C    WORK   COMPLEX(N)
C             a singly subscripted array of dimension at least N.
C    IWORK  INTEGER(N)
C             a singly subscripted array of dimension at least N.
C
C  Error Messages Printed ***
C
C    IND=-1  terminal   N is greater than LDA.
C    IND=-2  terminal   N is less than 1.
C    IND=-3  terminal   ITASK is less than 1.
C    IND=-4  terminal   The matrix A is computationally singular.
C                         A solution has not been computed.
C    IND=-5  terminal   ML is less than zero or is greater than
C                         or equal to N .
C    IND=-6  terminal   MU is less than zero or is greater than
C                         or equal to N .
C    IND=-10 warning    The solution has no apparent significance.
C                         The solution may be inaccurate or the matrix
C                         A may be poorly scaled.
C
C               NOTE-  The above terminal(*fatal*) error messages are
C                      designed to be handled by XERRWV in which
C                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
C                      for warning error messages from XERROR.  Unless
C                      the user provides otherwise, an error message
C                      will be printed followed by an abort.
C***REFERENCES  SUBROUTINE CNBFS WAS DEVELOPED BY GROUP C-3, LOS ALAMOS
C                 SCIENTIFIC LABORATORY, LOS ALAMOS, NM 87545.
C                 THE SLATEC SUBROUTINES USED BY CNBFS ARE SIMILAR
C                 (EXCEPT FOR MATRIX STORAGE FORMAT) TO THE CORRESPON-
C                 DING LINPACK SUBROUTINES (CGBCO AND CGBSL) DESCRIBED
C                 IN DETAIL IN THE *LINPACK USERS GUIDE* PUBLISHED BY
C                 THE SOCIETY FOR INDUSTRIAL AND APPLIED MATHEMATICS
C                 (SIAM) DATED 1979.
C***ROUTINES CALLED  CNBCO,CNBSL,R1MACH,XERROR,XERRWV
C***END PROLOGUE  CNBFS
C
      INTEGER LDA,N,ITASK,IND,IWORK(N),ML,MU
      COMPLEX ABE(LDA,*),V(N),WORK(N)
      REAL RCOND
      REAL R1MACH
C***FIRST EXECUTABLE STATEMENT  CNBFS
      IF (LDA.LT.N)  GO TO 101
      IF (N.LE.0)  GO TO 102
      IF (ITASK.LT.1) GO TO 103
      IF((ML.LT.0).OR.(ML.GE.N))GO TO 105
      IF((MU.LT.0).OR.(MU.GE.N))GO TO 106
      IF (ITASK.GT.1) GO TO 20
C
C     FACTOR MATRIX A INTO LU
      CALL CNBCO(ABE,LDA,N,ML,MU,IWORK,RCOND,WORK)
C
C     CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
      IF (RCOND.EQ.0.0)  GO TO 104
C
C     COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
      IND=-INT(ALOG10(R1MACH(4)/RCOND))
C
C     CHECK FOR IND GREATER THAN ZERO
      IF (IND.GT.0)  GO TO 20
      IND=-10
      CALL XERROR( 'CNBFS ERROR (IND=-10) -- SOLUTION MAY HAVE NO SIGNIF
     1ICANCE',58,-10,0)
C
C     SOLVE AFTER FACTORING
   20 CALL CNBSL(ABE,LDA,N,ML,MU,IWORK,V,0)
      RETURN
C
C     IF LDA.LT.N, IND=-1, TERMINAL XERRWV MESSAGE
  101 IND=-1
      CALL XERRWV( 'CNBFS ERROR (IND=-1) -- LDA=I1 IS LESS THAN N=I2',
     148,-1,1,2,LDA,N,0,0,0)
      RETURN
C
C     IF N.LT.1, IND=-2, TERMINAL XERRWV MESSAGE
  102 IND=-2
      CALL XERRWV( 'CNBFS ERROR (IND=-2) -- N=I1 IS LESS THAN 1',
     143,-2,1,1,N,0,0,0,0)
      RETURN
C
C     IF ITASK.LT.1, IND=-3, TERMINAL XERRWV MESSAGE
  103 IND=-3
      CALL XERRWV( 'CNBFS ERROR (IND=-3) -- ITASK=I1 IS LESS THAN 1',
     147,-3,1,1,ITASK,0,0,0,0)
      RETURN
C
C     IF SINGULAR MATRIX, IND=-4, TERMINAL XERRWV MESSAGE
  104 IND=-4
      CALL XERRWV( 'CNBFS ERROR (IND=-4) -- SINGULAR MATRIX A - NO SOLUT
     1ION',55,-4,1,0,0,0,0,0,0)
      RETURN
C
C     IF IMPROPER ML VALUE, IND=-5, TERMINAL XERRWV MESSAGE
  105 IND=-5
      CALL XERRWV( 'CNBFS ERROR (IND=-5) -- ML=I1 IS OUT OF RANGE',45,-5
     1,1,1,ML,0,0,0,0)
      RETURN
C
C     IF IMPROPER MU VALUE, IND=-6, TERMINAL XERRWV MESSAGE
  106 IND=-6
      CALL XERRWV( 'CNBFS ERROR (IND=-6) -- MU=I1 IS OUT OF RANGE',45,-6
     1,1,1,MU,0,0,0,0)
      RETURN
C
      END
