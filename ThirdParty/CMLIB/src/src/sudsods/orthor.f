      SUBROUTINE ORTHOR(A,N,M,NRDA,IFLAG,IRANK,ISCALE,DIAG,KPIVOT,
     1   SCALES,ROWS,RS)
C***BEGIN PROLOGUE  ORTHOR
C***REFER TO  SUDS
C
C   Reduction of the matrix A to lower triangular form by a sequence of
C   orthogonal HOUSEHOLDER transformations post-multiplying A
C
C     Written by H. A. Watts, SANDIA LABORATORIES, ALBUQUERQUE
C     and modeled after the ALGOL codes from the following
C
C    References
C     G. Golub, Numerical Methods For Solving Linear Least Squares
C     Problems, NUMER. MATH.,V.7,PP.206-216,1965.
C     P. Businger and G. Golub, Linear Least Squares Solutions By
C     Householder Transformations, NUMER. MATH.,V.7,PP.269-276,1965.
C
C **********************************************************************
C   INPUT
C **********************************************************************
C
C     A -- Contains the matrix to be decomposed, must be dimensioned
C           NRDA by N
C     N -- Number of rows in the matrix, N greater or equal to 1
C     M -- Number of columns in the matrix, M greater or equal to N
C     IFLAG -- Indicates the uncertainty in the matrix data
C             = 0 when the data is to be treated as exact
C             =-K when the data is assumed to be accurate to about
C                 K digits
C     ISCALE -- Scaling indicator
C               =-1 if the matrix is to be pre-scaled by
C               columns when appropriate.
C               Otherwise no scaling will be attempted
C     NRDA -- Row dimension of A, NRDA greater or equal to N
C     DIAG,KPIVOT,ROWS -- Arrays of length at least N used internally
C         ,RS,SCALES         (except for SCALES which is M)
C
C **********************************************************************
C   OUTPUT
C **********************************************************************
C
C     IFLAG - status indicator
C            =1 for successful decomposition
C            =2 if improper input is detected
C            =3 if rank of the matrix is less than N
C     A -- contains the reduced matrix in the strictly lower triangular
C          part and transformation information
C     IRANK -- contains the numerically determined matrix rank
C     DIAG -- contains the diagonal elements of the reduced
C             triangular matrix
C     KPIVOT -- contains the pivotal information,the column interchanges
C               performed on the original matrix are recorded here
C     SCALES -- contains the column scaling parameters
C
C **********************************************************************
C***ROUTINES CALLED  CSCALE,R1MACH,SDOT,XERROR
C***END PROLOGUE  ORTHOR
      DIMENSION A(NRDA,N),DIAG(N),KPIVOT(N),ROWS(N),RS(N),SCALES(M)
C
C END OF ABSTRACT
C
C **********************************************************************
C
C     MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED
C     BY THE FUNCTION R1MACH.
C
C **********************************************************************
C
C***FIRST EXECUTABLE STATEMENT  ORTHOR
      URO = R1MACH(4)
      IF (M .GE. N  .AND.  N .GE. 1  .AND.  NRDA .GE. N) GO TO 1
      IFLAG=2
      CALL XERROR ( 'ORTHOR-INVALID INPUT PARAMETERS.',32,2,1)
      RETURN
C
    1 ACC=10.*URO
      IF (IFLAG .LT. 0) ACC=AMAX1(ACC,10.**IFLAG)
      SRURO=SQRT(URO)
      IFLAG=1
      IRANK=N
C
C     COMPUTE NORM**2 OF JTH ROW AND A MATRIX NORM
C
      ANORM=0.
      DO 2 J=1,N
         KPIVOT(J)=J
         ROWS(J)=SDOT(M,A(J,1),NRDA,A(J,1),NRDA)
         RS(J)=ROWS(J)
         ANORM=ANORM+ROWS(J)
    2 CONTINUE
C
C     PERFORM COLUMN SCALING ON A WHEN SPECIFIED
C
      CALL CSCALE(A,NRDA,N,M,SCALES,DUM,ROWS,RS,ANORM,SCALES,ISCALE,1)
C
      ANORM=SQRT(ANORM)
C
C
C     CONSTRUCTION OF LOWER TRIANGULAR MATRIX AND RECORDING OF
C     ORTHOGONAL TRANSFORMATIONS
C
C
      DO 50 K=1,N
         MK=M-K+1
         IF (K .EQ. N) GO TO 25
         KP=K+1
C
C        SEARCHING FOR PIVOTAL ROW
C
         DO 10 J=K,N
            IF (ROWS(J) .GE. SRURO*RS(J)) GO TO 5
            ROWS(J)=SDOT(MK,A(J,K),NRDA,A(J,K),NRDA)
            RS(J)=ROWS(J)
    5       IF (J .EQ. K) GO TO 7
            IF (SIGMA .GE. 0.99*ROWS(J)) GO TO 10
    7       SIGMA=ROWS(J)
            JROW=J
   10    CONTINUE
         IF (JROW .EQ. K) GO TO 25
C
C        PERFORM ROW INTERCHANGE
C
         L=KPIVOT(K)
         KPIVOT(K)=KPIVOT(JROW)
         KPIVOT(JROW)=L
         ROWS(JROW)=ROWS(K)
         ROWS(K)=SIGMA
         RSS=RS(K)
         RS(K)=RS(JROW)
         RS(JROW)=RSS
         DO 20 L=1,M
            ASAVE=A(K,L)
            A(K,L)=A(JROW,L)
   20       A(JROW,L)=ASAVE
C
C        CHECK RANK OF THE MATRIX
C
   25    SIG=SDOT(MK,A(K,K),NRDA,A(K,K),NRDA)
         DIAGK=SQRT(SIG)
         IF (DIAGK .GT. ACC*ANORM) GO TO 30
C
C        RANK DEFICIENT PROBLEM
         IFLAG=3
         IRANK=K-1
         CALL XERROR(  'ORTHOR-RANK OF MATRIX IS LESS THAN THE NUMBER OF
     1 ROWS.',54,1,1)
         RETURN
C
C        CONSTRUCT AND APPLY TRANSFORMATION TO MATRIX A
C
   30    AKK=A(K,K)
         IF (AKK .GT. 0.) DIAGK=-DIAGK
         DIAG(K)=DIAGK
         A(K,K)=AKK-DIAGK
         IF (K .EQ. N) GO TO 50
         SAD=DIAGK*AKK-SIG
         DO 40 J=KP,N
            AS=SDOT(MK,A(K,K),NRDA,A(J,K),NRDA)/SAD
            DO 35 L=K,M
   35          A(J,L)=A(J,L)+AS*A(K,L)
   40       ROWS(J)=ROWS(J)-A(J,K)**2
   50 CONTINUE
C
C
      RETURN
      END
