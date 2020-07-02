      SUBROUTINE DPPSL(AP,N,B)
C***BEGIN PROLOGUE  DPPSL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2B1B
C***KEYWORDS  DOUBLE PRECISION,LINEAR ALGEBRA,LINPACK,MATRIX,PACKED,
C             POSITIVE DEFINITE,SOLVE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Solves the double precision SYMMETRIC POSITIVE DEFINITE
C            system A*X=B
C***DESCRIPTION
C
C     DPPSL solves the double precision symmetric positive definite
C     system A * X = B
C     using the factors computed by DPPCO or DPPFA.
C
C     On Entry
C
C        AP      DOUBLE PRECISION (N*(N+1)/2)
C                the output from DPPCO or DPPFA.
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal.  Technically this indicates
C        singularity, but it is usually caused by improper subroutine
C        arguments.  It will not occur if the subroutines are called
C        correctly and  INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DPPCO(AP,N,RCOND,Z,INFO)
C           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL DPPSL(AP,N,C(1,J))
C        10 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DDOT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DDOT
C***END PROLOGUE  DPPSL
      INTEGER N
      DOUBLE PRECISION AP(*),B(*)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,KK
C***FIRST EXECUTABLE STATEMENT  DPPSL
      KK = 0
      DO 10 K = 1, N
         T = DDOT(K-1,AP(KK+1),1,B(1),1)
         KK = KK + K
         B(K) = (B(K) - T)/AP(KK)
   10 CONTINUE
      DO 20 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/AP(KK)
         KK = KK - K
         T = -B(K)
         CALL DAXPY(K-1,T,AP(KK+1),1,B(1),1)
   20 CONTINUE
      RETURN
      END
