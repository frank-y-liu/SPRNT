      SUBROUTINE SPPSL(AP,N,B)
C***BEGIN PROLOGUE  SPPSL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2B1B
C***KEYWORDS  LINEAR ALGEBRA,LINPACK,MATRIX,PACKED,POSITIVE DEFINITE,
C             SOLVE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Solves the real SYMMETRIC POSITIVE DEFINITE system A*X=B
C            using the factors computed by SPPCO or SPPFA.
C***DESCRIPTION
C
C     SPPSL solves the real symmetric positive definite system
C     A * X = B
C     using the factors computed by SPPCO or SPPFA.
C
C     On Entry
C
C        AP      REAL (N*(N+1)/2)
C                the output from SPPCO or SPPFA.
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        B       REAL(N)
C                the right hand side vector.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal.  Technically, this indicates
C        singularity, but it is usually caused by improper subroutine
C        arguments.  It will not occur if the subroutines are called
C        correctly and  INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL SPPCO(AP,N,RCOND,Z,INFO)
C           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL SPPSL(AP,N,C(1,J))
C        10 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS SAXPY,SDOT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  SAXPY,SDOT
C***END PROLOGUE  SPPSL
      INTEGER N
      REAL AP(*),B(*)
C
      REAL SDOT,T
      INTEGER K,KB,KK
C***FIRST EXECUTABLE STATEMENT  SPPSL
      KK = 0
      DO 10 K = 1, N
         T = SDOT(K-1,AP(KK+1),1,B(1),1)
         KK = KK + K
         B(K) = (B(K) - T)/AP(KK)
   10 CONTINUE
      DO 20 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/AP(KK)
         KK = KK - K
         T = -B(K)
         CALL SAXPY(K-1,T,AP(KK+1),1,B(1),1)
   20 CONTINUE
      RETURN
      END
