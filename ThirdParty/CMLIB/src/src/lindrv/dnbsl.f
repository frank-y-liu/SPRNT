      SUBROUTINE DNBSL(ABE,LDA,N,ML,MU,IPVT,B,JOB)
C***BEGIN PROLOGUE  DNBSL
C***DATE WRITTEN   800728   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2A2
C***KEYWORDS  BAND,DOUBLE PRECISION,LINEAR EQUATIONS,NONSYMMETRIC
C***AUTHOR  VOORHEES, E., (LANL)
C***PURPOSE  Solves a double precision BAND system using
C            the factors computed by DNBCO or DNBFA.
C***DESCRIPTION
C
C     DNBSL solves the double precision band system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DNBCO or DNBFA.
C
C     On Entry
C
C        ABE     DOUBLE PRECISION(LDA, NC)
C                the output from DNBCO or DNBFA.
C                NC must be .GE. 2*ML+MU+1 .
C
C        LDA     INTEGER
C                the leading dimension of the array  ABE .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C
C        IPVT    INTEGER(N)
C                the pivot vector from DNBCO or DNBFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B .
C                = nonzero   to solve  TRANS(A)*X = B , where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA.  It will not occur if the subroutines are
C        called correctly and if DNBCO has set RCOND .GT. 0.0
C        or DNBFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DNBCO(ABE,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C             CALL DNBSL(ABE,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     SLATEC.  This version dated 07/28/80 .
C     E. A. Voorhees, Los Alamos Scientific Laboratory
C
C     Subroutines and Functions
C
C      DAXPY,DDOT
C     Fortran  MIN0
C***REFERENCES  SUBROUTINE DNBSL WAS DEVELOPED BY GROUP C-3, LOS ALAMOS
C                 SCIENTIFIC LABORATORY, LOS ALAMOS, NM 87545.
C***ROUTINES CALLED  DAXPY,DDOT
C***END PROLOGUE  DNBSL
      INTEGER LDA,N,ML,MU,IPVT(*),JOB
      DOUBLE PRECISION ABE(LDA,*),B(*)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,LB,LDB,LM,M,MLM,NM1
C***FIRST EXECUTABLE STATEMENT  DNBSL
      M=MU+ML+1
      NM1=N-1
      LDB=1-LDA
      IF(JOB.NE.0)GO TO 50
C
C       JOB = 0 , SOLVE  A * X = B
C       FIRST SOLVE L*Y = B
C
        IF(ML.EQ.0)GO TO 30
        IF(NM1.LT.1)GO TO 30
          DO 20 K=1,NM1
            LM=MIN0(ML,N-K)
            L=IPVT(K)
            T=B(L)
            IF(L.EQ.K)GO TO 10
              B(L)=B(K)
              B(K)=T
   10       CONTINUE
            MLM=ML-(LM-1)
            CALL DAXPY(LM,T,ABE(K+LM,MLM),LDB,B(K+1),1)
   20     CONTINUE
   30   CONTINUE
C
C       NOW SOLVE  U*X = Y
C
        DO 40 KB=1,N
          K=N+1-KB
          B(K)=B(K)/ABE(K,ML+1)
          LM=MIN0(K,M)-1
          LB=K-LM
          T=-B(K)
          CALL DAXPY(LM,T,ABE(K-1,ML+2),LDB,B(LB),1)
   40   CONTINUE
      GO TO 100
   50 CONTINUE
C
C       JOB = NONZERO, SOLVE TRANS(A) * X = B
C       FIRST SOLVE  TRANS(U)*Y = B
C
        DO 60 K = 1, N
          LM = MIN0(K,M) - 1
          LB = K - LM
          T = DDOT(LM,ABE(K-1,ML+2),LDB,B(LB),1)
          B(K) = (B(K) - T)/ABE(K,ML+1)
   60   CONTINUE
C
C       NOW SOLVE TRANS(L)*X = Y
C
        IF (ML .EQ. 0) GO TO 90
        IF (NM1 .LT. 1) GO TO 90
          DO 80 KB = 1, NM1
            K = N - KB
            LM = MIN0(ML,N-K)
            MLM = ML - (LM - 1)
            B(K) = B(K) + DDOT(LM,ABE(K+LM,MLM),LDB,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
              T = B(L)
              B(L) = B(K)
              B(K) = T
   70       CONTINUE
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
      RETURN
      END
