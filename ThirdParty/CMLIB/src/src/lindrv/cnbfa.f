      SUBROUTINE CNBFA(ABE,LDA,N,ML,MU,IPVT,INFO)
C***BEGIN PROLOGUE  CNBFA
C***DATE WRITTEN   800730   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2C2
C***KEYWORDS  BAND,COMPLEX,LINEAR EQUATIONS,NONSYMMETRIC
C***AUTHOR  VOORHEES, E., (LANL)
C***PURPOSE  Factors a COMPLEX BAND matrix by elimination.
C***DESCRIPTION
C
C    CNBFA factors a complex band matrix by elimination.
C
C     CNBFA is usually called by CNBCO, but it can be called
C     directly with a saving in time if RCOND is not needed.
C
C     On Entry
C
C        ABE     COMPLEX(LDA, NC)
C                contains the matrix in band storage.  The rows
C                of the original matrix are stored in the rows
C                of ABE and the diagonals of the original matrix
C                are stored in columns 1 through ML+MU+1 of ABE.
C                NC must be .GE. 2*ML+MU+1 .
C                See the comments below for details.
C
C        LDA     INTEGER
C                the leading dimension of the array ABE.
C                LDA must be .GE. N .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C                0 .LE. ML .LT. N .
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C                0 .LE. MU .LT. N .
C                More efficient if ML .LE. MU .
C
C     On Return
C
C        ABE     an upper triangular matrix in band storage
C                and the multipliers which were used to obtain it.
C                the factorization can be written  A = L*U  where
C                L is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                =0  normal value
C                =K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                condition for this subroutine, but it does
C                indicate that CNBSL will divide by zero if
C                called.  Use RCOND in CNBCO for a reliable
C                indication of singularity.
C
C     Band Storage
C
C           If  A  is a band matrix, the following program segment
C           will set up the input.
C
C                   ML = (band width below the diagonal)
C                   MU = (band width above the diagonal)
C                   DO 20 I = 1, N
C                      J1 = MAX0(1, I-ML)
C                      J2 = MIN0(N, I+MU)
C                      DO 10 J = J1, J2
C                         K = J - I + ML + 1
C                         ABE(I,K) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           This uses columns  1  through  ML+MU+1  of ABE .
C           Furthermore,  ML  additional columns are needed in
C           ABE  starting with column  ML+MU+2  for elements
C           generated during the triangularization.  The total
C           number of columns needed in  ABE  is  2*ML+MU+1 .
C
C     Example:  If the original matrix is
C
C           11 12 13  0  0  0
C           21 22 23 24  0  0
C            0 32 33 34 35  0
C            0  0 43 44 45 46
C            0  0  0 54 55 56
C            0  0  0  0 65 66
C
C      then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABE should contain
C
C            * 11 12 13  +     , * = not used
C           21 22 23 24  +     , + = used for pivoting
C           32 33 34 35  +
C           43 44 45 46  +
C           54 55 56  *  +
C           65 66  *  *  +
C
C     SLATEC.  This version dated 07/30/80 .
C     E. A. Voorhees, Los Alamos Scientific Laboratory
C
C     Subroutines and Functions
C
C     BLAS ICAMAX,CAXPY,CSCAL,CSWAP
C     Fortran  ABS,AIMAG,MIN0,REAL
C***REFERENCES  SUBROUTINE CNBFA WAS DEVELOPED BY GROUP C-3, LOS ALAMOS
C                 SCIENTIFIC LABORATORY, LOS ALAMOS, NM 87545.
C***ROUTINES CALLED  CAXPY,CSCAL,CSWAP,ICAMAX
C***END PROLOGUE  CNBFA
      INTEGER LDA,N,ML,MU,IPVT(*),INFO
      COMPLEX ABE(LDA,*)
C
      INTEGER ML1,MB,M,N1,LDB,I,J,K,L,LM,LM1,LM2,MP,ICAMAX
      COMPLEX T
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
C***FIRST EXECUTABLE STATEMENT  CNBFA
      ML1=ML+1
      MB=ML+MU
      M=ML+MU+1
      N1=N-1
      LDB=LDA-1
      INFO=0
C
C     SET FILL-IN COLUMNS TO ZERO
C
      IF(N.LE.1)GO TO 50
      IF(ML.LE.0)GO TO 7
      DO 6 J=1,ML
        DO 5 I=1,N
          ABE(I,M+J)=(0.0E0,0.0E0)
    5   CONTINUE
    6 CONTINUE
    7 CONTINUE
C
C     GAUSSIAN ELIMINATION WITH PARTIAL ELIMINATION
C
      DO 40 K=1,N1
        LM=MIN0(N-K,ML)
        LM1=LM+1
        LM2=ML1-LM
C
C     SEARCH FOR PIVOT INDEX
C
        L=-ICAMAX(LM1,ABE(LM+K,LM2),LDB)+LM1+K
        IPVT(K)=L
        MP=MIN0(MB,N-K)
C
C     SWAP ROWS IF NECESSARY
C
        IF(L.NE.K)CALL CSWAP(MP+1,ABE(K,ML1),LDA,ABE(L,ML1+K-L),LDA)
C
C     SKIP COLUMN REDUCTION IF PIVOT IS ZERO
C
        IF(CABS1(ABE(K,ML1)).EQ.0.0E0) GO TO 20
C
C     COMPUTE MULTIPLIERS
C
        T=-(1.0E0,0.0E0)/ABE(K,ML1)
        CALL CSCAL(LM,T,ABE(LM+K,LM2),LDB)
C
C     ROW ELIMINATION WITH COLUMN INDEXING
C
        DO 10 J=1,MP
          CALL CAXPY(LM,ABE(K,ML1+J),ABE(LM+K,LM2),LDB,ABE(LM+K,LM2+J),L
     1DB)
   10   CONTINUE
        GO TO 30
   20   CONTINUE
        INFO=K
   30   CONTINUE
   40 CONTINUE
   50 CONTINUE
      IPVT(N)=N
      IF(CABS1(ABE(N,ML1)).EQ.0.0E0) INFO=N
      RETURN
      END
