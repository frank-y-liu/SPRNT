      SUBROUTINE SNBCO(ABE,LDA,N,ML,MU,IPVT,RCOND,Z)
C***BEGIN PROLOGUE  SNBCO
C***DATE WRITTEN   800723   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2A2
C***KEYWORDS  BAND,LINEAR EQUATIONS,NONSYMMETRIC
C***AUTHOR  VOORHEES, E., (LANL)
C***PURPOSE  SNBCO factors a real BAND matrix by Gaussian
C            elimination and estimates the condition number.
C***DESCRIPTION
C
C     SNBCO factors a real band matrix by Gaussian
C     elimination and estimates the condition of the matrix.
C
C     If RCOND is not needed, SNBFA is slightly faster.
C     To solve  A*X = B , follow SNBCO by SNBSL.
C     To compute  INVERSE(A)*C , follow SNBCO by SNBSL.
C     To compute  DETERMINANT(A) , follow SNBCO by SNBDI.
C
C     On Entry
C
C        ABE     REAL(LDA, NC)
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
C                The factorization can be written  A = L*U , where
C                L is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        RCOND   REAL
C                an estimate of the reciprocal condition of  A .
C                For the system  A*X = B , relative perturbations
C                in  A  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                         1.0 + RCOND .EQ. 1.0
C                is true, then  A  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C
C        Z       REAL(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
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
C     SLATEC.  This version dated 07/23/80 .
C     E. A. Voorhees, Los Alamos Scientific Laboratory
C
C     Subroutines and Functions
C
C     SLATEC SNBFA
C      SDOT,SAXPY,SSCAL,SASUM
C     Fortran  ABS,AMAX1,MAX0,MIN0,SIGN
C***REFERENCES  SUBROUTINE SNBCO WAS DEVELOPED BY GROUP C-3, LOS ALAMOS
C                 SCIENTIFIC LABORATORY, LOS ALAMOS, NM 87545.
C***ROUTINES CALLED  SASUM,SAXPY,SDOT,SNBFA,SSCAL
C***END PROLOGUE  SNBCO
      INTEGER LDA,N,ML,MU,IPVT(*)
      REAL ABE(LDA,*),Z(*)
      REAL RCOND
C
      REAL SDOT,EK,T,WK,WKM
      REAL ANORM,S,SASUM,SM,YNORM
      INTEGER I,INFO,J,JU,K,KB,KP1,L,LDB,LM,LZ,M,ML1,MM,NL,NU
C***FIRST EXECUTABLE STATEMENT  SNBCO
      ML1=ML+1
      LDB = LDA - 1
      ANORM = 0.0E0
      DO 10 J = 1, N
        NU = MIN0(MU,J-1)
        NL = MIN0(ML,N-J)
        L = 1 + NU + NL
        ANORM = AMAX1(ANORM,SASUM(L,ABE(J+NL,ML1-NL),LDB))
   10 CONTINUE
C
C     FACTOR
C
      CALL SNBFA(ABE,LDA,N,ML,MU,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF  W WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0E0
      DO 20 J = 1, N
        Z(J) = 0.0E0
   20 CONTINUE
      M = ML + MU + 1
      JU = 0
      DO 100 K = 1, N
        IF (Z(K) .NE. 0.0E0) EK = SIGN(EK,-Z(K))
        IF (ABS(EK-Z(K)) .LE. ABS(ABE(K,ML1))) GO TO 30
          S = ABS(ABE(K,ML1))/ABS(EK-Z(K))
          CALL SSCAL(N,S,Z,1)
          EK = S*EK
   30   CONTINUE
        WK = EK - Z(K)
        WKM = -EK - Z(K)
        S = ABS(WK)
        SM = ABS(WKM)
        IF (ABE(K,ML1) .EQ. 0.0E0) GO TO 40
          WK = WK/ABE(K,ML1)
          WKM = WKM/ABE(K,ML1)
        GO TO 50
   40   CONTINUE
          WK = 1.0E0
          WKM = 1.0E0
   50   CONTINUE
        KP1 = K + 1
        JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
        MM = ML1
        IF (KP1 .GT. JU) GO TO 90
          DO 60 I = KP1, JU
            MM = MM + 1
            SM = SM + ABS(Z(I)+WKM*ABE(K,MM))
            Z(I) = Z(I) + WK*ABE(K,MM)
            S = S + ABS(Z(I))
   60     CONTINUE
          IF (S .GE. SM) GO TO 80
            T = WKM -WK
            WK = WKM
            MM = ML1
            DO 70 I = KP1, JU
              MM = MM + 1
              Z(I) = Z(I) + T*ABE(K,MM)
   70       CONTINUE
   80     CONTINUE
   90   CONTINUE
      Z(K) = WK
  100 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
        K = N + 1 - KB
        NL = MIN0(ML,N-K)
        IF (K .LT. N) Z(K) = Z(K) + SDOT(NL,ABE(K+NL,ML1-NL),-LDB,Z(K+1)
     1  ,1)
        IF (ABS(Z(K)) .LE. 1.0E0) GO TO 110
          S = 1.0E0/ABS(Z(K))
          CALL SSCAL(N,S,Z,1)
  110   CONTINUE
        L = IPVT(K)
        T = Z(L)
        Z(L) = Z(K)
        Z(K) = T
  120 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
C
      YNORM = 1.0E0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
        L = IPVT(K)
        T = Z(L)
        Z(L) = Z(K)
        Z(K) = T
        NL = MIN0(ML,N-K)
        IF (K .LT. N) CALL SAXPY(NL,T,ABE(K+NL,ML1-NL),-LDB,Z(K+1),1)
        IF (ABS(Z(K)) .LE. 1.0E0) GO TO 130
          S = 1.0E0/ABS(Z(K))
          CALL SSCAL(N,S,Z,1)
          YNORM = S*YNORM
  130   CONTINUE
  140 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
        K = N + 1 - KB
        IF (ABS(Z(K)) .LE. ABS(ABE(K,ML1))) GO TO 150
          S = ABS(ABE(K,ML1))/ABS(Z(K))
          CALL SSCAL(N,S,Z,1)
          YNORM = S*YNORM
  150   CONTINUE
        IF (ABE(K,ML1) .NE. 0.0E0) Z(K) = Z(K)/ABE(K,ML1)
        IF (ABE(K,ML1) .EQ. 0.0E0) Z(K) = 1.0E0
        LM = MIN0(K,M) - 1
        LZ = K - LM
        T = -Z(K)
        CALL SAXPY(LM,T,ABE(K-1,ML+2),-LDB,Z(LZ),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0E0
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
