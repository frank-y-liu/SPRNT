      SUBROUTINE CHIEV(A,LDA,N,E,V,LDV,WORK,JOB,INFO)
C***BEGIN PROLOGUE  CHIEV
C***DATE WRITTEN   800808   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D4A3
C***KEYWORDS  COMPLEX,COMPLEX HERMITIAN,EIGENANALYSIS,EIGENVALUE,
C             EIGENVECTOR,MATRIX
C***AUTHOR  KAHANER, D. K., (NBS)
C           MOLER, C. B., (U. OF NEW MEXICO)
C           STEWART, G. W., (U. OF MARYLAND)
C***PURPOSE  To compute the eigenvalues and, optionally, the eigen-
C            vectors of a COMPLEX HERMITIAN matrix.
C***DESCRIPTION
C
C     David Kahaner, Cleve Moler, G. W. Stewart,
C       N.B.S.         U.N.M.      N.B.S./U.MD.
C
C     Abstract
C      CHIEV computes the eigenvalues and, optionally,
C      the eigenvectors of a complex Hermitian matrix.
C
C     Call Sequence Parameters-
C       (the values of parameters marked with * (star) will be changed
C         by CHIEV.)
C
C        A*      COMPLEX(LDA,N)
C                complex Hermitian input matrix.
C                Only the upper triangle of A need be
C                filled in.  Elements on diagonal must be real.
C
C        LDA     INTEGER
C                set by the user to
C                the leading dimension of the complex array A.
C
C        N       INTEGER
C                set by the user to
C                the order of the matrices A and V, and
C                the number of elements in E.
C
C        E*      REAL(N)
C                on return from CHIEV E contains the eigenvalues of A.
C                See also INFO below.
C
C        V*      COMPLEX(LDV,N)
C                on return from CHIEV if the user has set JOB
C                = 0        V is not referenced.
C                = nonzero  the N eigenvectors of A are stored in the
C                first N columns of V.  See also INFO below.
C
C        LDV     INTEGER
C                set by the user to
C                the leading dimension of the array V if JOB is also
C                set nonzero.  In that case N must be .LE. LDV.
C                If JOB is set to zero LDV is not referenced.
C
C        WORK*   REAL(4N)
C                temporary storage vector.  Contents changed by CHIEV.
C
C        JOB     INTEGER
C                set by the user to
C                = 0        eigenvalues only to be calculated by CHIEV.
C                           Neither V nor LDV are referenced.
C                = nonzero  eigenvalues and vectors to be calculated.
C                           In this case A and V must be distinct arrays
C                           also if LDA .GT. LDV CHIEV changes all the
C                           elements of A thru column N.  If LDA < LDV
C                           CHIEV changes all the elements of V through
C                           column N.  If LDA = LDV only A(I,J) and V(I,
C                           J) for I,J = 1,...,N are changed by CHIEV.
C
C        INFO*   INTEGER
C                on return from CHIEV the value of INFO is
C                = 0  normal return, calculation successful.
C                = K  if the eigenvalue iteration fails to converge,
C                     eigenvalues (and eigenvectors if requested)
C                     1 through K-1 are correct.
C
C      Error Messages
C           No. 1  recoverable  N is greater than LDA
C           No. 2  recoverable  N is less than one.
C           No. 3  recoverable  JOB is nonzero and N is greater than LDV
C           No. 4  warning      LDA > LDV,  elements of A other than the
C                               N by N input elements have been changed
C           No. 5  warning      LDA < LDV,  elements of V other than the
C                               N by N output elements have been changed
C           No. 6  recoverable  nonreal element on diagonal of A.
C
C
C     Subroutines Used
C
C     EISPACK-  HTRIBK, HTRIDI, IMTQL2, TQLRAT
C     BLAS-  SCOPY, SCOPYM
C     SLATEC- XERROR
C***REFERENCES  (NONE)
C***ROUTINES CALLED  HTRIBK,HTRIDI,IMTQL2,SCOPY,SCOPYM,TQLRAT,XERROR
C***END PROLOGUE  CHIEV
      INTEGER I,INFO,J,JOB,K,L,LDA,LDV,M,MDIM,MIN0,N
      REAL A(*),E(*),WORK(*),V(*)
C***FIRST EXECUTABLE STATEMENT  CHIEV
      IF(N .GT. LDA) CALL XERROR( 'CHIEV-N .GT. LDA.',17,1,1)
      IF(N .GT. LDA) RETURN
      IF(N .LT. 1) CALL XERROR( 'CHIEV-N .LT. 1',14,2,1)
      IF(N .LT. 1) RETURN
      IF(N .EQ. 1 .AND. JOB .EQ. 0) GO TO 35
      MDIM = 2 * LDA
      IF(JOB .EQ. 0) GO TO 5
      IF(N .GT. LDV) CALL XERROR( 'CHIEV-JOB .NE. 0, AND N .GT. LDV.',
     1  33,3,1)
      IF(N .GT. LDV) RETURN
      IF(N .EQ. 1) GO TO 35
C
C       REARRANGE A IF NECESSARY WHEN LDA.GT.LDV AND JOB .NE.0
C
      MDIM = MIN0(MDIM,2 * LDV)
      IF(LDA.LT.LDV) CALL XERROR( 'CHIEV-LDA.LT.LDV,  ELEMENTS OF V OTHE
     1R THAN THE N BY N OUTPUT ELEMENTS HAVE BEEN CHANGED.',89,5,0)
      IF(LDA.LE.LDV) GO TO 5
      CALL XERROR(  'CHIEV-LDA.GT.LDV, ELEMENTS OF A OTHER THAN THE N BY
     1 N INPUT ELEMENTS HAVE BEEN CHANGED.',87,4,0)
      L = N - 1
      DO 4 J=1,L
         M = 1+J*2*LDV
         K = 1+J*2*LDA
         CALL SCOPY(2*N,A(K),1,A(M),1)
    4 CONTINUE
    5 CONTINUE
C
C     FILL IN LOWER TRIANGLE OF A, COLUMN BY COLUMN.
C
      DO 6 J = 1,N
       K = (J-1)*(MDIM+2)+1
       IF(A(K+1).NE.0.0) CALL XERROR
     1     ( 'CHIEV-NONREAL ELEMENT ON DIAGONAL OF A',38,6,1)
      IF(A(K+1) .NE.0.0) RETURN
       CALL SCOPY(N-J+1,A(K),MDIM,A(K),2)
       CALL SCOPYM(N-J+1,A(K+1),MDIM,A(K+1),2)
    6 CONTINUE
C
C     SEPARATE REAL AND IMAGINARY PARTS
C
      DO 10 J = 1, N
       K = (J-1) * MDIM +1
       L = K + N
       CALL SCOPY(N,A(K+1),2,WORK(1),1)
       CALL SCOPY(N,A(K),2,A(K),1)
       CALL SCOPY(N,WORK(1),1,A(L),1)
   10 CONTINUE
C
C    REDUCE A TO TRIDIAGONAL MATRIX.
C
      CALL HTRIDI(MDIM,N,A(1),A(N+1),E,WORK(1),WORK(N+1),
     1            WORK(2*N+1))
      IF(JOB .NE. 0) GOTO 15
C
C     EIGENVALUES ONLY.
C
      CALL TQLRAT(N,E,WORK(N+1),INFO)
      RETURN
C
C     EIGENVALUES AND EIGENVECTORS.
C
   15 DO 17 J = 1,N
       K = (J-1) * MDIM + 1
       M = K + N - 1
       DO 16 I = K,M
   16   V(I) = 0.
       I = K + J - 1
       V(I) = 1.
   17 CONTINUE
      CALL IMTQL2(MDIM,N,E,WORK(1),V,INFO)
      IF(INFO .NE. 0) RETURN
      CALL HTRIBK(MDIM,N,A(1),A(N+1),WORK(2*N+1),N,V(1),V(N+1))
C
C    CONVERT EIGENVECTORS TO COMPLEX STORAGE.
C
      DO 20 J = 1,N
       K = (J-1) * MDIM + 1
       I = (J-1) * 2 * LDV + 1
       L = K + N
       CALL SCOPY(N,V(K),1,WORK(1),1)
       CALL SCOPY(N,V(L),1,V(I+1),2)
       CALL SCOPY(N,WORK(1),1,V(I),2)
   20 CONTINUE
      RETURN
C
C     TAKE CARE OF N=1 CASE.
C
   35 IF(A(2) .NE. 0.) CALL XERROR
     1   ( 'CHIEV- NONREAL ELEMENT ON DIAGONAL OF A',39,6,1)
      IF(A(2) .NE. 0.) RETURN
      E(1) = A(1)
      INFO = 0
      IF(JOB .EQ. 0) RETURN
      V(1) = A(1)
      V(2) = 0.
      RETURN
      END
