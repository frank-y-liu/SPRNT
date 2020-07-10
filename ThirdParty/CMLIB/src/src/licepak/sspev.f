      SUBROUTINE SSPEV(A,N,E,V,LDV,WORK,JOB,INFO)
C***BEGIN PROLOGUE  SSPEV
C***DATE WRITTEN   800808   (YYMMDD)
C***REVISION DATE  840425   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C           Changed call to TRED3 to pass size of A
C           instead of 1.
C***CATEGORY NO.  D4A1
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK,PACKED,REAL,SYMMETRIC
C***AUTHOR  KAHANER, K. K., (NBS)
C           MOLER, C. B., (U. OF NEW MEXICO)
C           STEWART, G. W., (U. OF MARYLAND)
C***PURPOSE  To compute the eigenvalues and, optionally, the eigen-
C            vectors of a real SYMMETRIC matrix stored in packed form.
C***DESCRIPTION
C
C     LICEPACK.  This version dated 08/08/80.
C     David Kahaner, Cleve Moler, Pete Stewart
C          N.B.S.       U.N.M.     N.B.S./U.MD.
C
C     Abstract
C      SSPEV computes the eigenvalues and, optionally, the eigenvectors
C      of a real symmetric matrix stored in packed form.
C
C     Call Sequence Parameters-
C       (The values of parameters marked with * (star) will be  changed
C         by SSPEV.)
C
C        A*      REAL(N*(N+1)/2)
C                real symmetric packed input matrix.  Contains upper
C                triangle and diagonal of A, by column (elements
C                11, 12, 22, 13, 23, 33, ...).
C
C        N       INTEGER
C                set by the user to
C                the order of the matrix A.
C
C        E*      REAL(N)
C                on return from SSPEV, E contains the eigenvalues of A.
C                See also INFO below.
C
C        V*      REAL(LDV,N)
C                on return from SSPEV, if the user has set JOB
C                = 0        V is not referenced.
C                = nonzero  the N eigenvectors of A are stored in the
C                first N columns of V.  See also INFO below.
C
C        LDV     INTEGER
C                set by the user to
C                the leading dimension of the array V if JOB is also
C                set nonzero.  In that case, N must be .LE. LDV.
C                If JOB is set to zero, LDV is not referenced.
C
C        WORK*   REAL(2N)
C                temporary storage vector.  Contents changed by SSPEV.
C
C        JOB     INTEGER
C                set by the user to
C                = 0        eigenvalues only to be calculated by SSPEV.
C                           Neither V nor LDV are referenced.
C                = nonzero  eigenvalues and vectors to be calculated.
C                           In this case, A & V must be distinct arrays.
C
C       INFO*   INTEGER
C               on return from SSPEV, the value of INFO is
C               = 0 for normal return.
C               = K if the eigenvalue iteration fails to converge.
C                   Eigenvalues and vectors 1 through K-1 are correct.
C
C
C     Error Messages-
C          No. 1   recoverable  N is greater than LDV and JOB is nonzero
C          No. 2   recoverable  N is less than one
C
C     Subroutines Used
C
C      EISPACK- IMTQL2, TQLRAT, TRBAK3, TRED3
C      SLATEC- XERROR
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IMTQL2,TQLRAT,TRBAK3,TRED3,XERROR
C***END PROLOGUE  SSPEV
      INTEGER I,INFO,J,K,LDV,M,N
      REAL A(*),E(N),V(LDV,N),WORK(*)
C***FIRST EXECUTABLE STATEMENT  SSPEV
       IF(N .GT. LDV) CALL XERROR( 'SSPEV-N .GT. LDV.',17,1,1)
       IF(N .GT. LDV) RETURN
       IF(N .LT. 1) CALL XERROR( 'SSPEV-N .LT. 1',14,2,1)
       IF(N .LT. 1) RETURN
C
C       CHECK N=1 CASE
C
      E(1) = A(1)
      INFO = 0
      IF(N .EQ. 1) RETURN
C
      IF(JOB.NE.0) GO TO 20
C
C     EIGENVALUES ONLY
C
*
*  Modified length of A from 1 to N*(N+1)/2, so code can be
*  run with array subscript checking turned on.
      CALL TRED3(N,N*(N+1)/2,A,E,WORK(1),WORK(N+1))
      CALL TQLRAT(N,E,WORK(N+1),INFO)
      RETURN
C
C     EIGENVALUES AND EIGENVECTORS
C
*  Modified length of A from 1 to N*(N+1)/2, so code can be
*  run with array subscript checking turned on.
   20 CALL TRED3(N,N*(N+1)/2,A,E,WORK(1),WORK(1))
      DO 30 I = 1, N
        DO 25 J = 1, N
   25     V(I,J) = 0.
   30   V(I,I) = 1.
      CALL IMTQL2(LDV,N,E,WORK,V,INFO)
      M = N
      IF(INFO .NE. 0) M = INFO - 1
      CALL TRBAK3(LDV,N,1,A,M,V)
      RETURN
      END
