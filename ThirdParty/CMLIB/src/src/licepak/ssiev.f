      SUBROUTINE SSIEV(A,LDA,N,E,WORK,JOB,INFO)
C***BEGIN PROLOGUE  SSIEV
C***DATE WRITTEN   800808   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D4A1
C***KEYWORDS  EIGENVALUE,EIGENVECTOR,MATRIX,REAL,SYMMETRIC
C***AUTHOR  KAHANER, D. K., (NBS)
C           MOLER, C. B., (U. OF NEW MEXICO)
C           STEWART, G. W., (U. OF MARYLAND)
C***PURPOSE  To compute the eigenvalues and, optionally, the eigen-
C            vectors of a real SYMMETRIC matrix.
C***DESCRIPTION
C
C     LICEPACK.  This version dated 08/08/80.
C     David Kahaner, Cleve Moler, Pete Stewart
C          N.B.S.       U.N.M.     N.B.S./U.MD.
C
C     Abstract
C      SSIEV computes the eigenvalues and, optionally, the eigenvectors
C      of a real symmetric matrix.
C
C     Call Sequence Parameters-
C       (The values of parameters marked with * (star) will be  changed
C         by SSIEV.)
C
C       A*      REAL (LDA,N)
C               real symmetric input matrix.
C               Only the diagonal and upper triangle of A must be input,
C               as SSIEV copies the upper triangle to the lower.
C               That is, the user must define A(I,J), I=1,..N, and J=I,.
C               ..,N.
C               On return from SSIEV, if the user has set JOB
C               = 0        the lower triangle of A has been altered.
C               = nonzero  the N eigenvectors of A are stored in its
C               first N columns.  See also INFO below.
C
C       LDA     INTEGER
C               set by the user to
C               the leading dimension of the array A.
C
C       N       INTEGER
C               set by the user to
C               the order of the matrix A and
C               the number of elements in E.
C
C       E*      REAL (N)
C               on return from SSIEV, E contains the N
C               eigenvalues of A.  See also INFO below.
C
C       WORK*   REAL (2*N)
C               temporary storage vector.  Contents changed by SSIEV.
C
C       JOB     INTEGER
C               set by user on input
C               = 0         only calculate eigenvalues of A.
C               = nonzero   calculate eigenvalues and eigenvectors of A.
C
C       INFO*   INTEGER
C               on return from SSIEV, the value of INFO is
C               = 0 for normal return.
C               = K if the eigenvalue iteration fails to converge.
C                   eigenvalues and vectors 1 through K-1 are correct.
C
C
C     Error Messages-
C          No. 1   recoverable  N is greater than LDA
C          No. 2   recoverable  N is less than one
C
C     Subroutines Used
C
C     EISPACK- TRED1, TRED2, TQLRAT, IMTQL2
C     SLATEC-  XERROR
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IMTQL2,TQLRAT,TRED1,TRED2,XERROR
C***END PROLOGUE  SSIEV
      INTEGER INFO,JOB,LDA,N
      REAL A(LDA,N),E(N),WORK(*)
C***FIRST EXECUTABLE STATEMENT  SSIEV
       IF(N .GT. LDA) CALL XERROR( 'SSIEV-N .GT. LDA.',17,1,1)
       IF(N .GT. LDA) RETURN
       IF(N .LT. 1) CALL XERROR( 'SSIEV-N .LT. 1',14,2,1)
       IF(N .LT. 1) RETURN
C
C       CHECK N=1 CASE
C
      E(1) = A(1,1)
      INFO = 0
      IF(N .EQ. 1) RETURN
C
C     COPY UPPER TRIANGLE TO LOWER
C
      DO 10 J=1,N
      DO 10 I=1,J
         A(J,I)=A(I,J)
   10 CONTINUE
C
      IF(JOB.NE.0) GO TO 20
C
C     EIGENVALUES ONLY
C
      CALL TRED1(LDA,N,A,E,WORK(1),WORK(N+1))
      CALL TQLRAT(N,E,WORK(N+1),INFO)
      RETURN
C
C     EIGENVALUES AND EIGENVECTORS
C
   20 CALL TRED2(LDA,N,A,E,WORK,A)
      CALL IMTQL2(LDA,N,E,WORK,A,INFO)
      RETURN
      END
