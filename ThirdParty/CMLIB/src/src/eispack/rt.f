      SUBROUTINE RT(NM,N,A,W,MATZ,Z,FV1,IERR)
C***BEGIN PROLOGUE  RT
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4A5
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Compute eigenvalues and eigenvectors of a special real
C            tridiagonal matrix.
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     to find the eigenvalues and eigenvectors (if desired)
C     of a special REAL TRIDIAGONAL matrix.
C
C     On Input
C
C        NM  must be set to the row dimension of the two-dimensional
C        array parameters as declared in the calling program
C        dimension statement.
C
C        N  is the order of the matrix  A.
C
C        A  contains the special real tridiagonal matrix in its
C        first three columns.  The subdiagonal elements are stored
C        in the last  N-1  positions of the first column, the
C        diagonal elements in the second column, and the superdiagonal
C        elements in the first  N-1  positions of the third column.
C        elements  A(1,1)  and  A(N,3)  are arbitrary.
C
C        MATZ  is an integer variable set equal to zero if
C        only eigenvalues are desired.  Otherwise it is set to
C        any non-zero integer for both eigenvalues and eigenvectors.
C
C     On Output
C
C        W  contains the eigenvalues in ascending order.
C
C        Z  contains the eigenvectors if MATZ is not zero.
C
C        IERR  is an integer output variable set equal to an
C        error completion code described in section 2B of the
C        documentation.  The normal completion code is zero.
C
C        FV1  is a temporary storage array.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  FIGI,FIGI2,IMTQL1,IMTQL2
C***END PROLOGUE  RT
C
      INTEGER N,NM,IERR,MATZ
      REAL A(NM,3),W(N),Z(NM,N),FV1(N)
C
C***FIRST EXECUTABLE STATEMENT  RT
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  FIGI(NM,N,A,W,FV1,FV1,IERR)
      IF (IERR .GT. 0) GO TO 50
      CALL  IMTQL1(N,W,FV1,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  FIGI2(NM,N,A,W,FV1,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  IMTQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END
