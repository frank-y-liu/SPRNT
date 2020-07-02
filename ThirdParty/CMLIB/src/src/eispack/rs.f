      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
C***BEGIN PROLOGUE  RS
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4A1
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Computes eigenvalues and, optionally, eigenvectors of
C            real symmetric matrix.
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     to find the eigenvalues and eigenvectors (if desired)
C     of a REAL SYMMETRIC matrix.
C
C     On Input
C
C        NM  must be set to the row dimension of the two-dimensional
C        array parameters as declared in the calling program
C        dimension statement.
C
C        N  is the order of the matrix  A.
C
C        A  contains the real symmetric matrix.
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
C        FV1  and  FV2  are temporary storage arrays.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  TQL2,TQLRAT,TRED1,TRED2
C***END PROLOGUE  RS
C
      INTEGER N,NM,IERR,MATZ
      REAL A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
C
C***FIRST EXECUTABLE STATEMENT  RS
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END
