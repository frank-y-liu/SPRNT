      SUBROUTINE RSP(NM,N,NV,A,W,MATZ,Z,FV1,FV2,IERR)
C***BEGIN PROLOGUE  RSP
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4A1
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Compute eigenvalues and, optionally, eigenvectors of
C            real symmetric matrix packed into a one dimensional
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     to find the eigenvalues and eigenvectors (if desired)
C     of a REAL SYMMETRIC PACKED matrix.
C
C     On Input
C
C        NM  must be set to the row dimension of the two-dimensional
C        array parameters as declared in the calling program
C        dimension statement.
C
C        N  is the order of the matrix  A.
C
C        NV  is an integer variable set equal to the
C        dimension of the array  A  as specified for
C        A  in the calling program.  NV  must not be
C        less than  N*(N+1)/2.
C
C        A  contains the lower triangle of the real symmetric
C        packed matrix stored row-wise.
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
C        error completion code described in section 2B OF the
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
C***ROUTINES CALLED  TQL2,TQLRAT,TRBAK3,TRED3
C***END PROLOGUE  RSP
C
      INTEGER I,J,N,NM,NV,IERR,MATZ
      REAL A(NV),W(N),Z(NM,N),FV1(N),FV2(N)
C
C***FIRST EXECUTABLE STATEMENT  RSP
      IF (N .LE. NM) GO TO 5
      IERR = 10 * N
      GO TO 50
    5 IF (NV .GE. (N * (N + 1)) / 2) GO TO 10
      IERR = 20 * N
      GO TO 50
C
   10 CALL  TRED3(N,NV,A,W,FV1,FV2)
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 DO 40 I = 1, N
C
         DO 30 J = 1, N
            Z(J,I) = 0.0E0
   30    CONTINUE
C
         Z(I,I) = 1.0E0
   40 CONTINUE
C
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  TRBAK3(NM,N,NV,A,N,Z)
   50 RETURN
      END
