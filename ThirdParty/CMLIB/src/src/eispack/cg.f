      SUBROUTINE CG(NM,N,AR,AI,WR,WI,MATZ,ZR,ZI,FV1,FV2,FV3,IERR)
C***BEGIN PROLOGUE  CG
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4A4
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Computes the eigenvalues and, optionally, the eigenvectors
C            of a complex general matrix.
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     to find the eigenvalues and eigenvectors (if desired)
C     of a COMPLEX GENERAL matrix.
C
C     On INPUT
C
C        NM  must be set to the row dimension of the two-dimensional
C        array parameters as declared in the calling program
C        dimension statement.
C
C        N  is the order of the matrix  A=(AR,AI).
C
C        AR  and  AI  contain the real and imaginary parts,
C        respectively, of the complex general matrix.
C
C        MATZ  is an integer variable set equal to zero if
C        only eigenvalues are desired.  Otherwise it is set to
C        any non-zero integer for both eigenvalues and eigenvectors.
C
C     On OUTPUT
C
C        WR  and  WI  contain the real and imaginary parts,
C        respectively, of the eigenvalues.
C
C        ZR  and  ZI  contain the real and imaginary parts,
C        respectively, of the eigenvectors if MATZ is not zero.
C
C        IERR  is an integer output variable set equal to an
C        error completion code described in section 2B of the
C        documentation.  The normal completion code is zero.
C
C        FV1, FV2, and  FV3  are temporary storage arrays.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  CBABK2,CBAL,COMQR,COMQR2,CORTH
C***END PROLOGUE  CG
C
      INTEGER N,NM,IS1,IS2,IERR,MATZ
      REAL AR(NM,N),AI(NM,N),WR(N),WI(N),ZR(NM,N),ZI(NM,N)
      REAL FV1(N),FV2(N),FV3(N)
C
C***FIRST EXECUTABLE STATEMENT  CG
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 CALL  CBAL(NM,N,AR,AI,IS1,IS2,FV1)
      CALL  CORTH(NM,N,IS1,IS2,AR,AI,FV2,FV3)
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  COMQR(NM,N,IS1,IS2,AR,AI,WR,WI,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  COMQR2(NM,N,IS1,IS2,FV2,FV3,AR,AI,WR,WI,ZR,ZI,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  CBABK2(NM,N,IS1,IS2,FV1,N,ZR,ZI)
   50 RETURN
      END
