      SUBROUTINE DTRSL(T,LDT,N,B,JOB,INFO)
C***BEGIN PROLOGUE  DTRSL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2A3
C***KEYWORDS  DOUBLE PRECISION,LINEAR ALGEBRA,LINPACK,MATRIX,SOLVE,
C             TRIANGULAR
C***AUTHOR  STEWART, G. W., (U. OF MARYLAND)
C***PURPOSE  Solves systems of the form  T*X=B or  TRANS(T)*X=B  where T
C            is a TRIANGULAR matrix of order N.
C***DESCRIPTION
C
C     DTRSL solves systems of the form
C
C                   T * X = B
C     or
C                   TRANS(T) * X = B
C
C     where T is a triangular matrix of order N.  Here TRANS(T)
C     denotes the transpose of the matrix T.
C
C     On Entry
C
C         T         DOUBLE PRECISION(LDT,N)
C                   T contains the matrix of the system.  The zero
C                   elements of the matrix are not referenced, and
C                   the corresponding elements of the array can be
C                   used to store other information.
C
C         LDT       INTEGER
C                   LDT is the leading dimension of the array T.
C
C         N         INTEGER
C                   N is the order of the system.
C
C         B         DOUBLE PRECISION(N).
C                   B contains the right hand side of the system.
C
C         JOB       INTEGER
C                   JOB specifies what kind of system is to be solved.
C                   If JOB is
C
C                        00   solve T*X=B, T lower triangular,
C                        01   solve T*X=B, T upper triangular,
C                        10   solve TRANS(T)*X=B, T lower triangular,
C                        11   solve TRANS(T)*X=B, T upper triangular.
C
C     On Return
C
C         B         B contains the solution, if INFO .EQ. 0.
C                   Otherwise B is unaltered.
C
C         INFO      INTEGER
C                   INFO contains zero if the system is nonsingular.
C                   Otherwise INFO contains the index of
C                   the first zero diagonal element of T.
C
C     LINPACK.  This version dated 08/14/78 .
C     G. W. Stewart, University of Maryland, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DDOT
C     Fortran MOD
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DDOT
C***END PROLOGUE  DTRSL
      INTEGER LDT,N,JOB,INFO
      DOUBLE PRECISION T(LDT,*),B(*)
C
C
      DOUBLE PRECISION DDOT,TEMP
      INTEGER CASE,J,JJ
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 150
C
C        CHECK FOR ZERO DIAGONAL ELEMENTS.
C
C***FIRST EXECUTABLE STATEMENT  DTRSL
         DO 10 INFO = 1, N
C     ......EXIT
            IF (T(INFO,INFO) .EQ. 0.0D0) GO TO 150
   10    CONTINUE
         INFO = 0
C
C        DETERMINE THE TASK AND GO TO IT.
C
         CASE = 1
         IF (MOD(JOB,10) .NE. 0) CASE = 2
         IF (MOD(JOB,100)/10 .NE. 0) CASE = CASE + 2
         GO TO (20,50,80,110), CASE
C
C        SOLVE T*X=B FOR T LOWER TRIANGULAR
C
   20    CONTINUE
            B(1) = B(1)/T(1,1)
            IF (N .LT. 2) GO TO 40
            DO 30 J = 2, N
               TEMP = -B(J-1)
               CALL DAXPY(N-J+1,TEMP,T(J,J-1),1,B(J),1)
               B(J) = B(J)/T(J,J)
   30       CONTINUE
   40       CONTINUE
         GO TO 140
C
C        SOLVE T*X=B FOR T UPPER TRIANGULAR.
C
   50    CONTINUE
            B(N) = B(N)/T(N,N)
            IF (N .LT. 2) GO TO 70
            DO 60 JJ = 2, N
               J = N - JJ + 1
               TEMP = -B(J+1)
               CALL DAXPY(J,TEMP,T(1,J+1),1,B(1),1)
               B(J) = B(J)/T(J,J)
   60       CONTINUE
   70       CONTINUE
         GO TO 140
C
C        SOLVE TRANS(T)*X=B FOR T LOWER TRIANGULAR.
C
   80    CONTINUE
            B(N) = B(N)/T(N,N)
            IF (N .LT. 2) GO TO 100
            DO 90 JJ = 2, N
               J = N - JJ + 1
               B(J) = B(J) - DDOT(JJ-1,T(J+1,J),1,B(J+1),1)
               B(J) = B(J)/T(J,J)
   90       CONTINUE
  100       CONTINUE
         GO TO 140
C
C        SOLVE TRANS(T)*X=B FOR T UPPER TRIANGULAR.
C
  110    CONTINUE
            B(1) = B(1)/T(1,1)
            IF (N .LT. 2) GO TO 130
            DO 120 J = 2, N
               B(J) = B(J) - DDOT(J-1,T(1,J),1,B(1),1)
               B(J) = B(J)/T(J,J)
  120       CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
