      SUBROUTINE LLSIA(A,MDA,M,N,B,MDB,NB,RE,AE,KEY,MODE,NP,KRANK,
     1   KSURE,RNORM,W,LW,IWORK,LIW,INFO)
C***BEGIN PROLOGUE  LLSIA
C***DATE WRITTEN   810801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D9,D5
C***KEYWORDS  LINEAR LEAST SQUARES,QR FACTORIZATION
C***AUTHOR  MANTEUFFEL, T. A., (LANL)
C***PURPOSE  Solves LINEAR LEAST SQUARES problems by performing
C            a QR factorization of the matrix A using Householder
C            transformations. Emphasis is put on detecting possible
C            rank deficiency.
C***DESCRIPTION
C
C     LLSIA computes the least squares solution(s) to the problem AX=B
C     where A is an M by N matrix with M.GE.N and B is the M by NB
C     matrix of right hand sides.  User input bounds on the uncertainty
C     in the elements of A are used to detect numerical rank deficiency.
C     The algorithm employs a row and column pivot strategy to
C     minimize the growth of uncertainty and round-off errors.
C
C     LLSIA requires (MDA+6)*N + (MDB+1)*NB + M dimensioned space
C
C   ******************************************************************
C   *                                                                *
C   *         WARNING - All input arrays are changed on exit.        *
C   *                                                                *
C   ******************************************************************
C     SUBROUTINE LLSIA(A,MDA,M,N,B,MDB,NB,RE,AE,KEY,MODE,NP,
C    1   KRANK,KSURE,RNORM,W,LW,IWORK,LIW,INFO)
C
C     Input..
C
C     A(,)          Linear coefficient matrix of AX=B, with MDA the
C      MDA,M,N      actual first dimension of A in the calling program.
C                   M is the row dimension (no. of EQUATIONS of the
C                   problem) and N the col dimension (no. of UNKNOWNS).
C                   Must have MDA.GE.M and M.GE.N.
C
C     B(,)          Right hand side(s), with MDB the actual first
C      MDB,NB       dimension of B in the calling program. NB is the
C                   number of M by 1 right hand sides. Must have
C                   MDB.GE.M. If NB = 0, B is never accessed.
C
C   ******************************************************************
C   *                                                                *
C   *         Note - Use of RE and AE are what make this             *
C   *                code significantly different from               *
C   *                other linear least squares solvers.             *
C   *                However, the inexperienced user is              *
C   *                advised to set RE=0.,AE=0.,KEY=0.               *
C   *                                                                *
C   ******************************************************************
C     RE(),AE(),KEY
C     RE()          RE() is a vector of length N such that RE(I) is
C                   the maximum relative uncertainty in column I of
C                   the matrix A. The values of RE() must be between
C                   0 and 1. A minimum of 10*machine precision will
C                   be enforced.
C
C     AE()          AE() is a vector of length N such that AE(I) is
C                   the maximum absolute uncertainty in column I of
C                   the matrix A. The values of AE() must be greater
C                   than or equal to 0.
C
C     KEY           For ease of use, RE and AE may be input as either
C                   vectors or scalars. If a scalar is input, the algo-
C                   rithm will use that value for each column of A.
C                   The parameter key indicates whether scalars or
C                   vectors are being input.
C                        KEY=0     RE scalar  AE scalar
C                        KEY=1     RE vector  AE scalar
C                        KEY=2     RE scalar  AE vector
C                        KEY=3     RE vector  AE vector
C
C     MODE          The integer mode indicates how the routine
C                   is to react if rank deficiency is detected.
C                   If MODE = 0 return immediately, no solution
C                             1 compute truncated solution
C                             2 compute minimal length solution
C                   The inexperienced user is advised to set MODE=0
C
C     NP            The first NP columns of A will not be interchanged
C                   with other columns even though the pivot strategy
C                   would suggest otherwise.
C                   The inexperienced user is advised to set NP=0.
C
C     WORK()        A real work array dimensioned 5*N. However,if
C                   RE or AE have been specified as vectors, dimension
C                   WORK 4*N. If both RE and AE have been specified
C                   as vectors, dimension WORK 3*N.
C
C     LW            Actual dimension of WORK
C
C     IWORK()       Integer work array dimensioned at least N+M.
C
C     LIW           Actual dimension of IWORK.
C
C     INFO          Is a flag which provides for the efficient
C                   solution of subsequent problems involving the
C                   same A but different B.
C                   If INFO = 0 original call
C                      INFO = 1 subsequent calls
C                   On subsequent calls, the user must supply A, KRANK,
C                   LW, IWORK, LIW, and the first 2*N locations of WORK
C                   as output by the original call to LLSIA. MODE must
C                   be equal to the value of MODE in the original call.
C                   If MODE.LT.2, only the first N locations of WORK
C                   are accessed. AE, RE, KEY, and NP are not accessed.
C
C     Output..
C
C     A(,)          Contains the upper triangular part of the reduced
C                   matrix and the transformation information. It togeth
C                   with the first N elements of WORK (see below)
C                   completely specify the QR factorization of A.
C
C     B(,)          Contains the N by NB solution matrix for X.
C
C     KRANK,KSURE   The numerical rank of A,  based upon the relative
C                   and absolute bounds on uncertainty, is bounded
C                   above by KRANK and below by KSURE. The algorithm
C                   returns a solution based on KRANK. KSURE provides
C                   an indication of the precision of the rank.
C
C     RNORM()       Contains the Euclidean length of the NB residual
C                   vectors  B(I)-AX(I), I=1,NB.
C
C     WORK()        The first N locations of WORK contain values
C                   necessary to reproduce the Householder transformatio
C
C     IWORK()       The first N locations contain the order in
C                   which the columns of A were used. The next
C                   M locations contain the order in which the
C                   rows of A were used.
C
C     INFO          Flag to indicate status of computation on completion
C                  -1   Parameter error(s)
C                   0 - Rank deficient, no solution
C                   1 - Rank deficient, truncated solution
C                   2 - Rank deficient, minimal length solution
C                   3 - Numerical rank 0, zero solution
C                   4 - Rank .LT. NP
C                   5 - Full rank
C***REFERENCES  MANTEUFFEL, T., *AN INTERVAL ANALYSIS APPROACH TO RANK
C                 DETERMINATION IN LINEAR LEAST SQUARES PROBLEMS*,
C                 SANDIA LABORATORIES REPORT SAND80-0655, JUNE, 1980.
C***ROUTINES CALLED  R1MACH,U11LS,U12LS,XERROR
C***END PROLOGUE  LLSIA
      DIMENSION A(MDA,N),B(MDB,NB),RE(N),AE(N),RNORM(NB),W(LW)
      INTEGER IWORK(LIW)
C
C***FIRST EXECUTABLE STATEMENT  LLSIA
      IF(INFO.LT.0 .OR. INFO.GT.1) GO TO 514
      IT=INFO
      INFO=-1
      IF(NB.EQ.0 .AND. IT.EQ.1) GO TO 501
      IF(M.LT.1) GO TO 502
      IF(N.LT.1) GO TO 503
      IF(N.GT.M) GO TO 504
      IF(MDA.LT.M) GO TO 505
      IF(LIW.LT.M+N) GO TO 506
      IF(MODE.LT.0 .OR. MODE.GT.3) GO TO 515
      IF(NB.EQ.0) GO TO 4
      IF(NB.LT.0) GO TO 507
      IF(MDB.LT.M) GO TO 508
      IF(IT.EQ.0) GO TO 4
      GO TO 400
    4 IF(KEY.LT.0.OR.KEY.GT.3) GO TO 509
      IF(KEY.EQ.0 .AND. LW.LT.5*N) GO TO 510
      IF(KEY.EQ.1 .AND. LW.LT.4*N) GO TO 510
      IF(KEY.EQ.2 .AND. LW.LT.4*N) GO TO 510
      IF(KEY.EQ.3 .AND. LW.LT.3*N) GO TO 510
      IF(NP.LT.0 .OR. NP.GT.N) GO TO 516
C
      EPS=10.*R1MACH(4)
      N1=1
      N2=N1+N
      N3=N2+N
      N4=N3+N
      N5=N4+N
      N6=N5+N
C
      IF(KEY.EQ.1) GO TO 100
      IF(KEY.EQ.2) GO TO 200
      IF(KEY.EQ.3) GO TO 300
C
      IF(RE(1).LT.0.0) GO TO 511
      IF(RE(1).GT.1.0) GO TO 512
      IF(RE(1).LT.EPS) RE(1)=EPS
      IF(AE(1).LT.0.0) GO TO 513
      DO 20 I=1,N
      W(N4-1+I)=RE(1)
      W(N5-1+I)=AE(1)
   20 CONTINUE
      CALL U11LS(A,MDA,M,N,W(N4),W(N5),MODE,NP,KRANK,KSURE,
     1            W(N1),W(N2),W(N3),IWORK(N1),IWORK(N2))
      GO TO 400
C
  100 CONTINUE
      IF(AE(1).LT.0.0) GO TO 513
      DO 120 I=1,N
      IF(RE(I).LT.0.0) GO TO 511
      IF(RE(I).GT.1.0) GO TO 512
      IF(RE(I).LT.EPS) RE(I)=EPS
      W(N4-1+I)=AE(1)
  120 CONTINUE
      CALL U11LS(A,MDA,M,N,RE,W(N4),MODE,NP,KRANK,KSURE,
     1            W(N1),W(N2),W(N3),IWORK(N1),IWORK(N2))
      GO TO 400
C
  200 CONTINUE
      IF(RE(1).LT.0.0) GO TO 511
      IF(RE(1).GT.1.0) GO TO 512
      IF(RE(1).LT.EPS) RE(1)=EPS
      DO 220 I=1,N
      W(N4-1+I)=RE(1)
      IF(AE(I).LT.0.0) GO TO 513
  220 CONTINUE
      CALL U11LS(A,MDA,M,N,W(N4),AE,MODE,NP,KRANK,KSURE,
     1            W(N1),W(N2),W(N3),IWORK(N1),IWORK(N2))
      GO TO 400
C
  300 CONTINUE
      DO 320 I=1,N
      IF(RE(I).LT.0.0) GO TO 511
      IF(RE(I).GT.1.0) GO TO 512
      IF(RE(I).LT.EPS) RE(I)=EPS
      IF(AE(I).LT.0.0) GO TO 513
  320 CONTINUE
      CALL U11LS(A,MDA,M,N,RE,AE,MODE,NP,KRANK,KSURE,
     1            W(N1),W(N2),W(N3),IWORK(N1),IWORK(N2))
C
C     DETERMINE INFO
C
  400 IF(KRANK.NE.N) GO TO 402
          INFO=5
          GO TO 410
  402 IF(KRANK.NE.0) GO TO 404
          INFO=3
          GO TO 410
  404 IF(KRANK.GE.NP) GO TO 406
          INFO=4
          RETURN
  406 INFO=MODE
      IF(MODE.EQ.0) RETURN
  410 IF(NB.EQ.0) RETURN
C
C     SOLUTION PHASE
C
      N1=1
      N2=N1+N
      N3=N2+N
      IF(INFO.EQ.2) GO TO 420
      IF(LW.LT.N2-1) GO TO 510
      CALL U12LS(A,MDA,M,N,B,MDB,NB,MODE,KRANK,
     1            RNORM,W(N1),W(N1),IWORK(N1),IWORK(N2))
      RETURN
C
  420 IF(LW.LT.N3-1) GO TO 510
      CALL U12LS(A,MDA,M,N,B,MDB,NB,MODE,KRANK,
     1            RNORM,W(N1),W(N2),IWORK(N1),IWORK(N2))
      RETURN
C
C     ERROR MESSAGES
C
  501 CALL XERROR( ' LLSIA -- SOLUTION ONLY(INFO=1) BUT NO RIGHT HAND SI
     1DE(NB=0)',60,1,0)
      RETURN
  502 CALL XERROR( ' LLSIA -- M.LT.1',16,2,1)
      RETURN
  503 CALL XERROR( ' LLSIA -- N.LT.1',16,2,1)
      RETURN
  504 CALL XERROR( ' LLSIA -- N.GT.M',16,2,1)
      RETURN
  505 CALL XERROR( ' LLSIA -- MDA.LT.M',18,2,1)
      RETURN
  506 CALL XERROR( ' LLSIA -- LIW.LT.M+N',20,2,1)
      RETURN
  507 CALL XERROR( ' LLSIA -- NB.LT.0',17,2,1)
      RETURN
  508 CALL XERROR( ' LLSIA -- MDB.LT.M',18,2,1)
      RETURN
  509 CALL XERROR( ' LLSIA -- KEY OUT OF RANGE',26,2,1)
      RETURN
  510 CALL XERROR( ' LLSIA -- INSUFFICIENT WORK SPACE',33,8,1)
      INFO=-1
      RETURN
  511 CALL XERROR( ' LLSIA -- RE(I) .LT. 0',22,2,1)
      RETURN
  512 CALL XERROR( ' LLSIA -- RE(I) .GT. 1',22,2,1)
      RETURN
  513 CALL XERROR( ' LLSIA -- AE(I) .LT. 0',22,2,1)
      RETURN
  514 CALL XERROR( ' LLSIA -- INFO OUT OF RANGE',27,2,1)
      RETURN
  515 CALL XERROR( ' LLSIA -- MODE OUT OF RANGE',27,2,1)
      RETURN
  516 CALL XERROR( ' LLSIA -- NP OUT OF RANGE',25,2,1)
      RETURN
      END
