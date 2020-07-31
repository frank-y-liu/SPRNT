      SUBROUTINE LSI(W,MDW,MA,MG,N,PRGOPT,X,RNORM,MODE,WS,IP)
C***BEGIN PROLOGUE  LSI
C***REFER TO  LSEI
C
C     This is a companion subprogram to LSEI( ).
C     The documentation for LSEI( ) has more complete
C     usage instructions.
C     Written by R. J. Hanson, SLA.
C
C     Solve..
C              AX = B,  A  MA by N  (least squares equations)
C     subject to..
C
C              GX.GE.H, G  MG by N  (inequality constraints)
C
C     Input..
C
C      W(*,*) contains  (A B) in rows 1,...,MA+MG, cols 1,...,N+1.
C                       (G H)
C
C     MDW,MA,MG,N
C              contain (resp) var. dimension of W(*,*),
C              and matrix dimensions.
C
C     PRGOPT(*),
C              Program option vector.
C
C     OUTPUT..
C
C      X(*),RNORM
C
C              Solution vector(unless MODE=2), length of AX-B.
C
C      MODE
C              =0   Inequality constraints are compatible.
C              =2   Inequality constraints contradictory.
C
C      WS(*),
C              Working storage of dimension K+N+(MG+2)*(N+7),
C              where K=MAX(MA+MG,N).
C      IP(MG+2*N+1)
C              Integer working storage
C
C     Revised April 16, 1982.
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C           Changed call to HFTI to pass MA as 
C           WS dimension instead of 1.
C   000601  Modified AMAX1 to generic MAX (JEC)
C
C***ROUTINES CALLED  H12,HFTI,LPDP,SASUM,SAXPY,SCOPY,SDOT,SSCAL,SSWAP
C***END PROLOGUE  LSI
C
C     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
C     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
C     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
C     (START EDITING AT LINE WITH C++ IN COLS. 1-3.)
C     /REAL (12 BLANKS)/DOUBLE PRECISION/,/SASUM/DASUM/,/SDOT/DDOT/,
C     / SQRT/ DSQRT/,/SSWAP/DSWAP/,
C     /SCOPY/DCOPY/,/SSCAL/DSCAL/,/SAXPY/DAXPY/,/E0/D0/,/SRELPR/DRELPR/
C
C     SUBROUTINES CALLED
C
C     LPDP          THIS SUBPROGRAM MINIMIZES A SUM OF SQUARES
C                   OF UNKNOWNS SUBJECT TO LINEAR INEQUALITY
C                   CONSTRAINTS.  PART OF THIS PACKAGE.
C
C++
C     SDOT,SSCAL    SUBROUTINES FROM THE BLAS PACKAGE.
C     SAXPY,SASUM,  SEE TRANS. MATH. SOFT., VOL. 5, NO. 3, P. 308.
C     SCOPY,SSWAP
C
C     HFTI          SOLVES AN UNCONSTRAINED LINEAR LEAST SQUARES
C                   PROBLEM.  PART OF THIS PACKAGE.
C
C     H12           SUBROUTINE TO CONSTRUCT AND APPLY A HOUSEHOLDER
C                   TRANSFORMATION.
C
C     SUBROUTINE LSI(W,MDW,MA,MG,N,PRGOPT,X,RNORM,MODE,WS,IP)
C
      REAL             W(MDW,*), PRGOPT(*), RNORM, WS(*), X(*)
      INTEGER IP(*)
      REAL             ANORM, SRELPR, FAC, GAM, HALF, ONE, RB, TAU, TOL,
     1 XNORM, ZERO
      REAL             SASUM, SDOT, SQRT
      LOGICAL COV, SCLCOV
C
      DATA ZERO /0.E0/, SRELPR /0.E0/, ONE /1.E0/, HALF /.5E0/
C
C     COMPUTE MACHINE PRECISION=SRELPR ONLY WHEN NECESSARY.
C***FIRST EXECUTABLE STATEMENT  LSI
      IF (.NOT.(SRELPR.EQ.ZERO)) GO TO 30
c*** changed by RF Boisvert, 19-Feb-92  (fails on HP 9000 Series 300)
      srelpr = r1mach(4)
c      SRELPR = ONE
c   10 IF (ONE+SRELPR.EQ.ONE) GO TO 20
c      SRELPR = SRELPR*HALF
c      GO TO 10
c   20 SRELPR = SRELPR + SRELPR
   30 MODE = 0
      RNORM = ZERO
      M = MA + MG
      NP1 = N + 1
      KRANK = 0
      IF (N.LE.0 .OR. M.LE.0) GO TO 70
      ASSIGN 40 TO IGO994
      GO TO 500
C
C     PROCESS-OPTION-VECTOR
C
C     COMPUTE MATRIX NORM OF LEAST SQUARES EQUAS.
   40 ANORM = ZERO
      DO 50 J=1,N
        ANORM = MAX(ANORM,SASUM(MA,W(1,J),1))
   50 CONTINUE
C
C     SET TOL FOR HFTI( ) RANK TEST.
      TAU = TOL*ANORM
C
C     COMPUTE HOUSEHOLDER ORTHOGONAL DECOMP OF MATRIX.
      IF (N.GT.0) WS(1) = ZERO
      CALL SCOPY(N, WS, 0, WS, 1)
      CALL SCOPY(MA, W(1,NP1), 1, WS, 1)
      K = MAX0(M,N)
      MINMAN = MIN0(MA,N)
      N1 = K + 1
      N2 = N1 + N
      CALL HFTI(W, MDW, MA, N, WS, MA, 1, TAU, KRANK, RNORM, WS(N2),
     1 WS(N1), IP)
      FAC = ONE
      GAM = MA - KRANK
      IF (KRANK.LT.MA .AND. SCLCOV) FAC = RNORM**2/GAM
      ASSIGN 60 TO IGO990
      GO TO 80
C
C     REDUCE-TO-LPDP-AND-SOLVE
   60 CONTINUE
   70 IP(1) = KRANK
      IP(2) = N + MAX0(M,N) + (MG+2)*(N+7)
      RETURN
   80 CONTINUE
C
C     TO REDUCE-TO-LPDP-AND-SOLVE
      MAP1 = MA + 1
C
C     COMPUTE INEQ. RT-HAND SIDE FOR LPDP.
      IF (.NOT.(MA.LT.M)) GO TO 260
      IF (.NOT.(MINMAN.GT.0)) GO TO 160
      DO 90 I=MAP1,M
        W(I,NP1) = W(I,NP1) - SDOT(N,W(I,1),MDW,WS,1)
   90 CONTINUE
      DO 100 I=1,MINMAN
        J = IP(I)
C
C     APPLY PERMUTATIONS TO COLS OF INEQ. CONSTRAINT MATRIX.
        CALL SSWAP(MG, W(MAP1,I), 1, W(MAP1,J), 1)
  100 CONTINUE
C
C     APPLY HOUSEHOLDER TRANSFORMATIONS TO CONSTRAINT MATRIX.
      IF (.NOT.(0.LT.KRANK .AND. KRANK.LT.N)) GO TO 120
      DO 110 II=1,KRANK
        I = KRANK + 1 - II
        L = N1 + I
        CALL H12(2, I, KRANK+1, N, W(I,1), MDW, WS(L-1), W(MAP1,1),
     1   MDW, 1, MG)
  110 CONTINUE
C
C     COMPUTE PERMUTED INEQ. CONSTR. MATRIX TIMES R-INVERSE.
  120 DO 150 I=MAP1,M
        IF (.NOT.(0.LT.KRANK)) GO TO 140
        DO 130 J=1,KRANK
          W(I,J) = (W(I,J)-SDOT(J-1,W(1,J),1,W(I,1),MDW))/W(J,J)
  130   CONTINUE
  140   CONTINUE
  150 CONTINUE
C
C     SOLVE THE REDUCED PROBLEM WITH LPDP ALGORITHM,
C     THE LEAST PROJECTED DISTANCE PROBLEM.
  160 CALL LPDP(W(MAP1,1), MDW, MG, KRANK, N-KRANK, PRGOPT, X, XNORM,
     1 MDLPDP, WS(N2), IP(N+1))
      IF (.NOT.(MDLPDP.EQ.1)) GO TO 240
      IF (.NOT.(KRANK.GT.0)) GO TO 180
C
C     COMPUTE SOLN IN ORIGINAL COORDINATES.
      DO 170 II=1,KRANK
        I = KRANK + 1 - II
        X(I) = (X(I)-SDOT(II-1,W(I,I+1),MDW,X(I+1),1))/W(I,I)
  170 CONTINUE
C
C     APPLY HOUSEHOLDER TRANS. TO SOLN VECTOR.
  180 IF (.NOT.(0.LT.KRANK .AND. KRANK.LT.N)) GO TO 200
      DO 190 I=1,KRANK
        L = N1 + I
        CALL H12(2, I, KRANK+1, N, W(I,1), MDW, WS(L-1), X, 1, 1, 1)
  190 CONTINUE
  200 IF (.NOT.(MINMAN.GT.0)) GO TO 230
C
C     REPERMUTE VARIABLES TO THEIR INPUT ORDER.
      DO 210 II=1,MINMAN
        I = MINMAN + 1 - II
        J = IP(I)
        CALL SSWAP(1, X(I), 1, X(J), 1)
  210 CONTINUE
C
C     VARIABLES ARE NOW IN ORIG. COORDINATES.
C     ADD SOLN OF UNSCONSTRAINED PROB.
      DO 220 I=1,N
        X(I) = X(I) + WS(I)
  220 CONTINUE
C
C     COMPUTE THE RESIDUAL VECTOR NORM.
      RNORM = SQRT(RNORM**2+XNORM**2)
  230 GO TO 250
  240 MODE = 2
  250 GO TO 270
  260 CALL SCOPY(N, WS, 1, X, 1)
  270 IF (.NOT.(COV .AND. KRANK.GT.0)) GO TO 490
C
C     COMPUTE COVARIANCE MATRIX BASED ON THE ORTHOGONAL DECOMP.
C     FROM HFTI( ).
C
      KRM1 = KRANK - 1
      KRP1 = KRANK + 1
C
C     COPY DIAG. TERMS TO WORKING ARRAY.
      CALL SCOPY(KRANK, W, MDW+1, WS(N2), 1)
C
C     RECIPROCATE DIAG. TERMS.
      DO 280 J=1,KRANK
        W(J,J) = ONE/W(J,J)
  280 CONTINUE
      IF (.NOT.(KRANK.GT.1)) GO TO 310
C
C     INVERT THE UPPER TRIANGULAR QR FACTOR ON ITSELF.
      DO 300 I=1,KRM1
        IP1 = I + 1
        DO 290 J=IP1,KRANK
          W(I,J) = -SDOT(J-I,W(I,I),MDW,W(I,J),1)*W(J,J)
  290   CONTINUE
  300 CONTINUE
C
C     COMPUTE THE INVERTED FACTOR TIMES ITS TRANSPOSE.
  310 DO 330 I=1,KRANK
        DO 320 J=I,KRANK
          W(I,J) = SDOT(KRANK+1-J,W(I,J),MDW,W(J,J),MDW)
  320   CONTINUE
  330 CONTINUE
      IF (.NOT.(KRANK.LT.N)) GO TO 450
C
C     ZERO OUT LOWER TRAPEZOIDAL PART.
C     COPY UPPER TRI. TO LOWER TRI. PART.
      DO 340 J=1,KRANK
        CALL SCOPY(J, W(1,J), 1, W(J,1), MDW)
  340 CONTINUE
      DO 350 I=KRP1,N
        W(I,1) = ZERO
        CALL SCOPY(I, W(I,1), 0, W(I,1), MDW)
  350 CONTINUE
C
C     APPLY RIGHT SIDE TRANSFORMATIONS TO LOWER TRI.
      N3 = N2 + KRP1
      DO 430 I=1,KRANK
        L = N1 + I
        K = N2 + I
        RB = WS(L-1)*WS(K-1)
        IF (.NOT.(RB.LT.ZERO)) GO TO 420
C
C     IF RB.GE.ZERO, TRANSFORMATION CAN BE REGARDED AS ZERO.
        RB = ONE/RB
C
C     STORE UNSCALED RANK-ONE HOUSEHOLDER UPDATE IN WORK ARRAY.
        WS(N3) = ZERO
        CALL SCOPY(N, WS(N3), 0, WS(N3), 1)
        L = N1 + I
        K = N3 + I
        WS(K-1) = WS(L-1)
        DO 360 J=KRP1,N
          K = N3 + J
          WS(K-1) = W(I,J)
  360   CONTINUE
        DO 370 J=1,N
          L = N3 + I
          K = N3 + J
          WS(J) = SDOT(J-I,W(J,I),MDW,WS(L-1),1) + SDOT(N-J+1,W(J,J),1,
     1     WS(K-1),1)
          WS(J) = WS(J)*RB
  370   CONTINUE
        L = N3 + I
        GAM = SDOT(N-I+1,WS(L-1),1,WS(I),1)*RB
        GAM = GAM*HALF
        CALL SAXPY(N-I+1, GAM, WS(L-1), 1, WS(I), 1)
        DO 410 J=I,N
          IF (.NOT.(I.GT.1)) GO TO 390
          IM1 = I - 1
          K = N3 + J
          DO 380 L=1,IM1
            W(J,L) = W(J,L) + WS(K-1)*WS(L)
  380     CONTINUE
  390     K = N3 + J
          DO 400 L=I,J
            IL = N3 + L
            W(J,L) = W(J,L) + WS(J)*WS(IL-1) + WS(L)*WS(K-1)
  400     CONTINUE
  410   CONTINUE
  420   CONTINUE
  430 CONTINUE
C
C     COPY LOWER TRI. TO UPPER TRI. TO SYMMETRIZE THE COVARIANCE MATRIX.
      DO 440 I=1,N
        CALL SCOPY(I, W(I,1), MDW, W(1,I), 1)
  440 CONTINUE
C
C     REPERMUTE ROWS AND COLS.
  450 DO 470 II=1,MINMAN
        I = MINMAN + 1 - II
        K = IP(I)
        IF (.NOT.(I.NE.K)) GO TO 460
        CALL SSWAP(1, W(I,I), 1, W(K,K), 1)
        CALL SSWAP(I-1, W(1,I), 1, W(1,K), 1)
        CALL SSWAP(K-I-1, W(I,I+1), MDW, W(I+1,K), 1)
        CALL SSWAP(N-K, W(I,K+1), MDW, W(K,K+1), MDW)
  460   CONTINUE
  470 CONTINUE
C
C     PUT IN NORMALIZED RESIDUAL SUM OF SQUARES SCALE FACTOR
C     AND SYMMETRIZE THE RESULTING COVARIANCE MARIX.
      DO 480 J=1,N
        CALL SSCAL(J, FAC, W(1,J), 1)
        CALL SCOPY(J, W(1,J), 1, W(J,1), MDW)
  480 CONTINUE
  490 GO TO 540
  500 CONTINUE
C
C     TO PROCESS-OPTION-VECTOR
C
C     THE NOMINAL TOLERANCE USED IN THE CODE,
      TOL = SQRT(SRELPR)
      COV = .FALSE.
      SCLCOV = .TRUE.
      LAST = 1
      LINK = PRGOPT(1)
  510 IF (.NOT.(LINK.GT.1)) GO TO 520
      KEY = PRGOPT(LAST+1)
      IF (KEY.EQ.1) COV = PRGOPT(LAST+2).NE.ZERO
      IF (KEY.EQ.10) SCLCOV = PRGOPT(LAST+2).EQ.ZERO
      IF (KEY.EQ.5) TOL = MAX(SRELPR,PRGOPT(LAST+2))
      NEXT = PRGOPT(LINK)
      LAST = LINK
      LINK = NEXT
      GO TO 510
  520 GO TO 530
  530 GO TO IGO994, (40)
  540 GO TO IGO990, (60)
      END