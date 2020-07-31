      SUBROUTINE DBSPEV(T,AD,N,K,NDERIV,X,INEV,SVALUE,WORK)
C***BEGIN PROLOGUE  DBSPEV
C***DATE WRITTEN   800901   (YYMMDD)
C***REVISION DATE  840425   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  E3,K6
C***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
C             SPLINE
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  Calculates the value of the spline and its derivatives at X
C            from the B-representation .
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Reference
C         SIAM J. Numerical Analysis, 14, No. 3, June, 1977, pp.441-472.
C
C     Abstract    **** a double precision routine ****
C         DBSPEV is the BSPLEV routine of the reference.
C
C         DBSPEV calculates the value of the spline and its derivatives
C         at X from the B-representation (T,A,N,K) and returns them in
C         SVALUE(I),I=1,NDERIV, T(K) .LE. X .LE. T(N+1).  AD(I) can be
C         the B-spline coefficients A(I), I=1,N) if NDERIV=1.  Otherwise
C         AD must be computed before hand by a call to DBSPDR (T,A,N,K,
C         NDERIV,AD).  If X=T(I),I=K,N), right limiting values are
C         obtained.
C
C         To compute left derivatives or left limiting values at a
C         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
C
C         DBSPEV calls DINTRV, DBSPVN
C
C     Description of Arguments
C
C         Input      T,AD,X, are double precision
C          T       - knot vector of length N+K
C          AD      - vector of length (2*N-NDERIV+1)*NDERIV/2 containing
C                    the difference table from DBSPDR.
C          N       - number of B-spline coefficients
C                    N = sum of knot multiplicities-K
C          K       - order of the B-spline, K .GE. 1
C          NDERIV  - number of derivatives, 1 .LE. NDERIV .LE. K.
C                    NDERIV=1 gives the zero-th derivative =
C                    function value
C          X       - argument, T(K) .LE. X .LE. T(N+1)
C          INEV    - an initialization parameter which must be set
C                    to 1 the first time DBSPEV is called.
C
C         Output     SVALUE,WORK are double precision
C          INEV    - INEV contains information for efficient process-
C                    ing after the initial call and INEV must not
C                    be changed by the user.  Distinct splines require
C                    distinct INEV parameters.
C          SVALUE  - vector of length NDERIV containing the spline
C                    value in SVALUE(1) and the NDERIV-1 derivatives
C                    in the remaining components.
C          WORK    - work vector of length 3*K
C
C     Error Conditions
C         Improper input is a fatal error.
C***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
C                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
C                 JUNE 1977, PP. 441-472.
C***ROUTINES CALLED  DBSPVN,DINTRV,XERROR
C***END PROLOGUE  DBSPEV
C
C
      INTEGER I,ID,INEV,IWORK,JJ,K,KP1,KP1MN,L,LEFT,LL,MFLAG,
     1 N, NDERIV
      DOUBLE PRECISION AD, SVALUE, SUM, T, WORK, X
C     DIMENSION T(N+K)
      DIMENSION T(*), AD(*), SVALUE(NDERIV), WORK(*)
C***FIRST EXECUTABLE STATEMENT  DBSPEV
      IF(K.LT.1) GO TO 100
      IF(N.LT.K) GO TO 105
      IF(NDERIV.LT.1 .OR. NDERIV.GT.K) GO TO 115
      ID = NDERIV
      CALL DINTRV(T, N+1, X, INEV, I, MFLAG)
      IF (X.LT.T(K)) GO TO 110
      IF (MFLAG.EQ.0) GO TO 30
      IF (X.GT.T(I)) GO TO 110
   20 IF (I.EQ.K) GO TO 120
      I = I - 1
      IF (X.EQ.T(I)) GO TO 20
C
C *I* HAS BEEN FOUND IN (K,N) SO THAT T(I) .LE. X .LT. T(I+1)
C     (OR .LE. T(I+1), IF T(I) .LT. T(I+1) = T(N+1) ).
   30 KP1MN = K + 1 - ID
      KP1 = K + 1
      CALL DBSPVN(T, KP1MN, K, 1, X, I, WORK(1),WORK(KP1),IWORK)
      JJ = (N+N-ID+2)*(ID-1)/2
C     ADIF(LEFTPL,ID) = AD(LEFTPL-ID+1 + (2*N-ID+2)*(ID-1)/2)
C     LEFTPL = LEFT + L
   40 LEFT = I - KP1MN
      SUM = 0.0D0
      LL = LEFT + JJ + 2 - ID
      DO 50 L=1,KP1MN
        SUM = SUM + WORK(L)*AD(LL)
        LL = LL + 1
   50 CONTINUE
      SVALUE(ID) = SUM
      ID = ID - 1
      IF (ID.EQ.0) GO TO 60
      JJ = JJ-(N-ID+1)
      KP1MN = KP1MN + 1
      CALL DBSPVN(T, KP1MN, K, 2, X, I, WORK(1), WORK(KP1),IWORK)
      GO TO 40
C
   60 RETURN
C
C
  100 CONTINUE
      CALL XERROR( ' DBSPEV,  K DOES NOT SATISFY K.GE.1',35,2,1)
      RETURN
  105 CONTINUE
      CALL XERROR( ' DBSPEV,  N DOES NOT SATISFY N.GE.K',35,2,1)
      RETURN
  110 CONTINUE
      CALL XERROR( ' DBSPEV,  X IS NOT IN T(K).LE.X.LE.T(N+1)',41,2,1)
      RETURN
  115 CONTINUE
      CALL XERROR( ' DBSPEV,  NDERIV DOES NOT SATISFY 1.LE.NDERIV.LE.K',
     1 50, 2, 1)
      RETURN
  120 CONTINUE
      CALL XERROR( ' DBSPEV, A LEFT LIMITING VALUE CANNOT BE OBTAINED AT
     1 T(K)',57,2,1)
      RETURN
      END
