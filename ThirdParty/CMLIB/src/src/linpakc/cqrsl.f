      SUBROUTINE CQRSL(X,LDX,N,K,QRAUX,Y,QY,QTY,B,RSD,XB,JOB,INFO)
C***BEGIN PROLOGUE  CQRSL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D9,D2C1
C***KEYWORDS  COMPLEX,LINEAR ALGEBRA,LINPACK,MATRIX,
C             ORTHOGONAL TRIANGULAR,SOLVE
C***AUTHOR  STEWART, G. W., (U. OF MARYLAND)
C***PURPOSE  Applies the output of CQRDC to compute coordinate transfor-
C            mations, projections, and least squares solutions.
C***DESCRIPTION
C
C     CQRSL applies the output of CQRDC to compute coordinate
C     transformations, projections, and least squares solutions.
C     For K .LE. MIN(N,P), let XK be the matrix
C
C            XK = (X(JVPT(1)),X(JVPT(2)), ... ,X(JVPT(K)))
C
C     formed from columnns JVPT(1), ... ,JVPT(K) of the original
C     N x P matrix X that was input to CQRDC (if no pivoting was
C     done, XK consists of the first K columns of X in their
C     original order).  CQRDC produces a factored unitary matrix Q
C     and an upper triangular matrix R such that
C
C              XK = Q * (R)
C                       (0)
C
C     This information is contained in coded form in the arrays
C     X and QRAUX.
C
C     On Entry
C
C        X      COMPLEX(LDX,P).
C               X contains the output of CQRDC.
C
C        LDX    INTEGER.
C               LDX is the leading dimension of the array X.
C
C        N      INTEGER.
C               N is the number of rows of the matrix XK.  It must
C               have the same value as N in CQRDC.
C
C        K      INTEGER.
C               K is the number of columns of the matrix XK.  K
C               must not be greater than (N,P), where P is the
C               same as in the calling sequence to CQRDC.
C
C        QRAUX  COMPLEX(P).
C               QRAUX contains the auxiliary output from CQRDC.
C
C        Y      COMPLEX(N)
C               Y contains an N-vector that is to be manipulated
C               by CQRSL.
C
C        JOB    INTEGER.
C               JOB specifies what is to be computed.  JOB has
C               the decimal expansion ABCDE, with the following
C               meaning.
C
C                    If A .NE. 0, compute QY.
C                    If B,C,D, or E .NE. 0, compute QTY.
C                    If C .NE. 0, compute B.
C                    If D .NE. 0, compute RSD .
C                    If E .NE. 0, compute  XB.
C
C               Note that a request to compute B, RSD, or XB
C               automatically triggers the computation of QTY, for
C               which an array must be provided in the calling
C               sequence.
C
C     On Return
C
C        QY     COMPLEX(N).
C               QY contains Q*Y, if its computation has been
C               requested.
C
C        QTY    COMPLEX(N).
C               QTY contains CTRANS(Q)*Y, if its computation has
C               been requested.  Here CTRANS(Q) is the conjugate
C               transpose of the matrix Q.
C
C        B      COMPLEX(K)
C               B contains the solution of the least squares problem
C
C                    minimize NORM2(Y - XK*B),
C
C               if its computation has been requested.  (Note that
C               if pivoting was requested in CQRDC, the J-th
C               component of B will be associated with column JVPT(J)
C               of the original matrix X that was input into CQRDC.)
C
C        RSD    COMPLEX(N).
C               RSD contains the least squares residual Y - XK*B,
C               if its computation has been requested.  RSD is
C               also the orthogonal projection of Y onto the
C               orthogonal complement of the column space of XK.
C
C        XB     COMPLEX(N).
C               XB contains the least squares approximation XK*B,
C               if its computation has been requested.  XB is also
C               the orthogonal projection of Y onto the column space
C               of X.
C
C        INFO   INTEGER.
C               INFO is zero unless the computation of B has
C               been requested and R is exactly singular.  In
C               this case, INFO is the index of the first zero
C               diagonal element of R and B is left unaltered.
C
C     The parameters QY, QTY, B, RSD, and XB are not referenced
C     if their computation is not requested and in this case
C     can be replaced by dummy variables in the calling program.
C     To save storage, the user may in some cases use the same
C     array for different parameters in the calling sequence.  A
C     frequently occuring example is when one wishes to compute
C     any of B, RSD, or XB and does not need Y or QTY.  In this
C     case one may identify Y, QTY, and one of B, RSD, or XB, while
C     providing separate arrays for anything else that is to be
C     computed.  Thus the calling sequence
C
C          CALL CQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)
C
C     will result in the computation of B and RSD, with RSD
C     overwriting Y.  More generally, each item in the following
C     list contains groups of permissible identifications for
C     a single callinng sequence.
C
C          1. (Y,QTY,B) (RSD) (XB) (QY)
C
C          2. (Y,QTY,RSD) (B) (XB) (QY)
C
C          3. (Y,QTY,XB) (B) (RSD) (QY)
C
C          4. (Y,QY) (QTY,B) (RSD) (XB)
C
C          5. (Y,QY) (QTY,RSD) (B) (XB)
C
C          6. (Y,QY) (QTY,XB) (B) (RSD)
C
C     In any group the value returned in the array allocated to
C     the group corresponds to the last member of the group.
C
C     LINPACK.  This version dated 08/14/78 .
C     Stewart, G. W., University of Maryland, Argonne National Lab.
C
C     CQRSL uses the following functions and subprograms.
C
C     BLAS CAXPY,CCOPY,CDOTC
C     Fortran ABS,AIMAG,MIN0,MOD,REAL
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  CAXPY,CCOPY,CDOTC
C***END PROLOGUE  CQRSL
      INTEGER LDX,N,K,JOB,INFO
      COMPLEX X(LDX,*),QRAUX(*),Y(*),QY(*),QTY(*),B(*),RSD(*),XB(*)
C
      INTEGER I,J,JJ,JU,KP1
      COMPLEX CDOTC,T,TEMP
      LOGICAL CB,CQY,CQTY,CR,CXB
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
C     SET INFO FLAG.
C
C***FIRST EXECUTABLE STATEMENT  CQRSL
      INFO = 0
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
      CQY = JOB/10000 .NE. 0
      CQTY = MOD(JOB,10000) .NE. 0
      CB = MOD(JOB,1000)/100 .NE. 0
      CR = MOD(JOB,100)/10 .NE. 0
      CXB = MOD(JOB,10) .NE. 0
      JU = MIN0(K,N-1)
C
C     SPECIAL ACTION WHEN N=1.
C
      IF (JU .NE. 0) GO TO 40
         IF (CQY) QY(1) = Y(1)
         IF (CQTY) QTY(1) = Y(1)
         IF (CXB) XB(1) = Y(1)
         IF (.NOT.CB) GO TO 30
            IF (CABS1(X(1,1)) .NE. 0.0E0) GO TO 10
               INFO = 1
            GO TO 20
   10       CONTINUE
               B(1) = Y(1)/X(1,1)
   20       CONTINUE
   30    CONTINUE
         IF (CR) RSD(1) = (0.0E0,0.0E0)
      GO TO 250
   40 CONTINUE
C
C        SET UP TO COMPUTE QY OR QTY.
C
         IF (CQY) CALL CCOPY(N,Y,1,QY,1)
         IF (CQTY) CALL CCOPY(N,Y,1,QTY,1)
         IF (.NOT.CQY) GO TO 70
C
C           COMPUTE QY.
C
            DO 60 JJ = 1, JU
               J = JU - JJ + 1
               IF (CABS1(QRAUX(J)) .EQ. 0.0E0) GO TO 50
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -CDOTC(N-J+1,X(J,J),1,QY(J),1)/X(J,J)
                  CALL CAXPY(N-J+1,T,X(J,J),1,QY(J),1)
                  X(J,J) = TEMP
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
         IF (.NOT.CQTY) GO TO 100
C
C           COMPUTE CTRANS(Q)*Y.
C
            DO 90 J = 1, JU
               IF (CABS1(QRAUX(J)) .EQ. 0.0E0) GO TO 80
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -CDOTC(N-J+1,X(J,J),1,QTY(J),1)/X(J,J)
                  CALL CAXPY(N-J+1,T,X(J,J),1,QTY(J),1)
                  X(J,J) = TEMP
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        SET UP TO COMPUTE B, RSD, OR XB.
C
         IF (CB) CALL CCOPY(K,QTY,1,B,1)
         KP1 = K + 1
         IF (CXB) CALL CCOPY(K,QTY,1,XB,1)
         IF (CR .AND. K .LT. N) CALL CCOPY(N-K,QTY(KP1),1,RSD(KP1),1)
         IF (.NOT.CXB .OR. KP1 .GT. N) GO TO 120
            DO 110 I = KP1, N
               XB(I) = (0.0E0,0.0E0)
  110       CONTINUE
  120    CONTINUE
         IF (.NOT.CR) GO TO 140
            DO 130 I = 1, K
               RSD(I) = (0.0E0,0.0E0)
  130       CONTINUE
  140    CONTINUE
         IF (.NOT.CB) GO TO 190
C
C           COMPUTE B.
C
            DO 170 JJ = 1, K
               J = K - JJ + 1
               IF (CABS1(X(J,J)) .NE. 0.0E0) GO TO 150
                  INFO = J
C           ......EXIT
                  GO TO 180
  150          CONTINUE
               B(J) = B(J)/X(J,J)
               IF (J .EQ. 1) GO TO 160
                  T = -B(J)
                  CALL CAXPY(J-1,T,X(1,J),1,B,1)
  160          CONTINUE
  170       CONTINUE
  180       CONTINUE
  190    CONTINUE
         IF (.NOT.CR .AND. .NOT.CXB) GO TO 240
C
C           COMPUTE RSD OR XB AS REQUIRED.
C
            DO 230 JJ = 1, JU
               J = JU - JJ + 1
               IF (CABS1(QRAUX(J)) .EQ. 0.0E0) GO TO 220
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  IF (.NOT.CR) GO TO 200
                     T = -CDOTC(N-J+1,X(J,J),1,RSD(J),1)/X(J,J)
                     CALL CAXPY(N-J+1,T,X(J,J),1,RSD(J),1)
  200             CONTINUE
                  IF (.NOT.CXB) GO TO 210
                     T = -CDOTC(N-J+1,X(J,J),1,XB(J),1)/X(J,J)
                     CALL CAXPY(N-J+1,T,X(J,J),1,XB(J),1)
  210             CONTINUE
                  X(J,J) = TEMP
  220          CONTINUE
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE
      RETURN
      END
