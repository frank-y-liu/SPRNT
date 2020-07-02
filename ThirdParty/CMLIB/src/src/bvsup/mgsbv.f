      SUBROUTINE MGSBV(M,N,A,IA,NIV,IFLAG,S,P,IP,INHOMO,V,W,WCND)
C***BEGIN PROLOGUE  MGSBV
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***REFER TO  BVSUP
C
C **********************************************************************
C Orthogonalize a set of N real vectors and determine their rank
C
C **********************************************************************
C INPUT
C **********************************************************************
C   M = Dimension of vectors
C   N = No. of vectors
C   A = Array whose first N cols contain the vectors
C   IA = First dimension of array A (col length)
C   NIV = Number of independent vectors needed
C   INHOMO = 1 Corresponds to having a non-zero particular solution
C   V = Particular solution vector (not included in the pivoting)
C   INDPVT = 1 Means pivoting will not be used
C
C **********************************************************************
C OUTPUT
C **********************************************************************
C   NIV = No. of linear independent vectors in input set
C     A = Matrix whose first NIV cols. contain NIV orthogonal vectors
C         which span the vector space determined by the input vectors
C   IFLAG
C          = 0 success
C          = 1 incorrect input
C          = 2 rank of new vectors less than N
C   P = Decomposition matrix.  P is upper triangular and
C             (old vectors) = (new vectors) * P.
C         The old vectors will be reordered due to pivoting
C         The dimension of p must be .GE. N*(N+1)/2.
C             (  N*(2*N+1) when N .NE. NFCC )
C   IP = Pivoting vector. The dimension of IP must be .GE. N.
C             (  2*N when N .NE. NFCC )
C   S = Square of norms of incoming vectors
C   V = Vector which is orthogonal to the vectors of A
C   W = Orthogonalization information for the vector V
C   WCND = Worst case (smallest) norm decrement value of the
C          vectors being orthogonalized  (represents a test
C          for linear dependence of the vectors)
C **********************************************************************
C***ROUTINES CALLED  PRVEC,SDOT
C***COMMON BLOCKS    ML18JR,ML5MCO
C***END PROLOGUE  MGSBV
C
      DIMENSION A(IA,N),V(*),W(*),P(*),IP(*),S(*)
C
C
      COMMON /ML18JR/NXPTS,NIC,RE,AE,NOPG,MXNON,NDISK,NTAPE,NEQ,
     1              INDPVT,INTEG,TOL,NPS,NTP,NEQIVP,NUMORT,NFCC,ICOCO
C
      COMMON /ML5MCO/URO,SRU,EPS,LPAR,SQOVFL,TWOU,FOURU
C
C***FIRST EXECUTABLE STATEMENT  MGSBV
      IF(M .GT. 0  .AND.  N .GT. 0  .AND.  IA .GE. M) GO TO 10
      IFLAG=1
      RETURN
C
   10 JP=0
      IFLAG=0
      NP1=N+1
      Y=0.0
      M2=M/2
C
C     CALCULATE SQUARE OF NORMS OF INCOMING VECTORS AND SEARCH FOR
C     VECTOR WITH LARGEST MAGNITUDE
C
      J=0
      DO 30 I=1,N
      VL=SDOT(M,A(1,I),1,A(1,I),1)
      S(I)=VL
      IF (N .EQ. NFCC) GO TO 25
      J=2*I-1
      P(J)=VL
      IP(J)=J
   25 J=J+1
      P(J)=VL
      IP(J)=J
      IF(VL .LE. Y) GO TO 30
      Y=VL
      IX=I
   30 CONTINUE
      IF (INDPVT .NE. 1) GO TO 33
      IX=1
      Y=P(1)
   33 LIX=IX
      IF (N .NE. NFCC) LIX=2*IX-1
      P(LIX)=P(1)
      S(NP1)=0.
      IF (INHOMO .EQ. 1) S(NP1)=SDOT(M,V,1,V,1)
      WCND=1.
      NIVN=NIV
      NIV=0
C
      IF(Y .EQ. 0.0) GO TO 170
C **********************************************************************
      DO 140 NR=1,N
      IF (NIVN .EQ. NIV) GO TO 150
      NIV=NR
      IF(IX .EQ. NR) GO TO 80
C
C     PIVOTING OF COLUMNS OF P MATRIX
C
      NN=N
      LIX=IX
      LR=NR
      IF (N .EQ. NFCC) GO TO 40
      NN=NFCC
      LIX=2*IX-1
      LR=2*NR-1
   40 IF(NR .EQ. 1) GO TO 60
      KD=LIX-LR
      KJ=LR
      NRM1=LR-1
      DO 50 J=1,NRM1
      PSAVE=P(KJ)
      JK=KJ+KD
      P(KJ)=P(JK)
      P(JK)=PSAVE
   50 KJ=KJ+NN-J
      JY=JK+NMNR
      JZ=JY-KD
      P(JY)=P(JZ)
   60 IZ=IP(LIX)
      IP(LIX)=IP(LR)
      IP(LR)=IZ
      SV=S(IX)
      S(IX)=S(NR)
      S(NR)=SV
      IF (N .EQ. NFCC) GO TO 69
      IF (NR .EQ. 1) GO TO 67
      KJ=LR+1
      DO 65 K=1,NRM1
      PSAVE=P(KJ)
      JK=KJ+KD
      P(KJ)=P(JK)
      P(JK)=PSAVE
   65 KJ=KJ+NFCC-K
   67 IZ=IP(LIX+1)
      IP(LIX+1)=IP(LR+1)
      IP(LR+1)=IZ
C
C     PIVOTING OF COLUMNS OF VECTORS
C
   69 DO 70 L=1,M
      T=A(L,IX)
      A(L,IX)=A(L,NR)
   70 A(L,NR)=T
C
C     CALCULATE P(NR,NR) AS NORM SQUARED OF PIVOTAL VECTOR
C
   80 JP=JP+1
      P(JP)=Y
      RY=1.0/Y
      NMNR=N-NR
      IF (N .EQ. NFCC) GO TO 85
      NMNR=NFCC-(2*NR-1)
      JP=JP+1
      P(JP)=0.
      KP=JP+NMNR
      P(KP)=Y
   85 IF(NR .EQ. N  .OR.  NIVN .EQ. NIV) GO TO 125
C
C    CALCULATE ORTHOGONAL PROJECTION VECTORS AND SEARCH FOR LARGEST NORM
C
      Y=0.0
      IP1=NR+1
      IX=IP1
C     ****************************************
      DO 120 J=IP1,N
      DOT=SDOT(M,A(1,NR),1,A(1,J),1)
      JP=JP+1
      JQ=JP+NMNR
      IF (N .NE. NFCC) JQ=JQ+NMNR-1
      P(JQ)=P(JP)-DOT*DOT*RY
      P(JP)=DOT*RY
      DO 90 I = 1,M
   90 A(I,J)=A(I,J)-P(JP)*A(I,NR)
      IF (N .EQ. NFCC) GO TO 99
      KP=JP+NMNR
      JP=JP+1
      PJP=RY*PRVEC(M,A(1,NR),A(1,J))
      P(JP)=PJP
      P(KP)=-PJP
      KP=KP+1
      P(KP)=RY*DOT
      DO 95 K=1,M2
      L=M2+K
      A(K,J)=A(K,J)-PJP*A(L,NR)
   95 A(L,J)=A(L,J)+PJP*A(K,NR)
      P(JQ)=P(JQ)-PJP*PJP/RY
C
C     TEST FOR CANCELLATION IN RECURRENCE RELATION
C
   99 IF(P(JQ) .GT. S(J)*SRU) GO TO 100
      P(JQ)=SDOT(M,A(1,J),1,A(1,J),1)
  100 IF(P(JQ) .LE. Y) GO TO 120
      Y=P(JQ)
      IX=J
  120 CONTINUE
      IF (N .NE. NFCC) JP=KP
C     ****************************************
      IF(INDPVT .EQ. 1) IX=IP1
C
C     RECOMPUTE NORM SQUARED OF PIVOTAL VECTOR WITH SCALAR PRODUCT
C
      Y=SDOT(M,A(1,IX),1,A(1,IX),1)
      IF(Y  .LE.  EPS*S(IX))  GO TO 170
      WCND=AMIN1(WCND,Y/S(IX))
C
C     COMPUTE ORTHOGONAL PROJECTION OF PARTICULAR SOLUTION
C
  125 IF(INHOMO .NE. 1) GO TO 140
      LR=NR
      IF (N .NE. NFCC) LR=2*NR-1
      W(LR)=SDOT(M,A(1,NR),1,V,1)*RY
      DO 130 I=1,M
  130 V(I)=V(I)-W(LR)*A(I,NR)
      IF (N .EQ. NFCC) GO TO 140
      LR=2*NR
      W(LR)=RY*PRVEC(M,V,A(1,NR))
      DO 135 K=1,M2
      L=M2+K
      V(K)=V(K)+W(LR)*A(L,NR)
  135 V(L)=V(L)-W(LR)*A(K,NR)
  140 CONTINUE
C **********************************************************************
C
C     TEST FOR LINEAR DEPENDENCE OF PARTICULAR SOLUTION
C
  150 IF(INHOMO .NE. 1) RETURN
      VNORM=SDOT(M,V,1,V,1)
      IF (S(NP1) .NE. 0.) WCND=AMIN1(WCND,VNORM/S(NP1))
      IF(VNORM .GE. EPS*S(NP1)) RETURN
  170 IFLAG=2
      WCND=EPS
      RETURN
      END
