      SUBROUTINE RSPLIT(MM, M, N, A, RLAB, IR, KA, TH, IORD, DMIWRK,
     *                  IWORK, DMWORK, WORK)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      FINDS OPTIMAL SPLIT OF THE CASES
C
C   DESCRIPTION
C   -----------
C
C   1.  INITIALLY, THE FIRST CLUSTER CONSISTS OF ALL CASES WITHIN THE
C       BLOCK IR AND THE SECOND CLUSTER IS EMPTY.  THE REDUCTION IN THE
C       WITHIN-CLUSTER SUM OF SQUARES FOR MOVING EACH CASE FROM THE
C       FIRST CLUSTER TO THE SECOND IS CALCULATED.  THE CASE THAT
C       REDUCES THE SUM OF SQUARES THE MOST IS MOVED AND THIS CONTINUES
C       UNTIL ALL CASES ARE MOVED WITH EACH REDUCTION STORED.  THEN THE
C       SPLIT THAT HAD THE SMALLEST REDUCTION OF ALL IS RETURNED AS THE
C       OPTIMUM SPLIT.
C
C   INPUT PARAMETERS
C   ----------------
C
C   MM, M, N, A, RLAB, TH, IORD, DMIWRK, DMWORK -- SEE SUBROUTINE SPLIT2
C
C   IR    INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         NUMBER OF BLOCK TO BE SPLIT.
C
C   KA    INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         NUMBER OF BLOCKS.
C
C   IWORK INTEGER MATRIX WHOSE FIRST DIMENSION MUST BE DMIWRK AND SECOND
C            DIMENSION MUST BE AT LEAST KA.
C         THE MATRIX DEFINING THE BOUNDARIES OF THE BLOCKS.
C
C         IWORK(1,I) IS 1 + THE FIRST ROW IN BLOCK I
C         IWORK(2,I) IS 1 + THE LAST ROW IN BLOCK I
C         IWORK(3,I) IS 1 + THE FIRST COLUMN IN BLOCK I
C         IWORK(4,I) IS 1 + THE LAST COLUMN IN BLOCK I
C
C   WORK  REAL MATRIX WHOSE FIRST DIMENSION MUST BE DMWORK AND SECOND
C            DIMENSION MUST BE AT LEAST MAX(M,N).
C
C         WORK(1,I) = FIRST CASE IN CASE CLUSTER I
C         WORK(2,I) = LAST CASE IN CASE CLUSTER I
C         WORK(3,I) = REDUCTION IN SSQ DUE TO SPLITTING
C         WORK(4,I) = LAST CASE IN FIRST CLUSTER OF SPLIT OF I
C         WORK(5,I) = 1 IF CASE IS INCLUDED IN PRESENT VARIABLE SPLIT
C         WORK(6,I) = NUMBER OF VARIABLES IN I-TH ROW OF PRESENT
C                        VARIABLE SPLIT
C         WORK(7,I) = MEAN OF I-TH CASE, FIRST VARIABLE CLUSTER
C         WORK(8,I) = NUMBER OF VARIABLES SECOND CLUSTER
C         WORK(9,I) = MEAN OF I-TH CASE, SECOND CLUSTER
C
C         WORK(10-18,I) ARE SIMILAR WITH VARIABLES AND CASES REVERSED.
C
C   REFERENCE
C   ---------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGE 277.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      INTEGER DMWORK, DMIWRK
      DIMENSION A(MM,*),IWORK(DMIWRK,*),WORK(DMWORK,*)
      CHARACTER*4 RLAB(*), C
C
      XM=99999.
      DO 10 J=1,N
   10    WORK(14,J)=0.
C
C     LOOK FOR BLOCKS WITHIN THRESHOLD
C
      IL=INT(WORK(1,IR))
      IU=INT(WORK(2,IR))
      DO 40 K=1,KA
         IF(IWORK(1,K).EQ.IL+1.AND.IWORK(2,K).EQ.IU+1) THEN
            JL=IWORK(3,K)
            JU=IWORK(4,K)
            IF(JL.LT.0) JL=-JL
C
C     COMPUTE VARIANCES
C
            NC=0
            DO 30 J=JL-1,JU-1
               S1=0.
               S2=0.
               S3=0.
               DO 20 I=IL,IU
                  IF(A(I,J).NE.XM) THEN
                     S1=S1+1
                     S2=S2+A(I,J)
                     S3=S3+A(I,J)**2
                  ENDIF
   20          CONTINUE
               WORK(15,J)=S1
               IF(S1.NE.0) THEN
                  S3=S3/S1-(S2/S1)**2
                  WORK(16,J)=S2/S1
               ENDIF
               IF(S3.GT.TH) THEN
                  WORK(14,J)=1.
                  NC=1
               ENDIF
   30       CONTINUE
            IF(NC.EQ.0) IWORK(1,K)=-IWORK(1,K)
         ENDIF
   40 CONTINUE
C
C     FIND BEST CASE SPLIT
C
      DO 50 J=1,N
         WORK(17,J)=0.
   50    WORK(18,J)=0.
      DM=0.
      WORK(3,IR)=0.
      WORK(4,IR)=IL
      DO 100 I=IL,IU-1
         II=IU-I+IL
         ID=II
         DD=-R1MACH(2)
         DO 70 L=IL,II
            IF((IORD.NE.1.AND.IORD.NE.3).OR.L.EQ.II) THEN
               DL=0.
               DO 60 J=1,N
                  IF(WORK(14,J).NE.0.AND.A(L,J).NE.XM) THEN
                     DL=DL+(A(L,J)-WORK(16,J))**2*(WORK(15,J)+1)/
     *                     WORK(15,J)
                     DL=DL-(A(L,J)-WORK(18,J))**2*WORK(17,J)/
     *                     (WORK(17,J)+1)
                  ENDIF
   60          CONTINUE
               IF(DL.GT.DD) THEN
                  DD=DL
                  ID=L
               ENDIF
            ENDIF
   70    CONTINUE
C
C     INTERCHANGE ID AND II
C
         DO 80 J=1,N
            CC=A(II,J)
            A(II,J)=A(ID,J)
   80       A(ID,J)=CC
         C = RLAB(II)
         RLAB(II) = RLAB(ID)
         RLAB(ID) = C
C
C     UPDATE MEANS
C
         DO 90 J=1,N
            IF(WORK(14,J).NE.0.AND.A(II,J).NE.XM) THEN
               WORK(15,J)=WORK(15,J)-1.
               IF(WORK(15,J).NE.0.)WORK(16,J)=WORK(16,J)+
     *                            (WORK(16,J)-A(II,J))/WORK(15,J)
               WORK(17,J)=WORK(17,J)+1.
               WORK(18,J)=WORK(18,J)-(WORK(18,J)-A(II,J))/WORK(17,J)
            ENDIF
   90    CONTINUE
         DM=DM+DD
         IF(DM.GE.WORK(3,IR)) THEN
            WORK(3,IR)=DM
            WORK(4,IR)=II-1
         ENDIF
  100 CONTINUE
      RETURN
      END
