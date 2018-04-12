      SUBROUTINE CSPLIT(MM, M, N, A, CLAB, IR, KA, TH, IORD, DMIWRK,
     *                  IWORK, DMWORK, WORK)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      FINDS OPTIMAL SPLIT OF VARIABLES
C
C   DESCRIPTION
C   -----------
C
C   1.  INITIALLY, THE FIRST CLUSTER CONSISTS OF ALL VARIABLES WITHIN
C       THE BLOCK IR AND THE SECOND CLUSTER IS EMPTY.  THE REDUCTION IN
C       THE WITHIN-CLUSTER SUM OF SQUARES FOR MOVING EACH VARIABLE
C       FROM THE FIRST CLUSTER TO THE SECOND IS CALCULATED.  THE
C       VARIABLE THAT REDUCES THE SUM OF SQUARES THE MOST IS MOVED AND
C       THIS CONTINUES UNTIL ALL VARIABLES ARE MOVED WITH EACH
C       REDUCTION STORED.  THEN THE SPLIT THAT HAD THE SMALLEST
C       REDUCTION OF ALL IS RETURNED AS THE OPTIMUM SPLIT.
C
C   INPUT PARAMETERS
C   ----------------
C
C   MM, M, N, A, CLAB, TH, IORD, DMIWRK, DMWORK -- SEE SUBROUTINE SPLIT2
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
C                     VARIABLE SPLIT
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
C        SONS, INC., NEW YORK.  PAGE 276.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      INTEGER DMWORK, DMIWRK
      DIMENSION A(MM,*), IWORK(DMIWRK,*), WORK(DMWORK,*)
      CHARACTER*4 CLAB(*), C
C
      XM=99999.
      DO 10 I=1,M
   10    WORK(5,I)=0.
C
C     LOOK FOR BLOCKS WITHIN THRESHOLD
C
      JL=INT(WORK(10,IR))
      JU=INT(WORK(11,IR))
      DO 40 K=1,KA
         IF(IWORK(3,K).EQ.JL+1.AND.IWORK(4,K).EQ.JU+1) THEN
            IL=IWORK(1,K)
            IF(IL.LT.0) IL=-IL
            IU=IWORK(2,K)
C
C     COMPUTE VARIANCES
C
            NC=0
            DO 30 I=IL-1,IU-1
               S1=0.
               S2=0.
               S3=0.
               DO 20 J=JL,JU
                  IF(A(I,J).NE.XM) THEN
                     S1=S1+1
                     S2=S2+A(I,J)
                     S3=S3+A(I,J)**2
                  ENDIF
   20          CONTINUE
               WORK(6,I)=S1
               IF(S1.NE.0.) THEN
                  WORK(7,I)=S2/S1
                  S3=S3/S1-(S2/S1)**2
               ENDIF
               IF(S3.GT.TH) THEN
                  WORK(5,I)=1.
                  NC=1
               ENDIF
   30       CONTINUE
            IF(NC.EQ.0) IWORK(3,K)=-IWORK(3,K)
         ENDIF
   40 CONTINUE
C
C     FIND BEST VARIABLE SPLIT
C
      DO 50 I=1,M
         WORK(8,I)=0.
   50    WORK(9,I)=0.
      DM=0.
      WORK(12,IR)=0.
      WORK(13,IR)=JL
      DO 100 J=JL,JU-1
         JJ=JU-J+JL
         JD=JJ
         DD=-R1MACH(2)
         DO 70 L=JL,JJ
            IF(IORD.LT.2.OR.L.EQ.JJ) THEN
               DL=0.
               DO 60 I=1,M
                  IF(WORK(5,I).NE.0.AND.A(I,L).NE.XM) THEN
                    DL=DL+(A(I,L)-WORK(7,I))**2*(WORK(6,I)+1.)/WORK(6,I)
                    DL=DL-(A(I,L)-WORK(9,I))**2*WORK(8,I)/(WORK(8,I)+1.)
                  ENDIF
   60          CONTINUE
               IF(DL.GT.DD) THEN
                  DD=DL
                  JD=L
               ENDIF
            ENDIF
   70    CONTINUE
C
C     INTERCHANGE JD AND JJ
C
         DO 80 I=1,M
            CC=A(I,JJ)
            A(I,JJ)=A(I,JD)
   80       A(I,JD)=CC
         C = CLAB(JJ)
         CLAB(JJ) = CLAB(JD)
         CLAB(JD) = C
C
C     UPDATE MEANS
C
         DO 90 I=1,M
            IF(WORK(5,I).NE.0..AND.A(I,JJ).NE.XM) THEN
               WORK(6,I)=WORK(6,I)-1.
               IF(WORK(6,I).NE.0.)WORK(7,I)=WORK(7,I)+(WORK(7,I)-
     *                           A(I,JJ))/WORK(6,I)
               WORK(8,I)=WORK(8,I)+1.
               WORK(9,I)=WORK(9,I)-(WORK(9,I)-A(I,JJ))/WORK(8,I)
            ENDIF
   90    CONTINUE
         DM=DM+DD
         IF(DM.GE.WORK(12,IR)) THEN
            WORK(12,IR)=DM
            WORK(13,IR)=JJ-1
         ENDIF
  100 CONTINUE
      RETURN
      END
