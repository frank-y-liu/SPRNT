      SUBROUTINE DOT(MM, M, N, A, CLAB, RLAB, TITLE, KL, DMIWRK, IWORK,
     *               IWORK1, DMWORK, WORK, WORK1, CWORK1, CWORK2, IERR,
     *               OUNIT)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      CREATES A TREE OF CLUSTERS OF CASES FOR CATEGORICAL DATA BY
C      MINIMUM-MUTATION FITS (CHANGES OF VARIABLE VALUES BETWEEN A
C      CLUSTER AND ITS DIRECT ANCESTOR IN THE TREE ARE MINIMIZED.
C
C   DESCRIPTION
C   -----------
C
C   1.  THE MINIMUM MUTATION FIT ALGORITHM SETS UP A BINARY TREE OF
C       CLUSTERS IN POINTER FORM.  EACH CASE IS A SINGLETON CLUSTER AT
C       THE LEAF LEVEL OF THE TREE AND THE OTHER NODES LINK UP THE
C       CASES INTO CLUSTERS.  THE TREE IS BINARY SINCE EACH CLUSTER OF
C       CARDINALITY GREATER THAN 1 IS DIVIDED UP INTO TWO SUBCLUSTERS.
C       THE ANCESTOR OF A CLUSTER IN THE TREE IS THE SMALLEST CLUSTER
C       THAT CONTAINS IT.  VALUES ARE ASSIGNED TO THE ANCESTORS SO THAT
C       A MINIMUM NUMBER OF MUTATIONS (OR CHANGES OF A VARIABLE VALUE
C       BETWEEN A CLUSTER AND ITS ANCESTOR) OCCUR.
C
C   2.  THE DATA SHOULD BE REAL NUMBERS WHICH CORRESPOND TO NO MORE
C       THAN 36 CATAGORIES.  THE INITIAL TREE IS FORMED BY THE JOINING
C       ALGORITHM (SEE SUBROUTINE JOIN).  THE DATA MATRIX IS DESTROYED
C       DURING EXECUTION.
C
C   3.  THE OUTPUT IS WRITTEN ON FORTRAN UNIT OUNIT AND STARTS WITH THE
C       MATRIX OF CATEGORY VALUES.
C
C   4.  NEXT IS THE STANDARD TREE OUTPUT OF THE CLUSTERS FROM JOIN.  IT
C       LISTS THE CASES VERTICALLY AND HAS HORIZONTAL LINES EMANATING
C       FROM EACH CASE.  EACH CLUSTER WILL CORRESPOND TO A VERTICAL
C       LINE BETWEEN TWO HORIZONTAL LINES.  THE CASES BETWEEN AND
C       INCLUDED IN THE HORIZONTAL LINES ARE THE MEMBERS OF THE
C       CLUSTER.  THE DISTANCE FROM THE CASE NAMES TO THE VERTICAL
C       LINES CORRESPOND TO THE DISTANCE AT WHICH THE TWO CLUSTERS WERE
C       JOINED.
C
C   5.  NEXT IS THE DOT DIAGRAM OF THE CLUSTERS.  IT LISTS THE CASES
C       VERTICALLY, AND EACH COLUMN CORRESPONDS TO THE CLUSTER WHOSE
C       NUMBER IS AT THE BOTTOM OF THE DIAGRAM.  IF THE CASE IS A
C       MEMBER OF CLUSTER I, THEN THE ENTRY IN THE I-TH COLUMN OF THAT
C       CASE'S ROW IS THE LETTER O.  OTHERWISE, THE ENTRY IS A DOT.
C
C   6.  THE LAST DIAGRAM IS AN ARRAY WHOSE ROWS ARE THE VARIABLES AND
C       THE COLUMNS ARE THE CLUSTERS NUMBERED BY THE SAME NUMBERS AS IN
C       THE DIAGRAM IN 4.  THE ENTRY FOR CLUSTER I FOR VARIABLE J IS A
C       DOT IF THE VALUE FROM THE MINIMUM MUTATION FIT FOR CLUSTER I IS
C       THE SAME AS THE VALUE FROM THE MINIMUM MUTATION FIT FOR THE
C       DIRECT ANCESTOR OF CLUSTER I.  OTHERWISE, THE ENTRY IS THE
C       MINIMUM MUTATION VALUE FOR CLUSTER I.
C
C   7.  TO RECOVER THE CATEGORY VALUES FOR A CASE, START WITH THE
C       MINIMUM MUTATION VALUES FOR THE LAST CLUSTER AND GO THROUGH THE
C       CLUSTERS IN DECREASING ORDER CHANGING THE VALUES FOR A VARIABLE
C       IF AND ONLY IF THE CASE IS IN THE CLUSTER (IF THE LETTER O IS
C       IN THE COLUMN FOR THE CLUSTER) AND THE MINIMUM MUTATION VALUE
C       FOR THAT VARIABLE IN THAT CLUSTER IS NOT A DOT.  THE VALUE IS
C       THEN CHANGED TO THE MINIMUM MUTATION VALUE.  AFTER ALL CLUSTERS
C       HAVE BEEN PROCESSED, THE VALUES FOR THE VARIABLES WILL BE THE
C       ORIGINAL CATEGORY VALUES FOR THE CASE.
C
C   INPUT PARAMETERS
C   ----------------
C
C   MM    INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE FIRST DIMENSION OF THE MATRIX A.  MUST BE AT LEAST M.
C
C   M     INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE NUMBER OF CASES.
C
C   N     INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE NUMBER OF VARIABLES.
C
C   A     REAL MATRIX WHOSE FIRST DIMENSION MUST BE MM AND WHOSE SECOND
C            DIMENSION MUST BE AT LEAST N (DESTROYED DURING EXECUTION).
C         THE MATRIX OF CATEGORY VALUES.
C
C         A(I,J) IS THE CATEGORY FOR THE J-TH VARIABLE FOR THE I-TH
C                CASE.
C
C   CLAB  VECTOR OF 4-CHARACTER VARIABLES DIMENSIONED AT LEAST M.
C            (UNCHANGED ON OUTPUT).
C         THE LABELS OF THE VARIABLES.
C
C   RLAB  VECTOR OF 4-CHARACTER VARIABLES DIMENSIONED AT LEAST N.
C            (UNCHANGED ON OUTPUT).
C         THE LABELS OF THE CASES.
C
C   TITLE 10-CHARACTER VARIABLE (UNCHANGED ON OUTPUT).
C         TITLE OF THE DATA SET.
C
C   KL    INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         MAXIMUM CATEGORY VALUE.  MUST BE LESS THAN 36.
C
C         THE CATEGORY VALUES MUST BE BETWEEN 1 AND KL.  IF NOT, AN
C         ERROR MESSAGE WILL BE PRINTED AND EXECUTION WILL BE
C         TERMINATED.
C
C   DMIWRK INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE FIRST DIMENSION OF THE MATRIX IWORK. MUST BE AT LEAST 2*M.
C
C   IWORK INTEGER MATRIX WHOSE FIRST DIMENSION MUST BE DMIWRK AND SECOND
C            DIMENSION MUST BE AT LEAST M.
C         WORK MATRIX.
C
C   IWORK1 INTEGER VECTOR DIMENSIONED AT LEAST 2*M.
C         WORK VECTOR.
C
C   DMWORK INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE FIRST DIMENSION OF THE MATRIX WORK.  MUST BE AT LEAST M.
C
C   WORK  REAL MATRIX WHOSE FIRST DIMENSION MUST BE DMWORK AND WHOSE
C            SECOND DIMENSION MUST BE AT LEAST N.
C         WORK MATRIX.
C
C   WORK1 REAL VECTOR DIMENSIONED AT LEAST 2*M.
C         WORK VECTOR.
C
C   CWORK1 VECTOR OF 4-CHARACTER VARIABLES DIMENSIONED AT LEAST 2*M+1.
C         WORK VECTOR.
C
C   CWORK2 VECTOR OF 1-CHARACTER VARIABLES DIMENSIONED AT LEAST 2*M+10.
C         WORK VECTOR.
C
C   OUNIT INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         UNIT NUMBER FOR OUTPUT.
C
C   OUTPUT PARAMETER
C   ----------------
C
C   IERR  INTEGER SCALAR.
C         ERROR FLAG.
C
C         IERR = 0, NO ERRORS WERE DETECTED DURING EXECUTION
C
C         IERR = 1, EITHER THE FIRST AND LAST CASES OR THE CLUSTER
C                   DIAMETER FOR A CLUSTER IS OUT OF BOUNDS.  THE
C                   CLUSTER AND ITS VALUES ARE PRINTED ON UNIT OUNIT.
C                   EXECUTION WILL CONTINUE WITH QUESTIONABLE RESULTS
C                   FOR THAT CLUSTER.  ERROR FLAG SET IN JOINING
C                   ROUTINE.
C
C         IERR = 2, A CLUSTER BOUNDARY IS ILLEGAL.  THE CLUSTER AND THE
C                   BOUNDARY ARE PRINTED ON UNIT OUNIT.  EXECUTION WILL
C                   CONTINUE WITH QUESTIONABLE RESULTS FOR THAT CLUSTER.
C
C         IERR = 3, TWO CLUSTERS OVERLAP.  THE NUMBERS OF THE TWO
C                   CLUSTERS ARE PRINTED ON UNIT OUNIT.  EXECUTION WILL
C                   CONTINUE WITH QUESTIONABLE RESULTS FOR THAT CLUSTER.
C
C         IERR = 4, A VALUE OF THE VARIABLE IS GREATER THAN THE MAXIMUM
C                   VALUE.  THE VALUE AND THE MAXIMUM VALUE IS PRINTED
C                   ON UNIT OUNIT AND EXECUTION IS TERMINATED.  INCREASE
C                   THE VALUE OF KL OR CHANGE THE DATA VALUE.
C
C   REFERENCES
C   ----------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGES 233-250.
C
C     HARTIGAN, J. A. (1972).  "MINIMUM MUTATION FITS TO A GIVEN TREE"
C        BIOMETRICS.  VOL. 29, PAGES 53-65.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      INTEGER DMIWRK, DMWORK, OUNIT
      DIMENSION A(MM,*), IWORK1(*), WORK1(*), IWORK(DMIWRK,*),
     *          WORK(DMWORK,*)
      CHARACTER*4 CLAB(*), RLAB(*), CWORK1(*)
      CHARACTER*1 AA(36), CWORK2(*)
      CHARACTER*10 TITLE
      DATA AA/'.','1','2','3','4','5','6','7','8','9','A','B',
     *        'C','D','E','F','G','H','I','J','K','L','M','N',
     *        'O','P','Q','R','S','T','U','V','W','X','Y','Z'/
C
      IERR = 0
      KT = 2 * M - 1
C
C     ELIMINATE DECIMAL PARTS OF THE CATEGORY VALUES
C
      DO 10 I = 1 , M
         DO 10 J = 1 , N
  10        A(I,J) = REAL(INT(A(I,J)))
      IF (OUNIT .GT. 0) WRITE(OUNIT,1)
  1   FORMAT('1')
      CALL OUT(MM,M,N,A,CLAB,RLAB,TITLE,OUNIT)
C
C     COPY DATA
C
      DO 30 I = 1 , M
         DO 20 J = 1 , N
  20        WORK(I,J) = A(I,J)
  30     CWORK1(I) = RLAB(I)
C
C     GET TREE USING JOIN
C
      CALL JOIN(MM,M,N,A,RLAB,TITLE,DMIWRK,IWORK,WORK1,CWORK1(M+1),IERR,
     *          OUNIT)
C
C     PERMUTE TABLE TO AGREE WITH THE ORDERING PRODUCED BY JOIN
C     RESTORE DATA DESTROYED BY JOIN
C
      DO 60 I = 1 , M
         DO 50 K = 1 , M
            IF (CWORK1(K) .NE. CWORK1(M+I+1)) GOTO 50
            DO 40 J = 1 , N
  40           A(I,J) = WORK(K,J)
            RLAB(I) = CWORK1(K)
            GOTO 60
  50     CONTINUE
  60  CONTINUE
      CALL NBTOMT(M+1,M,DMIWRK,IWORK,KT,IWORK1,IERR,OUNIT)
      IF (OUNIT .GT. 0) WRITE(OUNIT,1)
C
C     CHECK RANGE OF A
C
      DO 70 I=1,M
         DO 70 J=1,N
            IF(A(I,J).LT.1..OR.A(I,J).GT.KL) THEN
            IF (OUNIT .GT. 0) WRITE(OUNIT,2) I,J,KL
    2       FORMAT(I5,' ,',I5,' TH DATA VALUE IS NOT IN THE RANGE (1,',
     *             I2,')')
            IERR = 4
            RETURN
            ENDIF
  70  CONTINUE
C
C     COUNT NUMBER IN EACH CLUSTER
C
      DO 80 I=1,KT
         IWORK(I,1)=0
   80    IF(I.LE.M) IWORK(I,1)=1
      DO 90 I=1,KT-1
         J=IWORK1(I)
         IWORK(J,1)=IWORK(J,1)+IWORK(I,1)
   90 CONTINUE
      IF (OUNIT .GT. 0) WRITE(OUNIT,3)
    3 FORMAT('0TREE')
C
C     PRINT TREE, A LINE AT A TIME.
C
      DO 110 I=1,M
         J=KT
         IWORK(J,1)=IWORK(J,1)-1
         CWORK2(KT)=AA(25)
         DO 100 L=2,KT
            LL=KT-L+1
            IF(IWORK(LL,1).GT.0) IWORK(LL,1)=-IWORK(LL,1)
            CWORK2(LL)=AA(1)
            IF(IWORK1(LL).EQ.J.AND.IWORK(LL,1).NE.0) THEN
               J=LL
               IWORK(LL,1)=-IWORK(LL,1)-1
               CWORK2(J)=AA(25)
            ENDIF
  100    CONTINUE
         IF (OUNIT .GT. 0) WRITE(OUNIT,7) J,RLAB(J),(CWORK2(L),L=1,KT)
  110 CONTINUE
C
C     GO THROUGH VARIABLES, ONE AT A TIME
C
      IKT = KT / 10
      INTKT = INT(KT)
      IF (OUNIT .GT. 0) THEN
         WRITE(OUNIT,4)((I,J=1,10),I=1,IKT)
    4    FORMAT(22X,4I1,1X,17(5I1,1X))
         WRITE(OUNIT,5)((I,I=0,9),J=1,IKT+1)
    5    FORMAT(10X,4I1,1X,19(5I1,1X))
         DO 120 III = 1 , INTKT
  120       CWORK2(III) = ' '
         IKT = IKT * 10 + 10 + (IKT+1) * 2
         DO 130 III = 1 , IKT
  130       CWORK2(III) = '-'
         WRITE(OUNIT,6) (CWORK2(III),III=1,IKT)
    6    FORMAT(10X,120(A1))
      ENDIF
      DO 160 JJ=1,N
C
C     FIND THE MINIMUM MUTATION FIT FOR EACH VARIABLE
C
         DO 140 I=1,M
  140       WORK1(I)=A(I,JJ)
         CALL MMFIT(M,WORK1,KT,IWORK1,KL,DMIWRK,IWORK,IERR,OUNIT)
         IF (IERR .EQ. 4) RETURN
         DO 150 I=1,KT
C
C     REPLACE WITH A ZERO IF THE MINIMUM MUTATION VALUES CORRESPOND FOR
C     A CLUSTER AND ITS DIRECT ANCESTOR
C
            IF(I.NE.KT) THEN
               II=IWORK1(I)
               IF(WORK1(I).EQ.WORK1(II)) WORK1(I)=0.
            ENDIF
            J=INT(WORK1(I)+1)
            IF(J.LT.1) J=1
            IF(J.GT.35) J=36
            CWORK2(I)=AA(J)
  150    CONTINUE
         IF (OUNIT .GT. 0) WRITE(OUNIT,7) JJ,CLAB(JJ),(CWORK2(I),I=1,KT)
    7    FORMAT(I4,A4,3X,20(3A1,1X,2A1))
  160 CONTINUE
      RETURN
      END