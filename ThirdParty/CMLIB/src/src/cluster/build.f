      SUBROUTINE BUILD(MM, M, N, A, CLAB, RLAB, TITLE, K, ITER, XMISS,
     1                 DMSUM1, DMSUM2, SUM, IWORK, WORK, CWORK, OUNIT)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      BUILDS CLUSTERS BY THE K-MEANS ALGORITHM, PRINTING THE RESULTS
C      FOR ALL INTERMEDIATE ITERATIONS
C
C   DESCRIPTION
C   -----------
C
C   1.  THE VARIABLES SHOULD BE SCALED SIMILARLY (CLUSTER SUBROUTINE
C       STAND CAN BE USED TO STANDARDIZE THE VARIABLES).
C
C   2.  THE ROUTINE ITERATES FROM 1 TO THE DESIRED NUMBER OF CLUSTERS.
C       THE FIRST ITERATION STARTS WITH THE CLUSTER OF ALL CASES AND
C       COMPUTES THE SUMMARY STATISTICS FOR EACH VARIABLE AND THE
C       DISTANCES FROM EACH CASE TO THE CLUSTER CENTER WITH ALL THE
C       CALCULATIONS BEING PRINTED.  THE SECOND ITERATION DIVIDES THE
C       CLUSTER INTO TWO CLUSTERS, MOVING CASES FROM ONE TO THE OTHER
C       UNTIL EITHER NO FURTHER MOVEMENTS DECREASE THE DISTANCES
C       BETWEEN EACH CASE AND THE CENTER OF ITS ASSIGNED CLUSTER OR THE
C       MAXIMUM NUMBER OF MOVEMENTS PER ITERATION HAS BEEN REACHED.
C       FOR THE THIRD AND SUBSEQUENT ITERATIONS, THE CLUSTER WITH THE
C       LARGEST VARIANCE IS SPLIT AND ITS CASES ARE ASSIGNED TO THE
C       CLUSTER WHOSE MEAN IS THE SMALLEST DISTANCE FROM THE CASE.  THE
C       MEANS ARE THEN UPDATED AND THE PROCESS OF REASSIGNING CASES TO
C       CLUSTERS CONTINUES UNTIL NO REASSIGNMENTS ARE MADE FOR AN
C       ITERATION.
C
C   3.  THE CLUSTERS AND THEIR STATISTICS WILL BE PRINTED OUT AFTER EACH
C       ITERATION ON FORTRAN UNIT OUNIT.
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
C            DIMENSION MUST BE AT LEAST N (UNCHANGED ON OUTPUT).
C         THE MATRIX OF DATA VALUES.
C
C         A(I,J) IS THE VALUE FOR THE J-TH VARIABLE FOR THE I-TH CASE.
C
C   CLAB  VECTOR OF 4-CHARACTER VARIABLES DIMENSIONED AT LEAST N.
C            (UNCHANGED ON OUTPUT).
C         THE LABELS OF THE VARIABLES.
C
C   RLAB  VECTOR OF 4-CHARACTER VARIABLES DIMENSIONED AT LEAST M.
C            (UNCHANGED ON OUTPUT).
C         THE LABELS OF THE CASES.
C
C   TITLE 10-CHARACTER VARIABLE (UNCHANGED ON OUTPUT).
C         TITLE OF THE DATA SET.
C
C   K     INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE NUMBER OF CLUSTERS DESIRED.
C
C   ITER  INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         MAXIMUM NUMBER OF MOVEMENTS ALLOWED PER ITERATION.
C
C   XMISS REAL SCALAR (UNCHANGED ON OUTPUT).
C         MISSING VALUE CODE.  IF A(I,J) = XMISS, THEN THE VALUE FOR THE
C         J-TH VARIABLE FOR THE I-TH CASE IS ASSUMED TO BE MISSING.
C
C   DMSUM1 INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE FIRST DIMENSION OF THE MATRIX SUM.  MUST BE AT LEAST 7.
C
C   DMSUM2 INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE SECOND DIMENSION OF THE MATRIX SUM. MUST BE AT LEAST N.
C
C   SUM   REAL MATRIX WHOSE FIRST DIMENSION MUST BE DMSUM1, WHOSE SECOND
C            DIMENSION MUST BE DMSUM2 AND WHOSE THIRD DIMENSION MUST
C            BE AT LEAST K+1.
C         WORK MATRIX.
C
C   IWORK INTEGER VECTOR DIMENSIONED AT LEAST M.
C         WORK VECTOR.
C
C   WORK  REAL VECTOR DIMENSIONED AT LEAST 2*N+2*M.
C         WORK VECTOR.
C
C   CWORK VECTOR OF 4-CHARACTER VARIABLES DIMENSIONED AT LEAST N.
C         WORK VECTOR.
C
C   OUNIT INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         UNIT NUMBER FOR OUTPUT.
C
C   REFERENCE
C   ---------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGES 84-108.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      INTEGER DMSUM1, DMSUM2, DCLUS, OUNIT
      DIMENSION SUM(DMSUM1,DMSUM2,*), A(MM,*), WORK(*), IWORK(*)
      CHARACTER*4 CLAB(*), RLAB(*), CWORK(*)
      CHARACTER*10 TITLE
C
C     SUM(1,J,I) IS THE VALUE OF THE J-TH VARIABLE AT THE CENTER OF
C                   CLUSTER I
C     SUM(2,J,I) IS THE NUMBER OF NON-MISSING OBSERVATIONS FOR THE J-TH
C                   VARIABLE IN CLUSTER I
C     SUM(3,J,I) IS THE MEAN OF THE J-TH VARIABLE IN CLUSTER I
C     SUM(4,J,I) IS THE STANDARD DEVIATION OF THE J-TH VARIABLE IN
C                   CLUSTER I
C     SUM(5,J,I) IS THE MINIMUM OF THE J-TH VARIABLE IN CLUSTER I
C     SUM(6,J,I) IS THE MAXIMUM OF THE J-TH VARIABLE IN CLUSTER I
C     SUM(7,J,I) IS THE SUM OF SQUARED DEVIATIONS FOR THE J-TH VARIABLE
C                   FROM THE MEAN OF CLUSTER I
C
C     THE K+1-ST ROW OF SUM STORES THE SAME CALCULATIONS AS ABOVE EXCEPT
C        FOR THE ENTIRE DATA SET RATHER THAN FOR AN INDIVIDUAL CLUSTER
C
      DCLUS = 2*N + M
      DO 10 I=1,7
         DO 10 J=1,N
            DO 10 KK=1,K+1
   10          SUM(I,J,KK)=0.
C
C     LOOP ONCE FOR EACH DESIRED CLUSTER
C
      DO 130 KK=1,K
         DO 60 NC=1,ITER
            ERR=0.
            DO 20 KKK=1,KK
               DO 20 J=1,N
                  IF(NC.EQ.1.OR.SUM(1,J,KKK).NE.SUM(3,J,KKK)) ERR=1.
   20       CONTINUE
C
C     IF NO CHANGES HAVE BEEN MADE, OUTPUT THE CLUSTER
C
            IF(ERR.EQ.0.) GO TO 70
            DO 30 KKK=1,KK
               DO 30 J=1,N
                  SUM(2,J,KKK)=0.
   30       SUM(1,J,KKK)=SUM(3,J,KKK)
            DO 50 I=1,M
               DO 40 J=1,N
   40             WORK(J)=A(I,J)
               IWORK(I)=NC
C
C     FIND BEST CLUSTER FOR CASE I
C
               CALL KMEANS(N, WORK, KK, XMISS, DMSUM1, DMSUM2, SUM,
     *                     IWORK(I), WORK(DCLUS+I))
   50       CONTINUE
   60    CONTINUE
   70    IF (OUNIT .GT. 0) CALL KOUT(MM, M, N, A, CLAB, RLAB, TITLE, KK,
     *                     DMSUM1, DMSUM2, SUM, IWORK, WORK(DCLUS+1),
     *                     WORK(N+1), WORK(M+N+1), CWORK, OUNIT)
 
C
C     CREATE A NEW CLUSTER BY SPLITTING VARIABLE WITH LARGEST WITHIN-
C     CLUSTER VARIANCE AT THAT VALUE OF THAT VARIABLE AT THE CENTER
C     OF THE CLUSTER
C
         SM=0.
         DO 80 J=1,N
            DO 80 KKK=1,KK
               IF(SUM(4,J,KKK).GE.SM) THEN
                  SM=SUM(4,J,KKK)
                  KM=KKK
               ENDIF
   80    CONTINUE
         KN=KK+1
         DO 90 JJ=1,N
            SUM(2,JJ,KM)=0.
            SUM(3,JJ,KM)=0.
            SUM(2,JJ,KN)=0.
   90       SUM(3,JJ,KN)=0.
         DO 110 I=1,M
            IF(IWORK(I).EQ.KM) THEN
               DO 100 JJ=1,N
                  IF(A(I,JJ).NE.XMISS) THEN
                     IF(A(I,JJ).GE.SUM(1,JJ,KM)) THEN
                        SUM(2,JJ,KN)=SUM(2,JJ,KN)+1
                        SUM(3,JJ,KN)=SUM(3,JJ,KN)+A(I,JJ)
                     ELSE
                        SUM(2,JJ,KM)=SUM(2,JJ,KM)+1
                        SUM(3,JJ,KM)=SUM(3,JJ,KM)+A(I,JJ)
                     ENDIF
                  ENDIF
  100          CONTINUE
            ENDIF
  110    CONTINUE
         DO 120 JJ=1,N
            IF(SUM(2,JJ,KN).NE.0.)SUM(3,JJ,KN)=SUM(3,JJ,KN)/SUM(2,JJ,KN)
            IF(SUM(2,JJ,KM).NE.0.)SUM(3,JJ,KM)=SUM(3,JJ,KM)/SUM(2,JJ,KM)
  120    CONTINUE
  130 CONTINUE
      RETURN
      END
