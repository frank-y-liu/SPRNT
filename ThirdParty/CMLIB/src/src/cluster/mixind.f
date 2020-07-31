      SUBROUTINE MIXIND(MM, M, N, A, CLAB, RLAB, TITLE, K, DMWORK,
     *                  WORK1, WORK2, OUNIT)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      FITS THE MIXTURE MODEL FROM K MULTIVARIATE NORMALS WHERE K IS
C      THE DESIRED NUMBER OF CLUSTERS.  THE VARIABLES ARE ASSUMED TO
C      HAVE VARIANCE CONSTANT OVER DIFFERENT CLUSTERS
C
C   DESCRIPTION
C   -----------
C
C   1.  THE DATA ARE ASSUMED TO BE A RANDOM SAMPLE OF SIZE M FROM A
C       MIXTURE OF K MULTIVARIATE NORMAL DISTRIBUTIONS IN N DIMENSIONS.
C       THE SUBROUTINE PREDICTS THE DISTRIBUTION THAT EACH OBSERVATION
C       WAS SAMPLED FROM AND HENCE GROUPS THE OBSERVATIONS INTO K
C       CLUSTERS.  SEE PAGE 113 OF THE REFERENCE FOR A FURTHER
C       DESCRIPTION OF THE MIXTURE ALGORITHM.
C
C   2.  THE ROUTINE BEGINS WITH THE CLUSTER OF ALL OBJECTS AND THEN
C       DIVIDES INTO TWO, THEN THREE, ..., THEN FINALLY K CLUSTERS.
C       THE RESULTS ARE PRINTED AFTER EACH DIVISION ON FORTRAN UNIT
C       OUNIT.  THE RESULTS CONSIST OF THE WITHIN-CLUSTER VARIANCES FOR
C       EACH VARIABLE, THEN SETS UP A COLUMN FOR EACH CLUSTER.  THE
C       MIXTURE PROBABILITY IS THE PROBABILITY THAT A NEW OBJECT WILL
C       BE GROUPED INTO THAT CLUSTER.  THEN THE MEANS OF THE VARIABLES
C       FOR THE CLUSTER ARE PRINTED, AS WELL AS THE PROBABILITIES THAT
C       EACH CASE BELONGS TO EACH CLUSTER.
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
C         THE NUMBER OF CLUSTERS.
C
C   DMWORK INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE LEADING DIMENSION OF THE MATRIX WORK1.  MUST BE AT LEAST
C            N+M+1.
C
C   WORK1 REAL MATRIX WHOSE FIRST DIMENSION MUST BE DMWORK AND WHOSE
C            SECOND DIMENSION MUST BE AT LEAST K.
C         WORK MATRIX.
C
C   WORK2 REAL VECTOR DIMENSIONED AT LEAST N.
C         WORK VECTOR.
C
C   OUNIT INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         UNIT NUMBER FOR OUTPUT.
C
C   REFERENCE
C   ---------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGES 113-129.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      INTEGER DMWORK, U, P, PMIX, OUNIT
      DIMENSION A(MM,*), WORK1(DMWORK,*), WORK2(*)
      CHARACTER*4 CLAB(*), RLAB(*)
      CHARACTER*10 TITLE
C
      U = 0
      P = U + N
      PMIX = P + M + 1
      XM=99999.
      TH=.0001
      DO 160 KK=1,K
C
C     IF NOT FIRST PASS, FIND FURTHEST CASE FROM PRESENT MEANS
C
         DM=0.
         IM=1
         IF(KK.NE.1) THEN
            DO 30 I=1,M
               DI=R1MACH(2)/N
               DO 20 KL=1,KK-1
                  DD=0.
                  XC=0.
                  DO 10 J=1,N
                     IF(A(I,J).NE.XM) THEN
                        XC=XC+1.
                        DD=DD+(A(I,J)-WORK1(U+J,KL))**2 /WORK2(J)
                        IF(DD.GT.DI*N) GO TO 20
                     ENDIF
   10             CONTINUE
                  IF(XC.EQ.0.) GO TO 30
                  DD=DD/XC
   20             IF(DD.LT.DI) DI=DD
               IF(DI.GE.DM) THEN
                  DM=DI
                  IM=I
               ENDIF
   30       CONTINUE
         ENDIF
C
C     BEGIN A NEW CLUSTER LABELED KK
C
         DO 40 J=1,N
   40       WORK1(U+J,KK)=A(IM,J)
         WORK1(PMIX,KK)=EXP(0.5*N)
         ITER=25
         DO 150 IT=1,ITER
C
C     UPDATE PROBABILITIES OF BELONGING
C
            DO 90 I=1,M
               PP=0.
               DO 60 KL=1,KK
                  DD=0.
                  DO 50 J=1,N
                     IF(A(I,J).NE.XM.AND.KK.NE.1)
     *                  DD=DD+(A(I,J)-WORK1(U+J,KL))**2/(WORK2(J)*2.)
   50             CONTINUE
                  IF(DD.GT.100.) DD=100.
                  WORK1(P+I,KL)=WORK1(PMIX,KL)*EXP(-DD)
   60          PP=PP+WORK1(P+I,KL)
               IF(PP.NE.0.) THEN
                  PN=0.
                  DO 70 KL=1,KK
                     WORK1(P+I,KL)=WORK1(P+I,KL)/PP
                     IF(WORK1(P+I,KL).LT.TH) WORK1(P+I,KL)=0.
                     PN =PN+WORK1(P+I,KL)
   70             CONTINUE
                  DO 80 KL=1,KK
   80                WORK1(P+I,KL)=WORK1(P+I,KL)/PN
               ENDIF
   90       CONTINUE
C
C     UPDATE MIXTURE PROBABILITIES
C
            DO 100 KL=1,KK
               WORK1(PMIX,KL)=0.
               DO 100 I=1,M
  100             WORK1(PMIX,KL)=WORK1(PMIX,KL)+WORK1(P+I,KL)/M
C
C     UPDATE CLUSTER ESTIMATES, EACH ONE A WEIGHTED MEAN
C
            DO 120 KL=1,KK
               DO 120 J=1,N
                  WORK1(U+J,KL)=0.
                  DO 110 I=1,M
  110                WORK1(U+J,KL)=WORK1(U+J,KL)+A(I,J)*WORK1(P+I,KL)
  120       IF(WORK1(PMIX,KL).NE.0.) WORK1(U+J,KL)=WORK1(U+J,KL)/
     *                                            (WORK1(PMIX,KL)*M)
            DO 140 J=1,N
               WORK2(J)=0.
               DO 130 I=1,M
                  DO 130 KL=1,KK
  130                WORK2(J)=WORK2(J)+(A(I,J)-WORK1(U+J,KL))**2*
     *                                  WORK1(P+I,KL)
  140       WORK2(J)=WORK2(J)/M
  150    CONTINUE
C
C     PRINT RESULTS OF ITERATION
C
         IF (OUNIT .GT. 0) CALL MIXOUT(MM,M,N,A,CLAB,RLAB,TITLE,KK,
     *                                 DMWORK,WORK1,WORK2,OUNIT)
  160 CONTINUE
      RETURN
      END