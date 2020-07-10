      SUBROUTINE CDIST(MM, M, N, A, BCLAB, BRLAB, I, NN, TH, JJ, DC)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      FINDS CLOSEST VARIABLE TO A GIVEN VARIABLE
C
C   DESCRIPTION
C   -----------
C
C   1.  THE DISTANCE FROM THE GIVEN VARIABLE TO EVERY VARIABLE WHICH
C       HASN'T BEEN JOINED IS CALCULATED AND THE CLOSEST VARIABLE AND
C       ITS DISTANCE ARE RETURNED.
C
C   INPUT PARAMETERS
C   ----------------
C
C   MM, M, N, A, BCLAB, BRLAB, TH -- SEE SUBROUTINE SPLIT2
C
C   I     INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE GIVEN VARIABLE.
C
C   NN    INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         NUMBER OF VARIABLES WHICH HAVE NOT BEEN JOINED.
C
C   OUTPUT PARAMETERS
C   -----------------
C
C   JJ    INTEGER SCALAR.
C         THE VARIABLE CLOSEST TO THE GIVEN VARIABLE.
C
C   DC    REAL SCALAR.
C         THE DISTANCE BETWEEN THE VARIABLE FOUND AND THE GIVEN
C            VARIABLE.
C
C   REFERENCE
C   ---------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGE 298.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      DIMENSION A(MM,*)
      INTEGER BCLAB(*), BRLAB(*)
C
      DC=R1MACH(2)/M
      TT=TH
      IF(TT.EQ.0.) TT=1.
      JJ=1
      LL=BCLAB(I)
      IF(LL.LT.0) RETURN
      DO 20 J=1,N
         IF(I.NE.J.AND.BCLAB(J).GE.0) THEN
            L=BCLAB(J)
C
C     COMPUTE THRESHOLD DISTANCE
C
            DN=0.
            DD=0.
            DO 10 K=1,M
               IF(BRLAB(K).GE.0) THEN
                  KK=BRLAB(K)
                  DN=DN+1.
                  DIF=AMAX1(A(KK,L),A(KK,LL))-AMIN1(A(K,I),A(K,J))
                  IF(DIF.GT.TH) DIF=TT
                  DD=DD+DIF
                  IF(DD.GE.DC*NN) GO TO 20
               ENDIF
   10       CONTINUE
            IF(DN.NE.0.) DD=DD/DN
            IF(DN.EQ.0.) DD=TH
C
C     FIND SMALLEST DISTANCE
C
            IF(DD.LT.DC) THEN
               DC=DD
               JJ=J
            ENDIF
         ENDIF
   20 CONTINUE
      RETURN
      END
