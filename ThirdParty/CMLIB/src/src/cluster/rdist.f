      SUBROUTINE RDIST(MM, M, N, A, BCLAB, BRLAB, I, NN, TH, II, DR)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      FINDS CLOSEST CASE TO A GIVEN CASE
C
C   DESCRIPTION
C   -----------
C
C   1.  THE DISTANCE FROM THE GIVEN CASE TO EVERY CASE WHICH HASN'T
C       BEEN JOINED IS CALCULATED AND THE CLOSEST CASE AND ITS DISTANCE
C       ARE RETURNED.
C
C   INPUT PARAMETERS
C   ----------------
C
C   MM, M, N, A, BCLAB, BRLAB, TH -- SEE SUBROUTINE SPLIT2
C
C   I     INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE GIVEN CASE.
C
C   NN    INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         NUMBER OF CASES WHICH HAVE NOT BEEN JOINED.
C
C   OUTPUT PARAMETERS
C   -----------------
C
C   II    INTEGER SCALAR.
C         THE CASE CLOSEST TO THE GIVEN OBJECT.
C
C   DR    REAL SCALAR.
C         THE DISTANCE BETWEEN THE CASE FOUND AND THE GIVEN CASE.
C
C   REFERENCE
C   ---------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGE 297.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      DIMENSION A(MM,*)
      INTEGER BCLAB(*), BRLAB(*)
C
      TT=TH
      IF(TT.EQ.0.) TT=1.
      DR=R1MACH(2)/M
      II=1
      LL=BRLAB(I)
      IF(LL.LT.0) RETURN
      DO 20 J=1,M
         IF(J.NE.I.AND.BRLAB(J).GE.0) THEN
C
C     COMPUTE THRESHOLD DISTANCE
C
            L=BRLAB(J)
            DN=0.
            DD=0.
            DO 10 K=1,N
               IF(BCLAB(K).GE.0) THEN
                  KK=BCLAB(K)
                  DIF=AMAX1(A(L,KK),A(LL,KK))-AMIN1(A(I,K),A(J,K))
                  DN=DN+1.
                  IF(DIF.GT.TH) DIF=TT
                  DD=DD+DIF
                  IF(DD.GE.DR*NN) GO TO 20
               ENDIF
   10       CONTINUE
            IF(DN.NE.0.) DD=DD/DN
            IF(DN.EQ.0.) DD=TH
C
C     FIND SMALLEST DISTANCE
C
            IF(DD.LT.DR) THEN
               DR=DD
               II=J
            ENDIF
         ENDIF
   20 CONTINUE
      RETURN
      END
