*SIDTS
      SUBROUTINE SIDTS
     +   (N,M,W,RHO,LDRHO,ALPHA,TT,LDTT,T,LDT,DTT,LDDTT)
C***BEGIN PROLOGUE  SIDTS
C***REFER TO  SODR,SODRC
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   860529   (YYMMDD)
C***REVISION DATE  870204   (YYMMDD)
C***CATEGORY NO.  G2E,I1B1
C***KEYWORDS  ORTHOGONAL DISTANCE REGRESSION,
C             NONLINEAR LEAST SQUARES,
C             ERRORS IN VARIABLES
C***AUTHOR  BOGGS, PAUL T.
C             OPTIMIZATION GROUP/SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS, GAITHERSBURG, MD 20899
C           BYRD, RICHARD H.
C             DEPARTMENT OF COMPUTER SCIENCE
C             UNIVERSITY OF COLORADO, BOULDER, CO 80309
C           DONALDSON, JANET R.
C             OPTIMIZATION GROUP/SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS, BOULDER, CO 80303-3328
C           SCHNABEL, ROBERT B.
C             DEPARTMENT OF COMPUTER SCIENCE
C             UNIVERSITY OF COLORADO, BOULDER, CO 80309
C             AND
C             OPTIMIZATION GROUP/SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS, BOULDER, CO 80303-3328
C***PURPOSE  SCALE MATRIX TT BY THE INVERSE OF DT, I.E., COMPUTE
C            DTT = T * INV(DT) WHERE DT = (W*D)**2 + ALPHA*TT**2,
C            W AND D ARE DEFINED BY EQ.2 OF THE PROLOGUE OF SODR
C            AND SODRC, AND TT IS THE SCALING MATRIX FOR THE DELTA'S,
C            ALSO DEFINED IN THE PROLOGUE OF SODR AND SODRC.
C***END PROLOGUE  SIDTS
C
C  VARIABLE DECLARATIONS (ALPHABETICALLY)
C
C  N.B.  THE LOCATIONS OF W, RHO AND TT ACCESSED DEPEND ON THE VALUE
C        OF THE FIRST ELEMENT OF EACH ARRAY AND THE LEADING DIMENSION
C        OF THE DOUBLY SUBSCRIPTED ARRAYS.
C
      REAL ALPHA
C        THE LEVENBERG-MARQUARDT PARAMETER.
      REAL DT
C        THE VALUE OF THE FACTOR DT = INV((W*D)**2+ALPHA*TT**2)
      REAL DTT(LDDTT,M)
C        THE ARRAY DTT = T * INV(DT) WHERE DT = (W*D)**2 + ALPHA*TT**2.
      INTEGER I
C        AN INDEXING VARIABLE.
      INTEGER J
C        AN INDEXING VARIABLE.
      INTEGER LDDTT
C        THE LEADING DIMENSION OF ARRAY DTT.
      INTEGER LDRHO
C        THE LEADING DIMENSION OF ARRAY RHO.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER LDT
C        THE LEADING DIMENSION OF ARRAY T.
      INTEGER LDTT
C        THE LEADING DIMENSION OF ARRAY TT.
      INTEGER M
C        THE NUMBER OF COLUMNS OF DATA IN THE INDEPENDENT VARIABLE.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      REAL ONE
C        THE VALUE 1.0E0.
      REAL RHO(LDRHO,M)
C        THE DELTA WEIGHTS.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      REAL T(LDT,M)
C        THE STEP FOR THE ESTIMATED DELTA'S.
      REAL TERM1
C        THE VALUE OF THE TERM (W(I)*RHO(I,J))**2
      REAL TERM2
C        THE VALUE OF THE TERM ALPHA*TT(I,J)**2
      REAL TT(LDTT,M)
C        THE SCALE USED FOR THE DELTA'S.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      REAL W(N)
C        THE OBSERVATIONAL ERROR WEIGHTS.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      REAL ZERO
C        THE VALUE 0.0E0.
C
C
      DATA ZERO,ONE/0.0E0,1.0E0/
C
C
C***FIRST EXECUTABLE STATEMENT  SIDTS
C
C
      IF (N.EQ.0 .OR. M.EQ.0) RETURN
C
      IF (W(1).GE.ZERO) THEN
         IF (RHO(1,1).GT.ZERO) THEN
            IF (LDRHO.GE.N) THEN
               IF (TT(1,1).GT.ZERO) THEN
                  IF (LDTT.GE.N) THEN
                     DO 1120 J=1,M
                        DO 1110 I=1,N
                           IF (W(I).NE.ZERO .OR. ALPHA.NE.ZERO) THEN
                              DTT(I,J) = T(I,J)/
     +                                   ((W(I)*RHO(I,J))**2 +
     +                                   ALPHA*TT(I,J)**2)
                           ELSE
                              DTT(I,J) = ZERO
                           END IF
 1110                   CONTINUE
 1120                CONTINUE
                  ELSE
                     DO 1140 J=1,M
                        TERM2 = ALPHA*TT(1,J)**2
                        DO 1130 I=1,N
                           IF (W(I).NE.ZERO .OR. ALPHA.NE.ZERO) THEN
                              DTT(I,J) = T(I,J)/
     +                                   ((W(I)*RHO(I,J))**2+TERM2)
                           ELSE
                              DTT(I,J) = ZERO
                           END IF
 1130                   CONTINUE
 1140                CONTINUE
                  END IF
               ELSE
                  TERM2 = ALPHA*TT(1,1)**2
                  DO 1160 J=1,M
                     DO 1150 I=1,N
                        IF (W(I).NE.ZERO .OR. ALPHA.NE.ZERO) THEN
                           DTT(I,J) = T(I,J)/((W(I)*RHO(I,J))**2+TERM2)
                        ELSE
                           DTT(I,J) = ZERO
                        END IF
 1150                CONTINUE
 1160             CONTINUE
               END IF
            ELSE
               IF (TT(1,1).GT.ZERO) THEN
                  IF (LDTT.GE.N) THEN
                     DO 1220 J=1,M
                        DO 1210 I=1,N
                           IF (W(I).NE.ZERO .OR. ALPHA.NE.ZERO) THEN
                              DTT(I,J) = T(I,J)/
     +                                   ((W(I)*RHO(1,J))**2 +
     +                                   ALPHA*TT(I,J)**2)
                           ELSE
                              DTT(I,J) = ZERO
                           END IF
 1210                   CONTINUE
 1220                CONTINUE
                  ELSE
                     DO 1240 J=1,M
                        TERM2 = ALPHA*TT(1,J)**2
                        DO 1230 I=1,N
                           IF (W(I).NE.ZERO .OR. ALPHA.NE.ZERO) THEN
                              DTT(I,J) = T(I,J)/
     +                                   ((W(I)*RHO(1,J))**2+TERM2)
                           ELSE
                              DTT(I,J) = ZERO
                           END IF
 1230                   CONTINUE
 1240                CONTINUE
                  END IF
               ELSE
                  TERM2 = ALPHA*TT(1,1)**2
                  DO 1260 J=1,M
                     DO 1250 I=1,N
                        IF (W(I).NE.ZERO .OR. ALPHA.NE.ZERO) THEN
                           DTT(I,J) = T(I,J)/((W(I)*RHO(1,J))**2+TERM2)
                        ELSE
                           DTT(I,J) = ZERO
                        END IF
 1250                CONTINUE
 1260             CONTINUE
               END IF
            END IF
         ELSE
            IF (TT(1,1).GT.ZERO) THEN
               IF (LDTT.GE.N) THEN
                  DO 1320 J=1,M
                     DO 1310 I=1,N
                        IF (W(I).NE.ZERO .OR. ALPHA.NE.ZERO) THEN
                           DTT(I,J) = T(I,J)/
     +                                ((W(I)*RHO(1,1))**2 +
     +                                ALPHA*TT(I,J)**2)
                        ELSE
                           DTT(I,J) = ZERO
                        END IF
 1310                CONTINUE
 1320             CONTINUE
               ELSE
                  DO 1340 J=1,M
                     TERM2 = ALPHA*TT(1,J)**2
                     DO 1330 I=1,N
                        IF (W(I).NE.ZERO .OR. ALPHA.NE.ZERO) THEN
                           DTT(I,J) = T(I,J)/((W(I)*RHO(1,1))**2+TERM2)
                        ELSE
                           DTT(I,J) = ZERO
                        END IF
 1330                CONTINUE
 1340             CONTINUE
               END IF
            ELSE
               TERM2 = ALPHA*TT(1,1)**2
               DO 1360 J=1,M
                  DO 1350 I=1,N
                     IF (W(I).NE.ZERO .OR. ALPHA.NE.ZERO) THEN
                        DTT(I,J) = T(I,J)/((W(I)*RHO(1,1))**2+TERM2)
                     ELSE
                        DTT(I,J) = ZERO
                     END IF
 1350             CONTINUE
 1360          CONTINUE
            END IF
         END IF
      ELSE
         IF (RHO(1,1).GT.ZERO) THEN
            IF (LDRHO.GE.N) THEN
               IF (TT(1,1).GT.ZERO) THEN
                  IF (LDTT.GE.N) THEN
                     DO 2120 J=1,M
                        DO 2110 I=1,N
                           DTT(I,J) = T(I,J)/
     +                                (RHO(I,J)**2 + ALPHA*TT(I,J)**2)
 2110                   CONTINUE
 2120                CONTINUE
                  ELSE
                     DO 2140 J=1,M
                        TERM2 = ALPHA*TT(1,J)**2
                        DO 2130 I=1,N
                           DTT(I,J) = T(I,J)/(RHO(I,J)**2+TERM2)
 2130                   CONTINUE
 2140                CONTINUE
                  END IF
               ELSE
                  TERM2 = ALPHA*TT(1,1)**2
                  DO 2160 J=1,M
                     DO 2150 I=1,N
                        DTT(I,J) = T(I,J)/(RHO(I,J)**2+TERM2)
 2150                CONTINUE
 2160             CONTINUE
               END IF
            ELSE
               IF (TT(1,1).GT.ZERO) THEN
                  IF (LDTT.GE.N) THEN
                     DO 2220 J=1,M
                        TERM1 = RHO(1,J)**2
                        DO 2210 I=1,N
                           DTT(I,J) = T(I,J)/(TERM1+ALPHA*TT(I,J)**2)
 2210                   CONTINUE
 2220                CONTINUE
                  ELSE
                     DO 2240 J=1,M
                        DT = ONE/(RHO(1,J)**2+ALPHA*TT(1,J)**2)
                        DO 2230 I=1,N
                           DTT(I,J) = T(I,J)*DT
 2230                   CONTINUE
 2240                CONTINUE
                  END IF
               ELSE
                  TERM2 = ALPHA*TT(1,1)**2
                  DO 2260 J=1,M
                     TERM1 = RHO(1,J)**2
                     DT = ONE/(TERM1+TERM2)
                     DO 2250 I=1,N
                        DTT(I,J) = T(I,J)*DT
 2250                CONTINUE
 2260             CONTINUE
               END IF
            END IF
         ELSE
            IF (TT(1,1).GT.ZERO) THEN
               IF (LDTT.GE.N) THEN
                  TERM1 = RHO(1,1)**2
                  DO 2320 J=1,M
                     DO 2310 I=1,N
                        DTT(I,J) = T(I,J)/(TERM1 + ALPHA*TT(I,J)**2)
 2310                CONTINUE
 2320             CONTINUE
               ELSE
                  TERM1 = RHO(1,1)**2
                  DO 2340 J=1,M
                     TERM2 = ALPHA*TT(1,J)**2
                     DT = ONE/(TERM1+TERM2)
                     DO 2330 I=1,N
                        DTT(I,J) = T(I,J)*DT
 2330                CONTINUE
 2340             CONTINUE
               END IF
            ELSE
               DT = ONE/(RHO(1,1)**2+ALPHA*TT(1,1)**2)
               DO 2360 J=1,M
                  DO 2350 I=1,N
                     DTT(I,J) = T(I,J)*DT
 2350             CONTINUE
 2360          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
      END
