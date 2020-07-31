      SUBROUTINE SHIFTB ( AI, IPIVOT, NROWI, NCOLI, LAST,
     *                    AI1, NROWI1, NCOLI1 )
C  SHIFTS THE ROWS IN CURRENT BLOCK, AI, NOT USED AS PIVOT ROWS, IF
C  ANY, I.E., ROWS IPIVOT(LAST+1),...,IPIVOT(NROWI), ONTO THE FIRST
C  MMAX = NROW-LAST ROWS OF THE NEXT BLOCK, AI1, WITH COLUMN LAST+J OF
C  AI  GOING TO COLUMN J , J=1,...,JMAX=NCOLI-LAST. THE REMAINING COL-
C  UMNS OF THESE ROWS OF AI1 ARE ZEROED OUT.
C
C                             PICTURE
C
C       ORIGINAL SITUATION AFTER         RESULTS IN A NEW BLOCK I+1
C       LAST = 2 COLUMNS HAVE BEEN       CREATED AND READY TO BE
C       DONE IN FACTRB (ASSUMING NO      FACTORED BY NEXT FACTRB CALL.
C       INTERCHANGES OF ROWS)
C                   1
C              X  X 1X  X  X           X  X  X  X  X
C                   1
C              0  X 1X  X  X           0  X  X  X  X
C  BLOCK I          1                       ---------------
C  NROWI = 4   0  0 1X  X  X           0  0 1X  X  X  0  01
C  NCOLI = 5        1                       1             1
C  LAST = 2    0  0 1X  X  X           0  0 1X  X  X  0  01
C  -------------------------------          1             1   NEW
C                   1X  X  X  X  X          1X  X  X  X  X1  BLOCK
C                   1                       1             1   I+1
C  BLOCK I+1        1X  X  X  X  X          1X  X  X  X  X1
C  NROWI1= 5        1                       1             1
C  NCOLI1= 5        1X  X  X  X  X          1X  X  X  X  X1
C  -------------------------------          1-------------1
C                   1
C
      INTEGER IPIVOT(NROWI),LAST, IP,J,JMAX,JMAXP1,M,MMAX
      REAL AI(NROWI,NCOLI),AI1(NROWI1,NCOLI1)
      MMAX = NROWI - LAST
      JMAX = NCOLI - LAST
      IF (MMAX .LT. 1 .OR. JMAX .LT. 1) RETURN
C              PUT THE REMAINDER OF BLOCK I INTO AI1
      DO 10 M=1,MMAX
         IP = IPIVOT(LAST+M)
         DO 10 J=1,JMAX
   10       AI1(M,J) = AI(IP,LAST+J)
      IF (JMAX .EQ. NCOLI1)             RETURN
C              ZERO OUT THE UPPER RIGHT CORNER OF AI1
      JMAXP1 = JMAX + 1
      DO 20 J=JMAXP1,NCOLI1
         DO 20 M=1,MMAX
   20       AI1(M,J) = 0.
                                        RETURN
      END
