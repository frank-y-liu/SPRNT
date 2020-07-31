 
      SUBROUTINE CSCALE(A,NRDA,NROW,NCOL,COLS,COLSAV,ROWS,ROWSAV,ANORM,
     1   SCALES,ISCALE,IC)
C***BEGIN PROLOGUE  CSCALE
C***REFER TO  SODS
C
C     This routine scales the matrix A by columns when needed
C***ROUTINES CALLED  SDOT
C***END PROLOGUE  CSCALE
      DIMENSION A(NRDA,NCOL),COLS(NCOL),COLSAV(NCOL),SCALES(NCOL),
     1          ROWS(NROW),ROWSAV(NROW)
C
      DATA TEN4,TEN20/1.E+4,1.E+20/
C
C***FIRST EXECUTABLE STATEMENT  CSCALE
      IF (ISCALE .NE. (-1)) GO TO 25
C
      IF (IC .EQ. 0) GO TO 10
      DO 5 K=1,NCOL
    5    COLS(K)=SDOT(NROW,A(1,K),1,A(1,K),1)
C
   10 ASCALE=ANORM/FLOAT(NCOL)
      DO 20 K=1,NCOL
         CS=COLS(K)
         IF ((CS .GT. TEN4*ASCALE) .OR. (TEN4*CS .LT. ASCALE)) GO TO 50
         IF ((CS .LT. 1./TEN20) .OR. (CS .GT. TEN20)) GO TO 50
   20 CONTINUE
C
   25 DO 30 K=1,NCOL
   30    SCALES(K)=1.
      RETURN
C
   50 ALOG2=ALOG(2.)
      ANORM=0.
      DO 100 K=1,NCOL
         CS=COLS(K)
         IF (CS .NE. 0.) GO TO 60
         SCALES(K)=1.
         GO TO 100
   60    P=ALOG(CS)/ALOG2
         IP=-0.5*P
         S=2.**IP
         SCALES(K)=S
         IF (IC .EQ. 1) GO TO 70
         COLS(K)=S*S*COLS(K)
         ANORM=ANORM+COLS(K)
         COLSAV(K)=COLS(K)
   70    DO 80 J=1,NROW
   80       A(J,K)=S*A(J,K)
  100 CONTINUE
C
      IF (IC .EQ. 0) RETURN
C
      DO 200 K=1,NROW
         ROWS(K)=SDOT(NCOL,A(K,1),NRDA,A(K,1),NRDA)
         ROWSAV(K)=ROWS(K)
  200    ANORM=ANORM+ROWS(K)
      RETURN
      END