      SUBROUTINE CSORT(X,Y,N,KFLAG,WORK)
C***BEGIN PROLOGUE   CSORT
C***REVISION  MARCH 26,1985
C***CATEGORY NO.  M1A3
C***KEYWORD(S) STRING SORTING,SORT,SINGLETON QUICKSORT,QUICKSORT
C***DATE WRITTEN  NOVEMBER,1976
C***AUTHOR  JONES R.E., SULLIVAN F., WISNIEWSKI J.A.
C***PURPOSE
C         CSORT SORTS A STRING ARRAY X AND OPTIONALLY MAKES THE SAME
C         INTERCHANGES IN ARRAY Y.  THE ARRAY X MAY BE SORTED IN
C         FORWARD OR REVERSE LEXICOGRAPHIC ORDER.
C         A SLIGHTLY MODIFIED QUICKSORT ALGORITHM IS USED.
C***DESCRIPTION
C
C     ORIGINALLY WRITTEN BY (FOR REAL DATA) BY R.E. JONES
C     MODIFIED BY JOHN A. WISNIEWSKI TO USE THE SINGLETON QUICKSORT
C       ALGORITHM. DATE 18 NOVEMBER 1976.
C     MODIFIED BY FRANCIS SULLIVAN FOR STRING DATA, 1 OCTOBER 1981.
C     Documentation slightly modified by D. Kahaner, 26Mar,1985
C
C     ABSTRACT
C         CSORT SORTS A STRING ARRAY X AND OPTIONALLY MAKES THE SAME
C         INTERCHANGES IN ARRAY Y.  THE ARRAY X MAY BE SORTED IN
C         FORWARD OR REVERSE LEXICOGRAPHIC ORDER.
C         A SLIGHTLY MODIFIED QUICKSORT ALGORITHM IS USED.
C
C     REFERENCE
C         SINGLETON,R.C., ALGORITHM 347, AN EFFICIENT ALGORITHM FOR
C         SORTING WITH MINIMAL STORAGE, CACM,12(3),1969,185-7.
C
C     DESCRIPTION OF PARAMETERS
C         X - ARRAY OF TYPE CHARACTER TO BE SORTED
C         Y - ARRAY OF TYPE CHARACTER TO BE (OPTIONALLY) CARRIED
C             ELEMENTS OF THE ARRAYS X AND Y CAN BE DIFFERENT LENGTHS.
C             THIS WILL BE DETERMINED BY THE DECLARATIONS IN THE USERS
C             CALLING PROGRAM.  FOR  EXAMPLE, TO SORT A 100 ELEMENT
C             ARRAY OF NAMES, EACH OF LENGTH 6, DECLARE X AS
C             CHARACTER X(100)*6
C             MOST COMPUTERS HAVE AN UPPER LIMIT ON THE LENGTH OF
C             CHARACTER VARIABLES, FOR  EXAMPLE THIS IS
C                 SPERRY 1100 FTN FORTRAN -- 255
C                 CYBER   205 FORTRAN 200 -- 65535
C                 CYBER 170/180 FORTRAN 5 -- 32767
C         N - NUMBER OF VALUES IN ARRAY X TO BE SORTED
C         KFLAG - CONTROL PARAMETER
C             =2   SORT X IN LEXICOGRAPHIC ORDER AND CARRY Y ALONG.
C             =1   SORT X IN LEXICOGRAPHIC ORDER (IGNORING Y)
C             =-1  SORT X IN REVERSE LEXICOGRAPHIC ORDER (IGNORING Y)
C             =-2  SORT X IN REVERSE LEXICOGRAPHIC ORDER AND CARRY Y ALO
C         WORK - ARRAY OF LENGTH 4 OF TYPE CHARACTER USED INTERNALLY
C                FOR SCRATCH STORAGE.  NOTE, EACH ELEMENT OF WORK MUST
C                HAVE A LENGTH SPECIFICATION AT LEAST AS GREAT AS THAT
C                OF THE LARGER OF X AND Y (SEE EXAMPLE BELOW).
C
C     E X A M P L E  O F  U S E
C
C      CHARACTER*2 X,Y,W
C      DIMENSION X(10),W(4)
C      DATA (X(I),I=1,10)/ '05','I ',' I','  ','Rs','9R','R9','89',',*'
C     1 ,'N"'/
C      CALL CSORT(X,Y,10,1,W)
C      PRINT 100,(X(I),I=1,10)
C  100 FORMAT(2X,A2)
C      STOP
C      END
C
C***REFERENCE(S)
C         SINGLETON,R.C., ALGORITHM 347, AN EFFICIENT ALGORITHM FOR
C         SORTING WITH MINIMAL STORAGE, CACM,12(3),1969,185-7.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE
C
      DIMENSION X(N),Y(N),IL(21),IU(21),WORK(4)
      CHARACTER X*(*),Y*(*),WORK*(*)
C***FIRST EXECUTABLE STATEMENT    CSORT
      NN = N
      IF (NN.GE.1) GO TO 10
      CALL XERROR ('CSORT- THE NUMBER OF VALUES TO BE SORTED WAS NOT POS
     1ITIVE.',58,1,1)
      RETURN
   10 KK = IABS(KFLAG)
      IF ((KK.EQ.1).OR.(KK.EQ.2)) GO TO 15
      CALL XERROR ('CSORT- THE SORT CONTROL PARAMETER, K, WAS NOT 2, 1,
     1-1, OR -2.',62,2,1)
      RETURN
C
   15 IF( LEN(WORK(1)).GE.LEN(X(1)) ) GO TO 16
      CALL XERROR('CSORT- ELEMENTS OF WORK ARRAY TOO SHORT FOR X ',
     1 46,3,1)
      RETURN
   16 IF( (KK.NE.2).OR.(LEN(WORK(1)).GE.LEN(Y(1))) ) GO TO 20
      CALL XERROR('CSORT- ELEMENTS OF WORK ARRAY TOO SHORT FOR Y ',
     1 46,4,1)
      RETURN
   20 IF (KFLAG.GE.1) GO TO 30
   30 GO TO (100,200),KK
C
C SORT X ONLY
C
  100 CONTINUE
      M=1
      I=1
      J=NN
      R=.375
  110 IF (I .EQ. J) GO TO 155
  115 IF (R .GT. .5898437) GO TO 120
      R=R+3.90625E-2
      GO TO 125
  120 R=R-.21875
  125 K=I
C                                  SELECT A CENTRAL ELEMENT OF THE
C                                  ARRAY AND SAVE IT IN LOCATION WORK(1)
      IJ = I + IFIX (FLOAT (J-I) * R)
      WORK(1)=X(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN WORK(1), INTERCHANGE WITH WORK(1
      IF (X(I) .LE. WORK(1)) GO TO 130
      X(IJ)=X(I)
      X(I)=WORK(1)
      WORK(1)=X(IJ)
  130 L=J
C                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
C                                  WORK(1), INTERCHANGE WITH WORK(1)
      IF (X(J) .GE. WORK(1)) GO TO 140
      X(IJ)=X(J)
      X(J)=WORK(1)
      WORK(1)=X(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN WORK(1), INTERCHANGE WITH WORK(1
      IF (X(I) .LE. WORK(1)) GO TO 140
      X(IJ)=X(I)
      X(I)=WORK(1)
      WORK(1)=X(IJ)
      GO TO 140
  135 WORK(2)=X(L)
      X(L)=X(K)
      X(K)=WORK(2)
C                                  FIND AN ELEMENT IN THE SECOND HALF OF
C                                  THE ARRAY WHICH IS SMALLER THAN WORK(
  140 L=L-1
      IF (X(L) .GT. WORK(1)) GO TO 140
C                                  FIND AN ELEMENT IN THE FIRST HALF OF
C                                  THE ARRAY WHICH IS GREATER THAN WORK(
  145 K=K+1
      IF (X(K) .LT. WORK(1)) GO TO 145
C                                  INTERCHANGE THESE ELEMENTS
      IF (K .LE. L) GO TO 135
C                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
C                                  THE ARRAY YET TO BE SORTED
      IF (L-I .LE. J-K) GO TO 150
      IL(M)=I
      IU(M)=L
      I=K
      M=M+1
      GO TO 160
  150 IL(M)=K
      IU(M)=J
      J=L
      M=M+1
      GO TO 160
C                                  BEGIN AGAIN ON ANOTHER PORTION OF
C                                  THE UNSORTED ARRAY
  155 M=M-1
      IF (M .EQ. 0) GO TO 300
      I=IL(M)
      J=IU(M)
  160 IF (J-I .GE. 1) GO TO 125
      IF (I .EQ. 1) GO TO 110
      I=I-1
  165 I=I+1
      IF (I .EQ. J) GO TO 155
      WORK(1)=X(I+1)
      IF (X(I) .LE. WORK(1)) GO TO 165
      K=I
  170 X(K+1)=X(K)
      K=K-1
      IF (WORK(1) .LT. X(K)) GO TO 170
      X(K+1)=WORK(1)
      GO TO 165
C
C SORT X AND CARRY Y ALONG
C
  200 CONTINUE
      M=1
      I=1
      J=NN
      R=.375
  210 IF (I .EQ. J) GO TO 255
  215 IF (R .GT. .5898437) GO TO 220
      R=R+3.90625E-2
      GO TO 225
  220 R=R-.21875
  225 K=I
C                                  SELECT A CENTRAL ELEMENT OF THE
C                                  ARRAY AND SAVE IT IN LOCATION WORK(1)
      IJ = I + IFIX (FLOAT (J-I) *R)
      WORK(1)=X(IJ)
      WORK(3)= Y(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN WORK(1), INTERCHANGE WITH WORK(1
      IF (X(I) .LE. WORK(1)) GO TO 230
      X(IJ)=X(I)
      X(I)=WORK(1)
      WORK(1)=X(IJ)
       Y(IJ)= Y(I)
       Y(I)=WORK(3)
      WORK(3)= Y(IJ)
  230 L=J
C                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
C                                  WORK(1), INTERCHANGE WITH WORK(1)
      IF (X(J) .GE. WORK(1)) GO TO 240
      X(IJ)=X(J)
      X(J)=WORK(1)
      WORK(1)=X(IJ)
       Y(IJ)= Y(J)
       Y(J)=WORK(3)
      WORK(3)= Y(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN WORK(1), INTERCHANGE WITH WORK(1
      IF (X(I) .LE. WORK(1)) GO TO 240
      X(IJ)=X(I)
      X(I)=WORK(1)
      WORK(1)=X(IJ)
       Y(IJ)= Y(I)
       Y(I)=WORK(3)
      WORK(3)= Y(IJ)
      GO TO 240
  235 WORK(2)=X(L)
      X(L)=X(K)
      X(K)=WORK(2)
      WORK(4)= Y(L)
       Y(L)= Y(K)
       Y(K)=WORK(4)
C                                  FIND AN ELEMENT IN THE SECOND HALF OF
C                                  THE ARRAY WHICH IS SMALLER THAN WORK(
  240 L=L-1
      IF (X(L) .GT. WORK(1)) GO TO 240
C                                  FIND AN ELEMENT IN THE FIRST HALF OF
C                                  THE ARRAY WHICH IS GREATER THAN WORK(
  245 K=K+1
      IF (X(K) .LT. WORK(1)) GO TO 245
C                                  INTERCHANGE THESE ELEMENTS
      IF (K .LE. L) GO TO 235
C                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
C                                  THE ARRAY YET TO BE SORTED
      IF (L-I .LE. J-K) GO TO 250
      IL(M)=I
      IU(M)=L
      I=K
      M=M+1
      GO TO 260
  250 IL(M)=K
      IU(M)=J
      J=L
      M=M+1
      GO TO 260
C                                  BEGIN AGAIN ON ANOTHER PORTION OF
C                                  THE UNSORTED ARRAY
  255 M=M-1
      IF (M .EQ. 0) GO TO 300
      I=IL(M)
      J=IU(M)
  260 IF (J-I .GE. 1) GO TO 225
      IF (I .EQ. 1) GO TO 210
      I=I-1
  265 I=I+1
      IF (I .EQ. J) GO TO 255
      WORK(1)=X(I+1)
      WORK(3)= Y(I+1)
      IF (X(I) .LE. WORK(1)) GO TO 265
      K=I
  270 X(K+1)=X(K)
       Y(K+1)= Y(K)
      K=K-1
      IF (WORK(1) .LT. X(K)) GO TO 270
      X(K+1)=WORK(1)
       Y(K+1)=WORK(3)
      GO TO 265
C
C CLEAN UP
C
  300 IF (KFLAG.GE.1) RETURN
C
C          ALTER ARRAY TO GET REVERSE ORDER, IF NECESSARY
C
      NN2=NN/2
      DO 310 I=1,NN2
        IR=NN-I+1
        WORK(1)=X(I)
        X(I)=X(IR)
        X(IR)=WORK(1)
        IF(KFLAG.GE.-1) GO TO 310
        WORK(3)=Y(I)
        Y(I)=Y(IR)
       Y(IR)=WORK(3)
  310 CONTINUE
      RETURN
      END
