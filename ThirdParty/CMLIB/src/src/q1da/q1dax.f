      SUBROUTINE Q1DAX(F,A,B,EPS,R,E,NINT,RST,W,NMAX,FMIN,FMAX,KF,IFLAG)
C***BEGIN PROLOGUE Q1DAX
C***DATE WRITTEN 821118
C***REVISION DATE 821118
C***CATEGORY NO. D1B1A1
C***KEYWORDS  ADAPTIVE QUADRATURE, AUTOMATIC QUADRATURE
C***AUTHOR  KAHANER, DAVID K., SCIENTIFIC COMPUING DIVISION, NBS.
C***PURPOSE  APPROXIMATES ONE DIMENSIONAL INTEGRALS OF USER DEFINED
C            FUNCTIONS, FLEXIBLE USAGE.
C
C***DESCRIPTION
C
C       Q1DAX IS A FLEXIBLE SUBROUTINE FOR THE AUTOMATIC EVALUATION
C             OF DEFINITE INTEGRALS OF A USER DEFINED FUNCTION
C             OF ONE VARIABLE.
C
C             FOR EASIER TO USE ROUTINES SEE Q1DA OR Q1DB.
C
C             CAPABILITIES OF Q1DAX (IN ADDITION TO THOSE OF Q1DA, Q1DB)
C             INCLUDE:
C                ABILITY TO RESTART A CALCULATION TO GREATER
C                  ACCURACY WITHOUT PENALTY...
C                ABILITY TO SPECIFY AN INITIAL PARTITION OF
C                  THE INTEGRATION INTERVAL...
C                ABILITY TO INCREASE THE WORK SPACE TO HANDLE
C                  MORE DIFFICULT PROBLEMS...
C                OUTPUT OF LARGEST/SMALLEST INTEGRAND VALUE FOR
C                  APPLICATIONS SUCH AS SCALING GRAPHS...
C
C       A R G U M E N T S   I N   T H E   C A L L   S E Q U E N C E
C
C       F      (INPUT) THE NAME OF YOUR INTEGRAND FUNCTION.
C                      THIS NAME MUST APPEAR IN AN EXTERNAL STATEMENT
C                      IN ANY PROGRAM WHICH CALLS Q1DAX.
C                      YOU MUST WRITE F IN THE FORM
C                           FUNCTION F(X)
C                              F=(EVALUATE INTEGRAND AT THE POINT X)
C                              RETURN
C                           END
C       A
C       B       (INPUT) ENDPOINTS OF INTEGRATION INTERVAL
C       EPS     (INPUT)  ACCURACY TO WHICH THE INTEGRAL IS TO BE CALCULA
C                          Q1DAX WILL TRY TO ACHIEVE RELATIVE ACCURACY,
C                          E.G. SET EPS=.01 FOR 2 DIGITS, .001 FOR 3, ET
C       R       (OUTPUT) THE ESTIMATE OF THE INTEGRAL
C       E       (OUTPUT) THE ESTIMATE OF THE ABSOLUTE ERROR IN R.
C       NINT    (INPUT
C                OUTPUT)
C                        AS AN INPUT QUANTITY, NINT MUST BE SET TO
C                          THE NUMBER OF SUBINTERVALS IN THE INITIAL
C                          PARTITION OF [A,B].  FOR MOST PROBLEMS
C                          THIS IS JUST 1, THE INTERVAL [A,B] ITSELF.
C                          NINT MUST BE LESS THAN NMAX, SEE BELOW.
C                          NINT IS USEFUL IF YOU WOULD LIKE TO HELP
C                          Q1DAX LOCATE A DIFFICULT SPOT ON [A,B].
C                          IN THIS REGARD NINT IS USED ALONG
C                          WITH THE ARRAY W (SEE BELOW).  IF YOU SET
C                          NINT=1 IT IS NOT NECESSARY TO BE CONCERNED
C                          WITH W, EXCEPT THAT IT MUST BE DIMENSIONED...
C                          AS AN EXAMPLE OF MORE GENERAL APPLICATIONS,
C                          IF [A,B]=[0,1] BUT THE INTEGRAND JUMPS AT 0.3
C                          IT WOULD BE WISE TO SET NINT=2 AND THEN SET
C                                  W(1,1)=0.0  (LEFT ENDPOINT)
C                                  W(2,1)=0.3  (SINGULAR POINT)
C                                  W(3,1)=1.0  (RIGHT ENDPOINT)
C                          IF YOU SET NINT GREATER THAN 1, BE SURE TO
C                          CHECK THAT YOU HAVE ALSO SET
C                                W(1,1)=A  AND  W(NINT+1,1)=B
C                        AS AN OUTPUT QUANTITY, NINT GIVES THE
C                          NUMBER OF SUBINTERVALS IN THE FINAL
C                          PARTITION OF [A,B].
C       RST     (INPUT) A LOGICAL VARIABLE (E.G. TRUE OR FALSE)
C                       SET RST=.FALSE. FOR INITIAL CALL TO Q1DAX
C                       SET RST=.TRUE. FOR A SUBSEQUENT CALL,
C                           E.G. ONE FOR WHICH MORE ACCURACY IS
C                           DESIRED (SMALLER EPS).  A RESTART ONLY
C                           MAKES SENSE IF THE PRECEDING CALL RETURNED
C                           WITH A VALUE OF IFLAG (SEE BELOW) LESS THAN 
C                           ON A RESTART YOU MAY NOT CHANGE THE VALUES O
C                           OTHER ARGUMENTS IN THE CALL SEQUENCE, EXCEPT
C       W(NMAX,6) (INPUT &
C                  SCRATCH)
C                       W IS AN ARRAY USED BY Q1DAX.
C                        YOU  M U S T  INCLUDE A DIMENSION STATEMENT IN
C                        YOUR CALLING PROGRAM TO ALLOCATE THIS STORAGE.
C                        THIS SHOULD BE OF THE FORM
C                                   DIMENSION W(NMAX,6)
C                        WHERE NMAX IS AN INTEGER. AN ADEQUATE VALUE OF
C                        NMAX IS 50.  IF YOU SET NINT>1 YOU MUST ALSO
C                        INITIALIZE W, SEE NINT ABOVE.
C       NMAX    (INPUT) AN INTEGER EQUAL TO THE FIRST SUBSCRIPT IN THE
C                        DIMENSION STATEMENT FOR THE ARRAY W.  THIS IS
C                        ALSO EQUAL TO THE MAXIMUM NUMBER OF SUBINTERVAL
C                        PERMITTED IN THE INTERNAL PARTITION OF [A,B].
C                        A VALUE OF 50 IS AMPLE FOR MOST PROBLEMS.
C       FMIN
C       FMAX    (OUTPUT) THE SMALLEST AND LARGEST VALUES OF THE INTEGRAN
C                          WHICH OCCURRED DURING THE CALCULATION.  THE
C                          ACTUAL INTEGRAND RANGE ON [A,B] MAY, OF COURS
C                          BE GREATER BUT PROBABLY NOT BY MORE THAN 10%.
C       KF      (OUTPUT) THE ACTUAL NUMBER OF INTEGRAND EVALUATIONS USED
C                          BY Q1DAX TO APPROXIMATE THIS INTEGRAL.  KF
C                          WILL ALWAYS BE AT LEAST 30.
C       IFLAG   (OUTPUT) TERMINATION FLAG...POSSIBLE VALUES ARE
C                  0   NORMAL COMPLETION, E SATISFIES
C                           E<EPS  AND  E<EPS*ABS(R)
C                  1   NORMAL COMPLETION, E SATISFIES
C                           E<EPS, BUT E>EPS*ABS(R)
C                  2   NORMAL COMPLETION, E SATISFIES
C                           E<EPS*ABS(R), BUT E>EPS
C                  3   NORMAL COMPLETION BUT EPS WAS TOO SMALL TO
C                        SATISFY ABSOLUTE OR RELATIVE ERROR REQUEST.
C                  4   ABORTED CALCULATION BECAUSE OF SERIOUS ROUNDING
C                        ERROR.  PROBABLY E AND R ARE CONSISTENT.
C                  5   ABORTED CALCULATION BECAUSE OF INSUFFICIENT STORA
C                        R AND E ARE CONSISTENT.  PERHAPS INCREASING NMA
C                        WILL PRODUCE BETTER RESULTS.
C                  6   ABORTED CALCULATION BECAUSE OF SERIOUS DIFFICULTI
C                        MEETING YOUR ERROR REQUEST.
C                  7   ABORTED CALCULATION BECAUSE EITHER EPS, NINT OR N
C                        HAS BEEN SET TO AN ILLEGAL VALUE.
C                  8   ABORTED CALCULATION BECAUSE YOU SET NINT>1 BUT FO
C                        TO SET W(1,1)=A  AND  W(NINT+1,1)=B
C
C     T Y P I C A L   P R O B L E M   S E T   U P
C
C      DIMENSION W(50,6)
C      LOGICAL RST
C      EXTERNAL F
C      A=0.0
C      B=1.0
C      W(1,1)=A
C      W(2,1)=.3      [SET INTERNAL PARTITION POINT AT .3]
C      W(3,1)=B
C      NINT=2         [INITIAL PARTITION HAS 2 INTERVALS]
C      RST=.FALSE.
C      EPS=.001
C      NMAX=50
C
C    1 CALL Q1DAX(F,A,B,EPS,R,E,NINT,RST,W,NMAX,FMIN,FMAX,KF,IFLAG)
C
C      IF(EPS.EQ. .0001 .OR. IFLAG.GE.3)STOP
C      RST=.TRUE
C      EPS=.0001      [ASK FOR ANOTHER DIGIT]
C      GO TO 1
C      END
C      FUNCTION F(X)
C      IF(X.LT. .3)
C     1  THEN
C          F=X**(0.2)*ALOG(X)
C        ELSE
C          F=SIN(X)
C      ENDIF
C      RETURN
C      END
C
C
C            R E M A R K
C
C               WHEN YOU USE Q1ADX WITH NINT=1, WE HAVE BUILT A SMALL
C            AMOUNT OF RANDOMIZATION INTO IT.  REPEATED CALLS DURING
C            THE SAME RUN WILL PRODUCE DIFFERENT, BUT HOPEFULLY
C            CONSISTENT, RESULTS.
C
C    E N D   O F   D O C U M E N T A T I O N
C***REFERENCES  (NONE)
C***ROUTINES CALLED (R1MACH,UNI,GL15T)
C***END PROLOGUE Q1DAX
C
      INTEGER C
      DIMENSION W(NMAX,6)
      EXTERNAL F
      LOGICAL RST
C
C***FIRST EXECUTABLE STATEMENT  Q1DAX
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
      IF(A.EQ.B) THEN
          R=0.
          E=0.
          NINT=0
          IFLAG=0
          KF=1
          FMIN=F(A)
          FMAX=FMIN
          GO TO 20
      ENDIF
      IF(RST) THEN
         IF(IFLAG.LT.3) THEN
           EB=AMAX1(100.*UFLOW,AMAX1(EPS,50.*EPMACH)*ABS(R))
           DO 19 I=1,NINT
               IF(ABS(W(I,3)).GT.(EB*(W(I,2)-W(I,1))/(B-A)))THEN
                                 W(I,3)=ABS(W(I,3))
               ELSE
                                 W(I,3)=-ABS(W(I,3))
               ENDIF
   19      CONTINUE
           GOTO 15
         ELSE
           GOTO 20
         ENDIF
      ENDIF
      KF=0
      IF(EPS .LE. 0. .OR. NINT .LE. 0 .OR. NINT .GE. NMAX) THEN
          IFLAG=7
          GO TO 20
      ENDIF
      IF(NINT.EQ.1)
     1  THEN
          W(1,1)=A
          W(2,2)=B
          W(1,5)=A
          W(1,6)=B
          W(2,5)=A
          W(2,6)=B
C          SELECT FIRST SUBDIVISION RANDOMLY
          W(1,2)=A+(B-A)/2.*(2*UNI(0)+7.)/8.
          W(2,1)=W(1,2)
          NINT=2
        ELSE
          IF(W(1,1).NE.A .OR. W(NINT+1,1).NE.B) THEN
               IFLAG=8
               GO TO 20
          ENDIF
          W(1,5)=A
          DO 89 I=1,NINT
             W(I,2)=W(I+1,1)
             W(I,5)=W(I,1)
             W(I,6)=W(I,2)
   89     CONTINUE
      ENDIF
C
      IFLAG = 0
      IROFF=0
      RABS=0.0
      DO 3 I=1,NINT
          CALL GL15T(F,W(I,1),W(I,2),DBLE(W(I,5)),DBLE(W(I,6)),
     1          W(I,4),W(I,3),RAB,RAV,FMN,FMX)
          KF=KF+15
          IF(I.EQ.1)
     1       THEN
               R=W(I,4)
               E=W(I,3)
               RABS=RABS+RAB
               FMIN=FMN
               FMAX=FMX
             ELSE
               R=R+W(I,4)
               E=E+W(I,3)
               RABS=RABS+RAB
               FMAX=AMAX1(FMAX,FMX)
               FMIN=AMIN1(FMIN,FMN)
          ENDIF
    3 CONTINUE
      DO 10 I=NINT+1,NMAX
          W(I,3) = 0.
   10 CONTINUE
   15 CONTINUE
C
C   MAIN SUBPROGRAM LOOP
C
      IF(100.*EPMACH*RABS.GE.ABS(R) .AND. E.LT.EPS)GO TO 20
      EB=AMAX1(100.*UFLOW,AMAX1(EPS,50.*EPMACH)*ABS(R))
      IF(E.LE.EB) GO TO 20
      IF (NINT.LT.NMAX)
     1 THEN
        NINT = NINT+1
        C = NINT
       ELSE
        C=0
   16   IF(C.EQ.NMAX) THEN
            IFLAG=5
            GO TO 20
        ENDIF
        C=C+1
        IF(W(C,3).GT.0.0) GO TO 16
      END IF
      LOC=ISAMAX(NINT,W(1,3),1)
      XM = W(LOC,1)+(W(LOC,2)-W(LOC,1))/2.
      IF ((AMAX1(ABS(W(LOC,1)),ABS(W(LOC,2)))).GT.
     1   ((1.+100.*EPMACH)*(ABS(XM)+0.1E+04*UFLOW)))
     2    THEN
            CALL GL15T(F,W(LOC,1),XM,DBLE(W(LOC,5)),DBLE(W(LOC,6)),
     1                    TR1,TE1,RAB,RAV,FMINL,FMAXL)
            KF=KF+15
            IF (TE1.LT.(EB*(XM-W(LOC,1))/(B-A))) TE1=-TE1
            CALL GL15T(F,XM,W(LOC,2),DBLE(W(LOC,5)),DBLE(W(LOC,6)),
     1                    TR2,TE2,RAB,RAV,FMINR,FMAXR)
            KF=KF+15
            FMIN=AMIN1(FMIN,FMINL,FMINR)
            FMAX=AMAX1(FMAX,FMAXL,FMAXR)
            IF (TE2.LT.(EB*(W(LOC,2)-XM)/(B-A))) TE2=-TE2
            TE = ABS(W(LOC,3))
            TR = W(LOC,4)
            W(C,3) = TE2
            W(C,4) = TR2
            W(C,1) = XM
            W(C,2) = W(LOC,2)
            W(C,5) = W(LOC,5)
            W(C,6) = W(LOC,6)
            W(LOC,3) = TE1
            W(LOC,4) = TR1
            W(LOC,2) = XM
            E = E-TE+(ABS(TE1)+ABS(TE2))
            R = R-TR+(TR1+TR2)
            IF(ABS(ABS(TE1)+ABS(TE2)-TE).LT.0.001*TE) THEN
                IROFF=IROFF+1
                IF(IROFF.GE.10) THEN
                     IFLAG=4
                     GO TO 20
                ENDIF
            ENDIF
          ELSE
            IF (EB.GT.W(LOC,3))
     1         THEN
                   W(LOC,3) = 0.
               ELSE
                   IFLAG=6
                   GO TO 20
            END IF
      END IF
      GO TO 15
C
C        ALL EXITS FROM HERE
C
   20 CONTINUE
      IF(IFLAG.GE.4)RETURN
      IFLAG=3
      T=EPS*ABS(R)
      IF(E.GT.EPS .AND. E.GT.T)RETURN
      IFLAG=2
      IF(E.GT.EPS .AND. E.LT.T)RETURN
      IFLAG=1
      IF(E.LT.EPS .AND. E.GT.T)RETURN
      IFLAG=0
      RETURN
      END
