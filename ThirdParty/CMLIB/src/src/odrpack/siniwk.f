*SINIWK
      SUBROUTINE SINIWK
     +   (N,M,NP,
     +   X,LDX,IFIXX,LDIFX,SCLD,LDSCLD,LDTT,
     +   BETA,IFIXB,SCLB,
     +   EPSMAC,SSTOL,PARTOL,MAXIT,TAUFAC,
     +   JOB,IPRINT,LUNERR,LUNRPT,
     +   WORK,LWORK,IWORK,LIWORK)
C***BEGIN PROLOGUE  SINIWK
C***REFER TO  SODR,SODRC
C***ROUTINES CALLED  SFLAGS,SIWINF,SSCLB,SSCLD,SWINF,SZERO
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
C***PURPOSE  INITIALIZE WORK VECTORS AS NECESSARY
C***END PROLOGUE  SINIWK
C
C  VARIABLE DECLARATIONS (ALPHABETICALLY)
C
      INTEGER ALPHAI
C        THE LOCATION IN ARRAY WORK OF
C        THE LEVENBERG-MARQUARDT PARAMETER.
      LOGICAL ANAJAC
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE JACOBIANS
C        ARE COMPUTED BY FINITE DIFFERENCES (ANAJAC=.FALSE.) OR NOT
C        (ANAJAC=.TRUE.).
      REAL BETA(NP)
C        THE FUNCTION PARAMETERS.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER BETACI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE CURRENT ESTIMATED VALUES OF THE UNFIXED BETA'S.
      INTEGER BETANI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE NEW ESTIMATED VALUES OF THE UNFIXED BETA'S.
      INTEGER BETASI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE SAVED ESTIMATED VALUES OF THE UNFIXED BETA'S.
      LOGICAL CHKJAC
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE USER-
C        SUPPLIED JACOBIANS ARE TO BE CHECKED (CHKJAC=.TRUE.) OR NOT
C        (CHKJAC=.FALSE.).
      INTEGER DDELTI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE ARRAY (W*D)**2 * DELTA.
      INTEGER DELTAI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE ESTIMATED ERRORS IN THE INDEPENDENT VARIABLES.
      INTEGER DELTNI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE NEW ESTIMATED ERRORS IN THE INDEPENDENT VARIABLES.
      INTEGER DELTSI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE SAVED ESTIMATED ERRORS IN THE INDEPENDENT VARIABLES.
      REAL EPSMAC
C        THE VALUE OF MACHINE PRECISION.
      INTEGER EPSMAI
C        THE LOCATION IN ARRAY WORK OF
C        THE VALUE OF MACHINE PRECISION.
      INTEGER FI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE (WEIGHTED) ESTIMATED VALUES OF EPSILON.
      INTEGER FJACBI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE JACOBIAN WITH RESPECT TO BETA.
      INTEGER FJACXI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE JACOBIAN WITH RESPECT TO X.
      INTEGER FNI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE NEW (WEIGHTED) ESTIMATED VALUES OF EPSILON.
      INTEGER FSI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE SAVED (WEIGHTED) ESTIMATED VALUES OF EPSILON.
      INTEGER I
C        AN INDEXING VARIABLE.
      INTEGER IFIXB(NP)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE INDIVIDUAL
C        ELEMENTS OF BETA ARE FIXED AT THEIR INPUT VALUES OR NOT.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER IFIXX(LDIFX,M)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE INDIVIDUAL
C        ELEMENTS OF DELTA ARE FIXED AT THEIR INPUT VALUES OR NOT.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      LOGICAL INITD
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE DELTA'S
C        ARE TO BE INITIALIZED TO ZERO (INITD=.TRUE.) OR WHETHER THEY
C        ARE TO BE INITIALIZED TO THE VALUES PASSED VIA THE FIRST N BY M
C        ELEMENTS OF ARRAY WORK (INITD=.FALSE.).
      INTEGER INT2I
C        THE LOCATION IN ARRAY IWORK OF
C        THE NUMBER OF INTERNAL DOUBLING STEPS.
      INTEGER IPRINI
C        THE LOCATION IN ARRAY IWORK OF
C        THE PRINT CONTROL VARIABLE.
      INTEGER IPRINT
C        THE PRINT CONTROL VARIABLE.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER IRANKI
C        THE LOCATION IN ARRAY IWORK OF
C        THE RANK DEFICIENCY OF THE JACOBIAN WRT BETA.
      LOGICAL ISODR
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE SOLUTION
C        IS TO BE FOUND BY ODR (ISODR=.TRUE.) OR BY OLS (ISODR=.FALSE.).
      INTEGER IWORK(LIWORK)
C        THE INTEGER WORK SPACE.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER J
C        AN INDEXING VARIABLE.
      INTEGER JOB
C        THE PROBLEM INITIALIZATION AND COMPUTATIONAL
C        METHOD CONTROL VARIABLE.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER JPVTI
C        THE STARTING LOCATION IN ARRAY IWORK OF
C        THE PIVOT VECTOR.
      INTEGER LDIFX
C        THE LEADING DIMENSION OF ARRAY IFIXX.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER LDSCLD
C        THE LEADING DIMENSION OF ARRAY SCLD.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER LDTT
C        THE LEADING DIMENSION OF ARRAY TT.
      INTEGER LDX
C        THE LEADING DIMENSION OF ARRAY X.
      INTEGER LIWORK
C        THE LENGTH OF VECTOR IWORK.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER LUNERI
C        THE LOCATION IN ARRAY IWORK OF
C        THE LOGICAL UNIT NUMBER USED FOR ERROR MESSAGES.
      INTEGER LUNERR
C        THE LOGICAL UNIT NUMBER USED FOR ERROR MESSAGES.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER LUNRPI
C        THE LOCATION IN ARRAY IWORK OF
C        THE LOGICAL UNIT NUMBER USED FOR COMPUTATION REPORTS.
      INTEGER LUNRPT
C        THE LOGICAL UNIT NUMBER USED FOR COMPUTATION REPORTS.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER LWORK
C        THE LENGTH OF VECTOR WORK.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER M
C        THE NUMBER OF COLUMNS OF DATA IN THE INDEPENDENT VARIABLE.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER MAXIT
C        THE MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER MAXITI
C        THE LOCATION IN ARRAY IWORK OF
C        THE MAXIMUM NUMBER OF ITERATIONS ALLOWED.
      INTEGER MSGB
C        THE STARTING LOCATION IN ARRAY IWORK OF
C        THE ERROR CHECKING RESULTS FOR THE JACOBIAN WRT BETA.
      INTEGER MSGX
C        THE STARTING LOCATION IN ARRAY IWORK OF
C        THE ERROR CHECKING RESULTS FOR THE JACOBIAN WRT X.
      INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      REAL NEGONE
C        THE VALUE -1.0E0.
      INTEGER NFEVI
C        THE LOCATION IN ARRAY IWORK OF
C        THE NUMBER OF FUNCTION EVALUATIONS.
      INTEGER NJEVI
C        THE LOCATION IN ARRAY IWORK OF
C        THE NUMBER OF JACOBIAN EVALUATIONS.
      INTEGER NP
C        THE NUMBER OF FUNCTION PARAMETERS.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER OLMAVI
C        THE LOCATION IN ARRAY WORK OF
C        THE AVERAGE NUMBER OF LEVENBERG-MARQUARDT STEPS PER ITERATION.
      INTEGER OMEGAI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE ARRAY (I-FJACX*INV(P)*TRANS(FJACX))**(-1/2)  WHERE
C        P = TRANS(FJACX)*FJACX + D**2 + ALPHA*TT**2
      REAL ONE
C        THE VALUE 1.0E0.
      INTEGER PARTLI
C        THE LOCATION IN ARRAY WORK OF
C        THE PARAMETER CONVERGENCE STOPPING CRITERIA.
      REAL PARTOL
C        THE PARAMETER CONVERGENCE STOPPING CRITERIA.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER QRAUXI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE ARRAY REQUIRED TO RECOVER THE ORTHOGONAL PART OF THE
C        Q-R DECOMPOSITION.
      INTEGER RCONDI
C        THE LOCATION IN ARRAY WORK OF
C        THE APPROXIMATE RECIPROCAL CONDITION OF TFJACB.
      LOGICAL RESTRT
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE CALL IS
C        A RESTART (RESTRT=TRUE) OR NOT (RESTRT=FALSE).
      INTEGER RNORMI
C        THE LOCATION IN ARRAY WORK OF
C        THE NORM OF THE WEIGHTED OBSERVATIONAL ERRORS.
      REAL SCLB(NP)
C        THE SCALE OF EACH BETA.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      REAL SCLD(LDSCLD,M)
C        THE SCALE OF EACH DELTA.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER SI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE STEP FOR THE ESTIMATED BETA'S.
      INTEGER SSFI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE SCALE USED FOR THE BETA'S.
      INTEGER SSI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE SCALE USED FOR THE ESTIMATED BETA'S.
      INTEGER SSSI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE ARRAY USED TO COMPUTED VARIOUS SUMS-OF-SQUARES.
      REAL SSTOL
C        THE SUM-OF-SQUARES CONVERGENCE STOPPING CRITERIA.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER SSTOLI
C        THE LOCATION IN ARRAY WORK OF
C        THE SUM-OF-SQUARES CONVERGENCE STOPPING CRITERIA.
      INTEGER STPI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE STEP USED TO COMPUTE FINITE DIFFERENCE DERIVATIVES.
      REAL TAUFAC
C        THE FACTOR USED TO COMPUTE THE INITIAL TRUST REGION DIAMETER.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER TAUFCI
C        THE LOCATION IN ARRAY WORK OF
C        THE FACTOR USED TO COMPUTE THE INITIAL TRUST REGION DIAMETER.
      INTEGER TAUI
C        THE LOCATION IN ARRAY WORK OF
C        THE TRUST REGION DIAMETER.
      INTEGER TFJACI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE ARRAY INV(DIAG(SQRT(OMEGA(I)),I=1,...,N))*FJACB.
      REAL THREE
C          THE VALUE 3.0E0.
      INTEGER TI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE STEP FOR THE ESTIMATED DELTA'S.
      INTEGER TTI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE SCALE USED FOR THE DELTA'S.
      REAL TWO
C          THE VALUE 2.0E0.
      INTEGER UI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE APPROXIMATE NULL VECTOR FOR TFJACB.
      REAL WORK(LWORK)
C        THE REAL WORK SPACE.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER WRK1I
C        THE STARTING LOCATION IN ARRAY WORK OF
C        A WORK ARRAY.
      INTEGER WRK2I
C        THE STARTING LOCATION IN ARRAY WORK OF
C        A WORK ARRAY.
      REAL X(LDX,M)
C        THE INDEPENDENT VARIABLE.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER XPLUSI
C        THE STARTING LOCATION IN ARRAY WORK OF
C        THE ARRAY X + DELTA.
      INTEGER YTI
C         THE STARTING LOCATION IN WORK OF
C         THE ARRAY -(DIAG(SQRT(OMEGA(I)),I=1,...,N)*(G1-V*INV(E)*D*G2).
      REAL ZERO
C          THE VALUE 0.0E0.
C
C
      DATA NEGONE,ZERO,ONE,TWO,THREE/-1.0E0,0.0E0,1.0E0,2.0E0,3.0E0/
C
C
C***FIRST EXECUTABLE STATEMENT  SINIWK
C
C
      CALL SFLAGS(JOB,RESTRT,INITD,ANAJAC,CHKJAC,ISODR)
C
C  SET STARTING LOCATIONS WITHIN INTEGER WORKSPACE
C
      CALL SIWINF
     +   (M,NP,
     +   MSGB,MSGX,JPVTI,MAXITI,IPRINI,NFEVI,NJEVI,INT2I,IRANKI,
     +   LUNERI,LUNRPI)
C
C  SET STARTING LOCATIONS WITHIN REAL WORK SPACE
C
      CALL SWINF
     +   (N,M,NP,
     +   DELTAI,FI,RNORMI,PARTLI,SSTOLI,TAUFCI,EPSMAI,OLMAVI,
     +   FJACBI,FJACXI,XPLUSI,BETACI,BETASI,BETANI,
     +   DELTSI,DELTNI,DDELTI,FSI,FNI,SI,SSSI,SSI,SSFI,TI,TTI,
     +   TAUI,ALPHAI,TFJACI,OMEGAI,YTI,UI,QRAUXI,WRK1I,WRK2I,STPI,
     +   RCONDI)
C
C  STORE VALUE OF MACHINE PRECISION IN WORK VECTOR
C
      WORK(EPSMAI) = EPSMAC
C
C  SET TOLERANCE FOR STOPPING CRITERIA BASED ON THE CHANGE IN THE
C  PARAMETERS
C
      IF (PARTOL.LT.EPSMAC .OR. PARTOL.GE.ONE) THEN
         WORK(PARTLI) = EPSMAC**(TWO/THREE)
      ELSE
         WORK(PARTLI) = PARTOL
      END IF
C
C  SET TOLERANCE FOR STOPPING CRITERIA BASED ON THE CHANGE IN THE
C  SUM OF SQUARES OF THE WEIGHTED OBSERVATIONAL ERRORS
C
      IF (SSTOL.LT.EPSMAC .OR. SSTOL.GE.ONE) THEN
         WORK(SSTOLI) = SQRT(EPSMAC)
      ELSE
         WORK(SSTOLI) = SSTOL
      END IF
C
C  SET FACTOR FOR COMPUTING TRUST RETION DIAMETER AT FIRST ITERATION
C
      IF (TAUFAC.LE.ZERO) THEN
         WORK(TAUFCI) = ONE
      ELSE
         WORK(TAUFCI) = TAUFAC
      END IF
C
C  SET MAXIMUM NUMBER OF ITERATIONS
C
      IF (MAXIT.LE.0) THEN
         IWORK(MAXITI) = 50
      ELSE
         IWORK(MAXITI) = MAXIT
      END IF
C
C  SET PRINT CONTROL
C
      IF (IPRINT.LT.0) THEN
         IWORK(IPRINI) = 2001
      ELSE
         IWORK(IPRINI) = IPRINT
      END IF
C
C  SET LOGICAL UNIT NUMBER FOR ERROR MESSAGES
C
      IF (LUNERR.LT.0) THEN
         IWORK(LUNERI) = 6
      ELSE
         IWORK(LUNERI) = LUNERR
      END IF
C
C  SET LOGICAL UNIT NUMBER FOR COMPUTATION REPORTS
C
      IF (LUNRPT.LT.0) THEN
         IWORK(LUNRPI) = 6
      ELSE
         IWORK(LUNRPI) = LUNRPT
      END IF
C
C  COMPUTE SCALING FOR BETA'S AND DELTA'S
C
      IF (SCLB(1).LE.ZERO) THEN
         CALL SSCLB(NP,BETA,EPSMAC,WORK(SSFI))
      ELSE
         CALL SCOPY(NP,SCLB,1,WORK(SSFI),1)
      END IF
      IF (SCLD(1,1).LE.ZERO) THEN
         LDTT = N
         CALL SSCLD(N,M,X,LDX,WORK(TTI),LDTT)
      ELSE
         IF (LDSCLD.EQ.1) THEN
            LDTT = 1
            CALL SCOPY(N,SCLD(1,1),1,WORK(TTI),1)
         ELSE
            LDTT = N
            DO 10 J=1,M
               CALL SCOPY(N,SCLD(1,J),1,WORK(TTI+(J-1)*LDTT),1)
   10       CONTINUE
         END IF
      END IF
C
C  INITIALIZE DELTA'S AS NECESSARY
C
      IF (ISODR) THEN
         IF (INITD) THEN
            CALL SZERO(N,M,WORK(DELTAI),N)
         ELSE
            IF (IFIXX(1,1).GE.0) THEN
               IF (LDIFX.EQ.1) THEN
                  DO 20 J=1,M
                     IF (IFIXX(1,J).EQ.0) THEN
                        CALL SZERO(N,1,WORK(DELTAI+(J-1)*N),N)
                     END IF
   20             CONTINUE
               ELSE
                  DO 40 J=1,M
                     DO 30 I=1,N
                        IF (IFIXX(I,J).EQ.0) THEN
                           WORK(DELTAI-1+I+(J-1)*N) = ZERO
                        END IF
   30                CONTINUE
   40             CONTINUE
               END IF
            END IF
         END IF
      ELSE
         CALL SZERO(N,M,WORK(DELTAI),N)
      END IF
C
      RETURN
      END
